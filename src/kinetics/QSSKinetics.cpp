/**
 * @file QSSKinetics.cpp
 * Kinetics manager applying the quasi-steady-state (QSS) approximation.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/QSSKinetics.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

namespace {

//! Total stoichiometric count of QSS species in a composition, and the QSS index
//! of the (last) QSS species found
std::pair<double, size_t> qssContent(const Composition& comp,
                                     const map<string, size_t>& qssIndex)
{
    double count = 0.0;
    size_t idx = npos;
    for (const auto& [name, nu] : comp) {
        auto it = qssIndex.find(name);
        if (it != qssIndex.end()) {
            count += nu;
            idx = it->second;
        }
    }
    return {count, idx};
}

}

void QSSKinetics::setupQSS(const AnyMap& phaseNode, const AnyMap& rootNode)
{
    if (nPhases() != 1) {
        throw InputFileError("QSSKinetics::setupQSS", phaseNode,
            "QSS kinetics requires exactly one phase");
    }
    ThermoPhase& trans = thermo(0);
    const auto& qssNode = phaseNode.at("qss-species");
    m_qssNames.clear();
    if (!(qssNode.is<vector<AnyValue>>() && qssNode.as<vector<AnyValue>>().empty())) {
        m_qssNames = qssNode.asVector<string>();
    }
    m_quadratic = phaseNode.getString("qss-quadratic-coupling", "error");
    if (m_quadratic != "error" && m_quadratic != "neglect") {
        throw InputFileError("QSSKinetics::setupQSS", phaseNode,
            "'qss-quadratic-coupling' must be 'error' or 'neglect', not '{}'",
            m_quadratic);
    }

    // Definitions of the QSS species: explicit list in the phase entry (used
    // when cloning or serializing), else the root 'species' section
    map<string, AnyMap> defs;
    if (phaseNode.hasKey("qss-species-definitions")) {
        for (const auto& item : phaseNode.at("qss-species-definitions").asVector<AnyMap>()) {
            AnyMap def = item;
            def.applyUnits(); // needed when built by Kinetics::clone
            defs[def.at("name").asString()] = def;
        }
    } else if (rootNode.hasKey("species")) {
        for (const auto& item : rootNode.at("species").asVector<AnyMap>()) {
            defs[item.at("name").asString()] = item;
        }
    }

    // Internal full phase: transported species followed by QSS species
    m_full = newThermoModel(trans.type());
    for (size_t k = 0; k < trans.nSpecies(); k++) {
        m_full->addSpecies(trans.species(k));
    }
    map<string, size_t> qssIndex;
    for (size_t i = 0; i < m_qssNames.size(); i++) {
        const string& name = m_qssNames[i];
        if (trans.speciesIndex(name, false) != npos) {
            throw InputFileError("QSSKinetics::setupQSS", phaseNode,
                "QSS species '{}' must not be listed in the phase species", name);
        }
        if (!defs.count(name)) {
            throw InputFileError("QSSKinetics::setupQSS", phaseNode,
                "No definition found for QSS species '{}'", name);
        }
        m_full->addSpecies(newSpecies(defs[name]));
        qssIndex[name] = i;
    }
    m_full->setName(trans.name() + "-qss-full");
    m_full->initThermo();
    m_full->setState_TP(trans.temperature(), trans.pressure());

    // Internal kinetics with all reactions of the mechanism
    AnyMap fullNode = phaseNode;
    fullNode["kinetics"] = "bulk";
    fullNode.erase("qss-species");
    fullNode.erase("qss-species-definitions");
    fullNode.erase("qss-quadratic-coupling");
    m_fullKin = newKinetics("bulk");
    m_fullKin->addThermo(m_full);
    m_fullKin->setParameters(fullNode);
    m_fullKin->init();
    addReactions(*m_fullKin, fullNode, rootNode);

    // Index maps
    m_trans2full.resize(trans.nSpecies());
    for (size_t k = 0; k < trans.nSpecies(); k++) {
        m_trans2full[k] = m_full->speciesIndex(trans.speciesName(k));
    }
    m_qssFull.resize(nQSS());
    for (size_t i = 0; i < nQSS(); i++) {
        m_qssFull[i] = m_full->speciesIndex(m_qssNames[i]);
    }

    // Expose the reactions of the full mechanism through this object
    size_t nRxn = m_fullKin->nReactions();
    m_reactions.clear();
    for (size_t j = 0; j < nRxn; j++) {
        m_reactions.push_back(m_fullKin->reaction(j));
    }
    m_perturb.assign(nRxn, 1.0);
    m_ropf.assign(nRxn, 0.0);
    m_ropr.assign(nRxn, 0.0);
    m_ropnet.assign(nRxn, 0.0);
    resizeReactions();

    // Linear structure of the QSS balance
    m_fwdQSS.assign(nRxn, npos);
    m_revQSS.assign(nRxn, npos);
    m_qssStoich.clear();
    fmt::memory_buffer quadratic;
    for (size_t j = 0; j < nRxn; j++) {
        auto& R = *m_reactions[j];
        for (const auto& [name, order] : R.orders) {
            if (qssIndex.count(name) && order != 1.0) {
                throw CanteraError("QSSKinetics::setupQSS", "Reaction '{}' has a "
                    "non-unity order for QSS species '{}'", R.equation(), name);
            }
        }
        if (R.usesThirdBody() && qssIndex.count(R.thirdBody()->name())) {
            throw CanteraError("QSSKinetics::setupQSS", "Reaction '{}' uses QSS "
                "species '{}' as explicit collider", R.equation(),
                R.thirdBody()->name());
        }
        auto [nf, qf] = qssContent(R.reactants, qssIndex);
        auto [nr, qr] = qssContent(R.products, qssIndex);
        if (nf == 1.0) {
            m_fwdQSS[j] = qf;
        } else if (nf > 1.0) {
            m_fwdQSS[j] = m_skip;
            fmt_append(quadratic, "\n    {} (forward)", R.equation());
        }
        if (R.reversible) {
            if (nr == 1.0) {
                m_revQSS[j] = qr;
            } else if (nr > 1.0) {
                m_revQSS[j] = m_skip;
                fmt_append(quadratic, "\n    {} (reverse)", R.equation());
            }
        }
        for (size_t i = 0; i < nQSS(); i++) {
            double nu = getValue(R.products, m_qssNames[i], 0.0)
                      - getValue(R.reactants, m_qssNames[i], 0.0);
            if (nu != 0.0) {
                m_qssStoich.emplace_back(i, j, nu);
            }
        }
    }
    if (quadratic.size() && m_quadratic == "error") {
        throw CanteraError("QSSKinetics::setupQSS", "Reaction rates quadratic in "
            "QSS concentrations (set 'qss-quadratic-coupling: neglect' to drop "
            "these terms from the QSS balance):{}", to_string(quadratic));
    }

    // QSS terms of each reaction direction, for mass-action scaling
    m_fwdTerms.assign(nRxn, {});
    m_revTerms.assign(nRxn, {});
    for (size_t j = 0; j < nRxn; j++) {
        for (const auto& [name, nu] : m_reactions[j]->reactants) {
            if (qssIndex.count(name)) {
                m_fwdTerms[j].emplace_back(qssIndex[name], nu);
            }
        }
        for (const auto& [name, nu] : m_reactions[j]->products) {
            if (qssIndex.count(name)) {
                m_revTerms[j].emplace_back(qssIndex[name], nu);
            }
        }
    }

    // Stoichiometric coefficients restricted to transported / QSS rows
    vector<size_t> transRow(m_full->nSpecies(), npos), qssRow(m_full->nSpecies(), npos);
    for (size_t k = 0; k < m_trans2full.size(); k++) {
        transRow[m_trans2full[k]] = k;
    }
    for (size_t i = 0; i < nQSS(); i++) {
        qssRow[m_qssFull[i]] = i;
    }
    auto restrict = [&](const Eigen::SparseMatrix<double>& full,
                        Eigen::SparseMatrix<double>& T, Eigen::SparseMatrix<double>& Q) {
        vector<Eigen::Triplet<double>> tt, tq;
        for (int c = 0; c < full.outerSize(); c++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(full, c); it; ++it) {
                size_t k = static_cast<size_t>(it.row());
                if (transRow[k] != npos) {
                    tt.emplace_back(static_cast<int>(transRow[k]), it.col(), it.value());
                } else if (qssRow[k] != npos) {
                    tq.emplace_back(static_cast<int>(qssRow[k]), it.col(), it.value());
                }
            }
        }
        T.resize(m_trans2full.size(), nRxn);
        T.setFromTriplets(tt.begin(), tt.end());
        Q.resize(nQSS(), nRxn);
        Q.setFromTriplets(tq.begin(), tq.end());
    };
    restrict(m_fullKin->reactantStoichCoeffs(), m_reacT, m_reacQ);
    restrict(m_fullKin->productStoichCoeffs(), m_prodT, m_prodQ);

    m_cTrans.resize(trans.nSpecies());
    m_cFull.resize(m_full->nSpecies());
    m_cQSS.resize(nQSS());
    m_work.resize(std::max(m_full->nSpecies(), nRxn));
    m_A.resize(nQSS(), nQSS());
    m_b.resize(nQSS());
}

bool QSSKinetics::addReaction(shared_ptr<Reaction> r, bool resize)
{
    throw NotImplementedError("QSSKinetics::addReaction",
        "Reactions of a QSS kinetics manager are defined by its input file");
}

pair<size_t, size_t> QSSKinetics::checkDuplicates(bool throw_err, bool fix)
{
    return m_fullKin->checkDuplicates(throw_err, fix);
}

void QSSKinetics::setMultiplier(size_t i, double f)
{
    Kinetics::setMultiplier(i, f);
    m_fullKin->setMultiplier(i, f);
    m_cacheT = -1.0;
}

void QSSKinetics::updateQSS()
{
    ThermoPhase& trans = thermo(0);
    double T = trans.temperature();
    trans.getConcentrations(m_cTrans.data());
    if (T == m_cacheT && m_cTrans == m_cacheC) {
        return;
    }
    double ctot = 0.0;
    for (size_t k = 0; k < m_cTrans.size(); k++) {
        m_cFull[m_trans2full[k]] = m_cTrans[k];
        ctot += m_cTrans[k];
    }
    size_t nq = nQSS();

    // Probe evaluation: every QSS species at a small concentration eps. A rate
    // linear in one QSS concentration C_q is then (probe rate) * C_q / eps.
    double eps = m_epsRel * std::max(ctot, SmallNumber);
    for (size_t i = 0; i < nq; i++) {
        m_cFull[m_qssFull[i]] = eps;
    }
    m_full->setTemperature(T);
    m_full->setConcentrationsNoNorm(m_cFull.data());
    m_fullKin->getFwdRatesOfProgress(m_ropf.data());
    m_fullKin->getRevRatesOfProgress(m_ropr.data());

    std::fill(m_cQSS.begin(), m_cQSS.end(), 0.0);
    if (nq) {
        m_A.setZero();
        m_b.setZero();
        for (const auto& [i, j, nu] : m_qssStoich) {
            size_t qf = m_fwdQSS[j];
            if (qf == npos) {
                m_b(i) += nu * m_ropf[j];
            } else if (qf != m_skip) {
                m_A(i, qf) += nu * m_ropf[j] / eps;
            }
            size_t qr = m_revQSS[j];
            if (qr == npos) {
                m_b(i) -= nu * m_ropr[j];
            } else if (qr != m_skip) {
                m_A(i, qr) -= nu * m_ropr[j] / eps;
            }
        }

        // Solve A C = -b on the QSS species that are actually consumed
        vector<size_t> active;
        for (size_t i = 0; i < nq; i++) {
            if (m_A(i, i) < 0.0) {
                active.push_back(i);
            }
        }
        if (!active.empty()) {
            size_t na = active.size();
            // Rows scaled by the consumption coefficient (unit diagonal): the
            // entries otherwise span many orders of magnitude
            Eigen::MatrixXd Aa(na, na);
            Eigen::VectorXd ba(na);
            for (size_t r = 0; r < na; r++) {
                double scale = -1.0 / m_A(active[r], active[r]);
                ba(r) = -m_b(active[r]) * scale;
                for (size_t c = 0; c < na; c++) {
                    Aa(r, c) = m_A(active[r], active[c]) * scale;
                }
            }
            auto lu = Aa.partialPivLu();
            Eigen::VectorXd x = lu.solve(ba);
            x += lu.solve(ba - Aa * x); // one step of iterative refinement
            // Not clipped at zero: tiny negative values (round-off near
            // equilibrium) are harmless, while clipping makes the rates
            // non-smooth and stalls stiff integrators
            for (size_t r = 0; r < na; r++) {
                m_cQSS[active[r]] = std::isfinite(x(r)) ? x(r) : 0.0;
            }
        }

        // Rates of progress at the QSS state, by mass action from the probe
        for (size_t j = 0; j < m_ropf.size(); j++) {
            for (const auto& [i, nu] : m_fwdTerms[j]) {
                m_ropf[j] *= std::pow(m_cQSS[i] / eps, nu);
            }
            for (const auto& [i, nu] : m_revTerms[j]) {
                m_ropr[j] *= std::pow(m_cQSS[i] / eps, nu);
            }
        }
    }
    for (size_t j = 0; j < m_ropf.size(); j++) {
        m_ropnet[j] = m_ropf[j] - m_ropr[j];
    }
    m_cacheT = T;
    m_cacheC = m_cTrans;
}

void QSSKinetics::getQSSConcentrations(double* conc)
{
    updateQSS();
    std::copy(m_cQSS.begin(), m_cQSS.end(), conc);
}

void QSSKinetics::getQSSResiduals(double* wdot)
{
    updateQSS();
    Eigen::Map<const Eigen::VectorXd> rop(m_ropnet.data(), m_ropnet.size());
    Eigen::Map<Eigen::VectorXd> out(wdot, nQSS());
    out = m_prodQ * rop - m_reacQ * rop;
}

void QSSKinetics::getNetProductionRates(double* wdot)
{
    updateQSS();
    Eigen::Map<const Eigen::VectorXd> rop(m_ropnet.data(), m_ropnet.size());
    Eigen::Map<Eigen::VectorXd> out(wdot, m_kk);
    out = m_prodT * rop - m_reacT * rop;
}

void QSSKinetics::getCreationRates(double* cdot)
{
    updateQSS();
    Eigen::Map<const Eigen::VectorXd> ropf(m_ropf.data(), m_ropf.size());
    Eigen::Map<const Eigen::VectorXd> ropr(m_ropr.data(), m_ropr.size());
    Eigen::Map<Eigen::VectorXd> out(cdot, m_kk);
    out = m_prodT * ropf + m_reacT * ropr;
}

void QSSKinetics::getDestructionRates(double* ddot)
{
    updateQSS();
    Eigen::Map<const Eigen::VectorXd> ropf(m_ropf.data(), m_ropf.size());
    Eigen::Map<const Eigen::VectorXd> ropr(m_ropr.data(), m_ropr.size());
    Eigen::Map<Eigen::VectorXd> out(ddot, m_kk);
    out = m_reacT * ropf + m_prodT * ropr;
}

void QSSKinetics::getFwdRatesOfProgress(double* fwdROP)
{
    updateQSS();
    std::copy(m_ropf.begin(), m_ropf.end(), fwdROP);
}

void QSSKinetics::getRevRatesOfProgress(double* revROP)
{
    updateQSS();
    std::copy(m_ropr.begin(), m_ropr.end(), revROP);
}

void QSSKinetics::getNetRatesOfProgress(double* netROP)
{
    updateQSS();
    std::copy(m_ropnet.begin(), m_ropnet.end(), netROP);
}

void QSSKinetics::getFwdRateConstants(double* kfwd)
{
    updateQSS();
    m_fullKin->getFwdRateConstants(kfwd);
}

void QSSKinetics::getRevRateConstants(double* krev, bool doIrreversible)
{
    updateQSS();
    m_fullKin->getRevRateConstants(krev, doIrreversible);
}

void QSSKinetics::getEquilibriumConstants(double* kc)
{
    updateQSS();
    m_fullKin->getEquilibriumConstants(kc);
}

void QSSKinetics::getParameters(AnyMap& phaseNode) const
{
    phaseNode["qss-species"] = m_qssNames;
    vector<AnyMap> defs;
    for (const auto& name : m_qssNames) {
        defs.push_back(m_full->species(name)->parameters(m_full.get()));
    }
    phaseNode["qss-species-definitions"] = std::move(defs);
    if (m_quadratic != "error") {
        phaseNode["qss-quadratic-coupling"] = m_quadratic;
    }
    AnyMap fullParams = m_fullKin->parameters();
    if (fullParams.hasKey("skip-undeclared-third-bodies")) {
        phaseNode["skip-undeclared-third-bodies"] = true;
    }
}

}

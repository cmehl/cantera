/**
 * @file QSSKinetics.h
 * Kinetics manager applying the quasi-steady-state (QSS) approximation.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_QSSKINETICS_H
#define CT_QSSKINETICS_H

#include "Kinetics.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera
{

//! Kinetics manager for a gas phase where some species are in quasi-steady state.
/*!
 * The phase attached to this kinetics manager contains only the *transported*
 * species. The QSS species are not part of the phase: their concentrations are
 * obtained at every evaluation from the algebraic balance
 * @f$ \dot\omega_i = 0 @f$ for each QSS species @f$ i @f$, given the
 * temperature and the concentrations of the transported species.
 *
 * Internally, a private ideal-gas phase containing transported + QSS species
 * and a BulkKinetics object holding all reactions are used to evaluate rates.
 * Each rate of every reaction is assumed to be at most linear in the QSS
 * concentrations (no two QSS molecules on one side of a reaction direction).
 * The QSS production rates are then affine in the QSS concentrations,
 * @f$ \dot\omega_{QSS} = A C_{QSS} + b @f$. A single evaluation of the rates
 * of progress with all QSS concentrations at a small probe value
 * @f$ \epsilon @f$ gives A and b; once the linear system is solved, the rates
 * of progress at the QSS state follow by mass action, scaling each probe rate
 * by @f$ \prod_i (C_i / \epsilon)^{\nu_i} @f$ over its QSS species, and
 * production rates follow from the stoichiometric matrices. The dependence of
 * third-body and falloff rates on the (small) QSS concentrations is neglected.
 *
 * YAML phase definition:
 * @code{.yaml}
 * - name: gas
 *   thermo: ideal-gas
 *   species: [H2, O2, ...]        # transported species only
 *   kinetics: qss
 *   qss-species: [CH2(S), HCO]    # defined in the 'species' section, or in
 *                                 # 'qss-species-definitions'
 *   qss-quadratic-coupling: error # or 'neglect'
 *   transport: mixture-averaged
 * @endcode
 *
 * Reaction-level quantities (rates of progress, rate constants) refer to the
 * reactions of the full mechanism, evaluated at the QSS state. Composition
 * derivatives (used by preconditioned or "mole" reactor Jacobians) are not
 * available.
 *
 * @ingroup kineticsmgr
 */
class QSSKinetics : public Kinetics
{
public:
    QSSKinetics() = default;

    string kineticsType() const override {
        return "qss";
    }

    //! Build the internal full mechanism. Called by newKinetics() after the
    //! transported phase has been added.
    void setupQSS(const AnyMap& phaseNode, const AnyMap& rootNode);

    bool addReaction(shared_ptr<Reaction> r, bool resize=true) override;
    //! Duplicate detection needs the QSS species: delegated to the full mechanism
    pair<size_t, size_t> checkDuplicates(bool throw_err=true, bool fix=false) override;
    void setMultiplier(size_t i, double f) override;

    //! @name Rates evaluated at the QSS state
    //! @{
    void getFwdRatesOfProgress(double* fwdROP) override;
    void getRevRatesOfProgress(double* revROP) override;
    void getNetRatesOfProgress(double* netROP) override;
    void getCreationRates(double* cdot) override;
    void getDestructionRates(double* ddot) override;
    void getNetProductionRates(double* wdot) override;
    void getFwdRateConstants(double* kfwd) override;
    void getRevRateConstants(double* krev, bool doIrreversible=false) override;
    void getEquilibriumConstants(double* kc) override;
    //! @}

    //! Names of the QSS species
    const vector<string>& qssSpeciesNames() const {
        return m_qssNames;
    }

    //! Number of QSS species
    size_t nQSS() const {
        return m_qssNames.size();
    }

    //! QSS concentrations [kmol/m^3] at the current state of the phase
    void getQSSConcentrations(double* conc);

    //! Net production rates of the QSS species at the current state
    //! (residual of the QSS balance, ideally zero) [kmol/m^3/s]
    void getQSSResiduals(double* wdot);

    //! Internal full phase (transported + QSS species)
    shared_ptr<ThermoPhase> fullPhase() const {
        return m_full;
    }

    void getParameters(AnyMap& phaseNode) const override;

protected:
    //! Solve the QSS concentrations and the rates of progress at the QSS state
    //! for the current state of the transported phase (cached)
    void updateQSS();

    shared_ptr<ThermoPhase> m_full; //!< transported + QSS species
    shared_ptr<Kinetics> m_fullKin; //!< all reactions, on m_full

    vector<string> m_qssNames;
    vector<size_t> m_trans2full; //!< index in m_full of each transported species
    vector<size_t> m_qssFull; //!< index in m_full of each QSS species

    //! For each reaction: index (in the QSS list) of the QSS species appearing
    //! linearly in the forward / reverse rate, npos if none, or #m_skip if the
    //! direction is quadratic in QSS concentrations and neglected
    vector<size_t> m_fwdQSS, m_revQSS;
    //! Non-zero net stoichiometric coefficients of QSS species:
    //! (QSS index, reaction index, coefficient)
    vector<std::tuple<size_t, size_t, double>> m_qssStoich;
    //! QSS species (QSS index, stoichiometric coefficient) on the reactant /
    //! product side of each reaction, used to scale probe rates
    vector<vector<std::pair<size_t, double>>> m_fwdTerms, m_revTerms;

    //! Reactant / product stoichiometric coefficients restricted to the
    //! transported species (rows) and to the QSS species (rows)
    Eigen::SparseMatrix<double> m_reacT, m_prodT, m_reacQ, m_prodQ;

    string m_quadratic = "error";
    double m_epsRel = 1e-10; //!< probe QSS concentration relative to total

    vector<double> m_cTrans, m_cFull, m_cQSS, m_work;
    vector<double> m_cacheC; //!< transported concentrations of the cached state
    double m_cacheT = -1.0; //!< temperature of the cached state
    Eigen::MatrixXd m_A;
    Eigen::VectorXd m_b;

    static const size_t m_skip = npos - 1;
};

}

#endif

"""Tests for the quasi-steady-state kinetics manager (``kinetics: qss``)."""

import io
from pathlib import Path

import numpy as np
import pytest
from pytest import approx
from ruamel.yaml import YAML

import cantera as ct


def _h2o2_path():
    for d in ct.get_data_directories():
        p = Path(d) / "h2o2.yaml"
        if p.exists():
            return p
    raise FileNotFoundError("h2o2.yaml")


def qss_yaml(qss, **extra):
    """h2o2.yaml with the first phase turned into QSS kinetics."""
    yml = YAML()
    data = yml.load(_h2o2_path().read_text())
    phase = data["phases"][0]
    data["phases"] = [phase]
    phase["species"] = [s for s in phase["species"] if s not in qss]
    phase["kinetics"] = "qss"
    phase["qss-species"] = list(qss)
    phase.update(extra)
    buf = io.StringIO()
    yml.dump(data, buf)
    return buf.getvalue()


STATE = (1200.0, 5 * ct.one_atm, "H2:2, O2:1, N2:3.76, H:1e-3, OH:1e-3, O:1e-4, HO2:1e-5")


class TestQSSKinetics:

    @pytest.fixture(scope="class")
    def bulk(self):
        return ct.Solution("h2o2.yaml", "ohmech")

    def test_load(self):
        gas = ct.Solution(yaml=qss_yaml(["H2O2"]))
        assert gas.kinetics_model == "qss"
        assert gas.qss_species_names == ["H2O2"]
        assert "H2O2" not in gas.species_names
        assert gas.n_species == 9
        assert gas.n_reactions == 29

    def test_no_qss_species_matches_bulk(self, bulk):
        gas = ct.Solution(yaml=qss_yaml([]))
        gas.TPX = STATE
        bulk.TPX = STATE
        assert gas.net_production_rates == approx(bulk.net_production_rates,
                                                  rel=1e-12, abs=1e-20)
        assert gas.net_rates_of_progress == approx(bulk.net_rates_of_progress,
                                                   rel=1e-12, abs=1e-20)

    def test_qss_balance(self, bulk):
        gas = ct.Solution(yaml=qss_yaml(["H2O2"]))
        gas.TPX = STATE
        c = gas.qss_concentrations
        assert c[0] > 0
        # QSS species balance is satisfied
        scale = np.abs(gas.destruction_rates).max()
        assert abs(gas.qss_residuals[0]) < 1e-10 * scale

        # Same rates as the full mechanism with H2O2 at its QSS concentration
        # (up to the third-body approximation: QSS species are not counted in
        # the collider concentration)
        bulk.TP = gas.T, gas.P
        conc = np.zeros(bulk.n_species)
        for k, name in enumerate(gas.species_names):
            conc[bulk.species_index(name)] = gas.concentrations[k]
        conc[bulk.species_index("H2O2")] = c[0]
        bulk.concentrations = conc
        ref = bulk.net_production_rates
        idx = [bulk.species_index(name) for name in gas.species_names]
        assert gas.net_production_rates == approx(ref[idx], rel=1e-4,
                                                  abs=1e-9 * np.abs(ref).max())

    def test_quadratic_coupling(self):
        with pytest.raises(ct.CanteraError, match="quadratic"):
            ct.Solution(yaml=qss_yaml(["HO2"]))
        gas = ct.Solution(yaml=qss_yaml(["HO2"], **{"qss-quadratic-coupling": "neglect"}))
        gas.TPX = STATE[0], STATE[1], "H2:2, O2:1, N2:3.76, H:1e-3, OH:1e-3, O:1e-4"
        assert gas.qss_concentrations[0] > 0

    def test_invalid_species(self):
        with pytest.raises(ct.CanteraError, match="No definition"):
            ct.Solution(yaml=qss_yaml(["H2O2"]).replace("- H2O2", "- XYZ"))

    def test_bulk_has_no_qss(self, bulk):
        assert bulk.qss_species_names == []
        with pytest.raises(ValueError):
            bulk.qss_concentrations

    def test_serialization(self, tmp_path):
        gas = ct.Solution(yaml=qss_yaml(["H2O2"]))
        gas.TPX = STATE
        wdot = gas.net_production_rates
        fname = tmp_path / "qss.yaml"
        gas.write_yaml(fname)
        gas2 = ct.Solution(fname)
        gas2.TPX = STATE
        assert gas2.kinetics_model == "qss"
        assert gas2.qss_species_names == ["H2O2"]
        assert gas2.net_production_rates == approx(wdot, rel=1e-12, abs=1e-20)

    def test_reactor_clone(self):
        gas = ct.Solution(yaml=qss_yaml(["H2O2"]))
        temps = []
        for clone in (True, False):
            gas.TPX = STATE
            r = ct.IdealGasConstPressureReactor(gas, clone=clone)
            net = ct.ReactorNet([r])
            net.advance(1e-3)
            temps.append(r.T)
        assert temps[0] == approx(temps[1], rel=1e-8)
        assert temps[0] > 2000

    def test_ignition_close_to_full(self, bulk):
        def ignition(gas):
            gas.TPX = 1100, ct.one_atm, "H2:2, O2:1, N2:3.76"
            r = ct.IdealGasConstPressureReactor(gas, clone=False)
            net = ct.ReactorNet([r])
            t, T = [0.0], [r.T]
            while net.time < 2e-3:
                net.step()
                t.append(net.time)
                T.append(r.T)
            return t[int(np.argmax(np.gradient(T, t)))]

        gas = ct.Solution(yaml=qss_yaml(["H2O2"]))
        assert ignition(gas) == approx(ignition(bulk), rel=0.05)

    def test_free_flame(self, bulk):
        def flame_speed(gas):
            gas.TPX = 300, ct.one_atm, "H2:1.5, O2:1, N2:3.76"
            f = ct.FreeFlame(gas, width=0.02)
            f.set_refine_criteria(ratio=3, slope=0.1, curve=0.2)
            f.solve(loglevel=0, auto=True)
            return f.velocity[0]

        gas = ct.Solution(yaml=qss_yaml(["H2O2"]))
        assert flame_speed(gas) == approx(flame_speed(bulk), rel=0.02)

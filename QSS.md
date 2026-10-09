# QSS kinetics in Cantera (mechred fork)

> **This is a modified version of Cantera 3.2.0, not an official Cantera
> release.** The QSS kinetics model described below was added by IFP Energies
> nouvelles (IFPEN) for mechred, a reduction tool for
> chemical kinetic mechanisms. Cantera is distributed under the BSD 3-clause
> license (`License.txt`); the additions of this branch are under the same
> license. Copyright of the additions: 2026 IFP Energies nouvelles (IFPEN).
> They were developed with the assistance of an AI coding assistant (Claude,
> Anthropic) and reviewed by the author. Everything outside the files listed
> in *Changes from Cantera 3.2.0* below is unchanged Cantera 3.2.0.

This branch (`qss`, based on Cantera 3.2.0) adds a kinetics model, `qss`, which
applies the quasi-steady-state (QSS) approximation to selected species. A QSS
species is not transported: its concentration is solved from ω̇ = 0 at each
rate evaluation. Every Cantera solver that uses net production rates works
unchanged with it, including 0D reactors and the free, counterflow and
burner flame solvers.

## Mechanism file

Take a normal mechanism and change three things in the phase entry:

1. Remove the QSS species from the phase `species` list.
2. Set `kinetics: qss`.
3. List the QSS species under `qss-species`.

```yaml
phases:
- name: gas
  thermo: ideal-gas
  elements: [H, O, C, N]
  species: [H2, H, O, O2, OH, H2O, HO2, H2O2, CH2, CH3, CH4, CO, CO2, ...]
  kinetics: qss
  qss-species: [CH3O, CH2(S), C2H5, HCO, CH2CHO, C2H3, CH]
  qss-quadratic-coupling: error   # optional: error (default) | neglect
  transport: mixture-averaged
species:
- name: H2
  ...
- name: CH3O        # QSS species keep their definitions in 'species'
  ...
reactions:
- ...               # unchanged: all reactions, including QSS species
```

- QSS species definitions are read from the root `species` section, or from a
  `qss-species-definitions` list in the phase entry. Cantera writes that list
  itself when a QSS Solution is saved with `write_yaml`.
- `qss-quadratic-coupling`: a reaction direction with two QSS molecules on its
  reactant side, such as `2 CH2 <=> ...` or `CH2O + HO2 <=> H2O2 + HCO` with
  H2O2 and HCO both QSS, makes the QSS balance nonlinear.
  - `error` (default) refuses such a QSS set.
  - `neglect` drops those terms from the QSS balance. They are still used
    when computing the transported species' rates.
- Not supported: QSS species with non-unity reaction orders, and QSS species
  used as explicit colliders such as `(+CH2)`.

## Python API

```python
gas = ct.Solution("mech_qss.yaml")
gas.kinetics_model        # 'qss'
gas.species_names         # transported species only
gas.qss_species_names     # ['CH3O', 'CH2(S)', ...]
gas.qss_concentrations    # QSS concentrations at the current state [kmol/m^3]
gas.qss_residuals         # QSS net production rates (≈ 0) [kmol/m^3/s]
gas.net_production_rates  # transported species, at the QSS state
gas.n_reactions           # all reactions of the mechanism
gas.net_rates_of_progress # per reaction, at the QSS state
```

Bulk kinetics returns `qss_species_names == []`.

## How it works

`QSSKinetics` ([include/cantera/kinetics/QSSKinetics.h](include/cantera/kinetics/QSSKinetics.h),
[src/kinetics/QSSKinetics.cpp](src/kinetics/QSSKinetics.cpp)) is attached to
the phase of transported species. Internally, it holds an ideal-gas phase with
the transported and QSS species, and a `BulkKinetics` object with all the
reactions. Each evaluation, cached until the state changes, goes like this:

1. **Probe:** compute the rates of progress with every QSS concentration set to
   a small probe value ε (10⁻¹⁰ of the total concentration). Each rate is
   linear in at most one QSS concentration, so this single evaluation gives
   the QSS balance in linear form, ω̇_QSS = A·C_QSS + b.
2. **Solve:** solve A·C = −b by LU decomposition, with rows scaled to a unit
   diagonal and one step of iterative refinement. QSS species that nothing
   consumes at this state get zero. Results are not clipped at zero: clipping
   makes the rates non-smooth, and stiff integrators stall on it.
3. **Rates:** compute the rates of progress at the QSS state by mass action.
   Each probe rate is scaled by ∏(C_i/ε)^ν_i over its QSS species. Production
   rates of the transported species then follow from the stoichiometric
   matrices.

The QSS species are not counted in third-body and falloff collider
concentrations. QSS mole fractions are small (about 10⁻⁴ at most), so the
effect is of order 10⁻⁴ on those rates.

## Limitations

- Composition derivatives (`*_ddX`, `*_ddC`) are not available. They are only
  needed by "mole" reactors with preconditioners, so use the standard reactors
  (`IdealGasReactor`, `IdealGasConstPressureReactor`). Finite-difference
  Jacobians, which flames and the standard reactors use by default, work.
- Only a single ideal-gas phase is supported.
- Reactions cannot be added or modified at run time (`add_reaction` raises).

## Changes from Cantera 3.2.0

| file | change |
|---|---|
| `include/cantera/kinetics/QSSKinetics.h`, `src/kinetics/QSSKinetics.cpp` | new kinetics manager |
| `src/kinetics/KineticsFactory.cpp` | registers `qss`; `newKinetics` hands reactions to the QSS manager |
| `include/cantera/kinetics/Kinetics.h`, `src/kinetics/Kinetics.cpp` | virtual `getParameters()` hook called by `parameters()`, so that `clone` and `write_yaml` keep the QSS settings |
| `interfaces/cython/cantera/kinetics.pxd`, `kinetics.pyx` | Python properties `qss_species_names`, `qss_concentrations`, `qss_residuals` |
| `test/python/test_qss.py` | tests |

## Building (IFPEN workstation)

The Python 3.12 used here was built with GCC 13.3, so Cantera must be built
with the same compiler. With the system GCC 11, linking the Python-extension
shim fails on link-time-optimised objects.

`cantera.conf` uses the bundled dependencies:

```python
python_cmd = '/work/mehlc/0_CODES_GITS/MechanismReduction/venv/bin/python'
prefix = '/work/mehlc/0_CODES_GITS/cantera-qss/install'
python_package = 'y'
f90_interface = 'n'
system_sundials = 'n'
system_yamlcpp = 'n'
system_fmt = 'n'
system_eigen = 'n'
system_blas_lapack = 'n'
hdf_support = 'n'
optimize = True
debug = False
```

The build environment also needs `scons cython ruamel.yaml jinja2 pytest`.

```bash
module load GCC/13.3.0
scons build -j16            # or: bash -l ./build_qss.sh
scons install               # installs libraries into install/
pip install --force-reinstall --no-deps build/python/dist/cantera-3.2.0-*.whl
```

The compiled module finds `libcantera_shared` through its rpath
(`install/lib`). Running `ldd` from a shell without the GCC 13 module reports a
missing `GLIBCXX_3.4.32`; the Python interpreter loads the newer `libstdc++`
itself, so imports work.

Tests:

```bash
cd test/python
pytest test_qss.py test_kinetics.py test_composite.py test_reactor.py
```

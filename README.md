openmm-validation
=================

An extensive validation suite for the OpenMM molecular simulation library.

**THIS PROJECT IS UNDER DEVELOPMENT AND NOT READY FOR USE**

## Manifest

* `src/` - code for various test systems
* `tests/` - nose-automated validation tests

## Testing philosophy

Tests for OpenMM are broken down into several categories:
* **Unit tests** are integrated into the main [OpenMM](http://github.com/simtk/openmm) source tree, and test individual components (such as `Force`s an `Integrator`s) separately to ensure all features function as intended.
* **Integration tests** ensure that all components of OpenMM function together in the correct way, and are themselves split into several categories:
  * **Statistical mechanics tests** ensure that the statistical mechanical properties of simple systems are reproduced as expeted.
  * **Stability tests** ensure that simulations of typical systems of interest are stable and robust.
  * **Numerical tests** ensure that the various components of OpenMM are numerically well-behaved, and agree between `Platform`s.
* **Comparison tests** ensure that implementations of forcefields and force terms from other software packages (such as [AMBER](http://ambermd.org), [CHARMM](http://www.charmm.org), and [Tinker](http://dasher.wustl.edu/tinker/)) agree to acceptable error tolerance.

## Unresolved questions

### Stochastic tests
Some tests compute the expectation of a quantity from a simulation of finite length.
* What criteria should we use to judge if the test has passed or failed?

## Detailed test description

### Unit tests

Unit tests are integrated into the main [OpenMM](http://github.com/simtk/openmm) source tree, and test individual components (such as `Force`s an `Integrator`s) separately to ensure all features function as intended.

The following is a partial list of unit tests currently implemented in OpenMM. New unit tests are constantly being added.

#### Andersen thermostat
* **testTemperature.** Ensure average of instantaneous kinetic temperature for 8 charged LJ particles matches control temperature when an Andersen thermostat is used over 150 ps simulation.
* **testConstraints.** Ensure average of instantaneous kinetic temperature for 8 charged LJ particles containing several bond constraints matches control temperature when an Andersen thermostat is used over 150 ps simulation.
* **testRandomSeed.** Ensure that setting the random seed to the same value produces identical particle positions following integration from the same initial positions and momenta, while changing the random seed produces a different set of particle positions.

#### Brownian integrator
* **testSingleBond.** Compare the time-dependent positions and velocities of a single harmonic bond to the analytical solution for a damped harmonic oscillator.
* **testTemperature.** Ensure average of potential energy for a chain of 8 harmonically bonded particles matches matches analytically expected result. **NOTE: This test differs from the Andersen thermostat test of the same name.**
* **testConstraints.** Simulate a chain of 8 charged LJ particles and ensure that constraints remain satisfied. **NOTE: This test differs from the Andersen thermostat test of the same name.**
* **testConstrainedMasslessParticles.** Ensure that attempting to constrain a massless (fixed) particle throws an exception, while making both particles massless (fixed) does not throw an exception despite the presence of a constraint.
* **testRandomSeed.** Ensure that setting the random seed to the same value produces identical particle positions following integration from the same initial positions and momenta, while changing the random seed produces a different set of particle positions.

#### CMAPTorsionForce
* **testCMAPTorsions.** Compare a 5-particle system with two torsions to ensure that a traditional `PeriodicTorsionForce` and the corresponding `CMAPTorsionForce` give the same energies (rel tol 0.001) and forces (rel tol 0.05) for 50 random configurations.

#### CMMotionRemover
* **testMotionRemoval.** Ensure that the center of mass of a system of 8 charged LJ particles (also containing a harmonic bond) is preserved when `CMMotionRemover` is in use. Both `LangevinIntegrator` and `VerletIntegrator` are tested.

#### Checkpoints
* `testCheckpoint.** Ensure that trajectories generated after saving a checkpoint and after resuming from the checkpoint generate identical positions, velocities, and parameters.
* `testSetState.** Perform the same test as `testCheckpoint` but utilize the `Context::set{State|Parameter|PeriodicBoxVectors|Positions}` interface instead.

### Integration tests

#### Statistical mechanics tests
Statistical mechanics tests ensure that the statistical mechanical properties of simple systems are reproduced as expeted.

* CheckEnsemble from Michael Shirts: [code](https://github.com/shirtsgroup/checkensemble) [DOI](http://pubs.acs.org/doi/abs/10.1021/ct300688p)

##### Standard battery of test systems
We have implemented a battery of test systems that span a variety of complexities:
https://github.com/choderalab/repex/blob/master/repex/testsystems.py

This includes:
* `HarmonicOscillator`: A 3D harmonic oscillator.
* `PowerOscillator`: A single particle confined to an x^d well.
* `Diatom`: A free diatomic molecule in a harmonic well, with or without bond constraint.
* `DiatomicFluid`: A periodic box of diatomic particles, with or without bond constraint.
* `DipolarFluid`: A diatomic fluid with charges to create dipolar particles.
* `ConstraintCoupledHarmonicOscillator`: A pair of particles in 3D harmonic oscillator wells, coupled by a constraint
* `HarmonicOscillatorArray`: A 1D array of noninteracting particles in 3D harmonic oscillator wells.
* `SodiumChlorideCrystal`: FCC crystal of sodium chloride
* `LennardJonesCluster`: A non-periodic rectilinear grid of Lennard-Jones particles in a harmonic restraining potential.
* `LennardJonesFluid`: periodic rectilinear grid of Lennard-Jones particles.
* `CustomLennardJonesFluidMixture`: A periodic rectilinear grid of Lennard-Jones particled, but implemented via CustomBondForce and NonbondedForce.
* `WCAFluid`: Weeks-Chandler-Andersen system
* `IdealGas`: An 'ideal gas' of noninteracting particles in a periodic box.
* `WaterBox`: A water box test system (multiple water models possible, with variants including discharged versions, flexible and constrained waters, etc.)
* `AlanineDipeptideVacuum`: Alanine dipeptide in vacuum.
* `AlanineDipeptideImplicit`: Alanine dipeptide in implicit solvent
* `AlanineDipeptideExplicit`: Alanine dipeptide in explicit solvent.
* `LysozymeImplicit`: T4 lysozyme L99A with p-xylene ligand in OBC GBSA
* `MethanolBox`: A periodic box of methanol
* `MolecularIdealGas`: Molecular ideal gas (methanol box).
* `CustomGBForceSystem`: A system of particles with a CustomGBForce.
* `AMOEBAIonBox`: A single Ca2 ion in a water box.
* Lots of alchemically-modified systems

#### Stability tests
Stability tests ensure that simulations of typical systems of interest are stable and robust.

##### Old validation suite tests
Tests are run by running the `runAllTests.sh` script, which runs the systems in `systems/` through four different kinds of tests in different subsets:
* **Consistency between platforms**: `TestForces.py`
* **Energy-force consistency**: `TestEnergyForces.py`
* **Energy conservation**: `TestVerletEnergyConservation.py`
* **Thermostability** (really just a thermostat test): `TestLangevinThermostability.py`

##### Test systems from old OpenMM validation suite
* 1PLX
* bpti
* ala10
* dnaDickerson
* ubiquitin
* spectrin
* sv582
* sv583
* 6TNA
* CheY
* lambda repressor, d14a
* 1yrf
* 1ubq
* 1not



#### Numerical tests
Numerical tests ensure that the various components of OpenMM are numerically well-behaved, and agree between `Platform`s.

### Comparison tests

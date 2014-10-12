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

### Numerical tests
* Is it appropriate to make comparisons using relative errors in forces and energies, or absolute errors?
* If absolute errors, what should the criteria be for energies and forces?
  * For energies, `kT` is a meaningful scale. For example, we may want errors less than `0.01*kT`, to ensure equilibrium probabilities are not corrupted substantially.
  * What would the criteria be for forces?
* How should finite-difference tests of energies and forces be conducted?

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

#### Checkpoint
* **testCheckpoint.** Ensure that trajectories generated after saving a checkpoint and after resuming from the checkpoint generate identical positions, velocities, and parameters.
* **testSetState.** Perform the same test as `testCheckpoint` but utilize the `Context::set{State|Parameter|PeriodicBoxVectors|Positions}` interface instead.

#### CustomAngleForce
* **testAngles.** Compare a four-particle system with two angles to `HarmonicAngleForce` over 10 random geometries to ensure energies and forces (rel tol `1e-5`) match. Change parameters `theta0` and `K` and test again.

#### CustomBondForce
* **testBonds.** For a single geometry of three atoms and two bonds, test forces and energies are correct (rel tol `1e-5`) for a harmonic bond. Change `r0` and `K` and test again.
* **testManyParameters.** For a two particle system connected by a bond whose energy is linear in `r`, check that the energy and forces (rel tol `1e-5`) is correct when 9 user-specified parameters are summed to give the multiplicative constant.

#### CustomCompoundBondForce
* **testBond.** For a four-particle system, test 10 random geometries to ensure a `CustomCompoundBondForce` combination of bond-angle-torsion terms matches that expected from a standard combination of bond, angle, and torsion forces.
* **testPositionDependence.** For a two-particle system interacting via a harmonic bond and position-dependent field per particle, ensure energy and force (rel tol `1e-5`) are correct.

#### CustomExternalForce
* **testForce.** For a harmonic osccilator in the y-axis, test that a three-particle system obtains the correct energies and forces (rel tol `1e-5`) for two specific geometries and two sets of parameters.
* **testManuParameters.** For a single particle, test a three-dimensional harmonic oscillator with separate x, y, and z-axis parameters obtains the correct energy and forces (rel tol `1e-5`).

#### CustomGBForce
* **testOBC.** Compare a `CustomGBForce` implementation of the OBC GBSA scheme (with ACE surface area term) with the `GBSAOBCForce` implementation for a system of 70 diatomic molecules, testing two parameter sets and a single random configuration of particles.
* **testMembrane.** Construct a `CustomGBForce` for an implicit membrane model and check whether a small step (0.1A) in the direction of the gradient causes the energy to change by the expected amount (rel tol `1e-3`).
* **testTabulatedFunction.** Test that a tabulated function with 2.5A spacing is reproduced by `CustomGBForce` in both forces (rel tol 0.1) and potential (rel tol 0.02).
* **testMultipleChainRules.** Test that successive application of the chain rule works correctly for a simple chain of linear functions, comparing energies (rel tol 0.02) and forces (rel tol `1e-4`).
* **testPositionDependence.** For a pair of particles, test that a `CustomGBForce` that depends both on interatomic distance and (x,y) position in space achieves correct energy (rel tol 0.02) and force (rel tol `1e-4`) for 5 random configurations. A small step in the direction of the gradient is also tested for expected energy change (rel tol `1e-3`).
* **testExclusions.** A system of two particles is tested for correct energy (rel tol `1e-4`) and forces (rel tol `1e-4`) when cycling through implemented particle exclusion options (`ParticlePair`, `ParticlePairNoExclusions`) for computed values and energy terms. A small step (0.01A) in the direction of the gradient is also tested for expected energy change (abs tol `1e-3` kJ/mol).

#### CustomHbondForce
* **testHbond.** For a system of 5 particles, implement a `CustomHbondForce` that depends on a distance, angle, and torsion, and check that the energy (rel tol 1e-5) and force (rel tol 1e-5) are reproduced for 10 random configurations.
* **testExclusions.** For a system of 3 particles, test a simple distance-dependent `CustomHbondForce` containing a single exclusion to ensure that the energy (rel tol 1e-5) and force (rel tol 1e-5) are reproduced for a single conformation.
* **testCutoff.** For a system of 3 particles, test a simple distance-dependent `CustomHbondForce` with a distance cutoff to ensure that the energy (rel tol 1e-5) and force (rel tol 1e-5) are reproduced for a single conformation.
* **testCustomFunctions.** Test that a linearly interpolated function of distance is correctly reproduced for a single geometry of three particles, testing energy (rel tol 1e-5) and force (rel tol 1e-5).

#### CustomIntegrator
* **testSingleBond.** Test a leapfrog Verlet integrator with a two-particle system joined by a harmonic bond by comparing the trajectory over 1000 steps of 10 fs timestep dynamics with the analytical solution, testing positions (rel tol 1e-4), velocities (rel tol 1e-4), and potential (rel tol 1e-4).
* **testConstraints.** For an 8-particle chain, test a leapfrog Verlet integrator with SHAKE constraints (tol 1e-5) starting from deterministic positions and random velocities over 1000 timesteps of 2 fs to ensure bond lengths remain constrained (rel tol 2e-5). Check that energy is approximately equal to initial energy (rel tol 0.01).
* **testVelocityConstraints.** For a 10-particle system, with three particles constrained by SHAKE (tol 1e-5), three particles constrained by SETTLE, and the rest by CCMA (tol 1e-5), integrate using a velocity Verlet integrator with SHAKE/RATTLE starting from a fixed initial geometry and random initial velocities for 1000 steps of 2 fs, checking that we preserve distances (rel tol 2e-5), velocities along constraints are zero (rel tol 2e-5), and that the total energy is approximately preserved (rel tol 0.01).
* **testConstrainedMasslessParticles.** For a system with two particles connected by a constraint, in which one particle has its mass set to zero (fixing it in space), ensure an Exception is thrown. When both particles have masses set to zero (fixing them in space), no exception is thrown. Integration should keep both velocities set to zero.
* **testWithThermostat.** Test an 8-particle alternately charged LJ system with a leapfrog Verlet integrator and Andersen thermostat to make sure the average instantaneous kinetic energy over 150 ps is correct for the thermostat temperature (rel tol 0.1).
* **testMonteCarlo.** Test a simple Metropolis Monte Carlo scheme using Gaussian proposal with a system of two harmonically bonded particles to ensure the binned distribution of distances (bin size 0.4, 100 bins) is approximately satisfied (rel tol 0.01).
* **testSum.** Test that specifying a kinetic energy expression for a leapfrog Verlet integrator gives the expected kinetic energy for a 200-particle system (rel tol 1e-5).
* **testParameter.**
* **testRandomDistributions.**
* **testPerDofVariables.**
* **testForceGroups.**
* **testRespa.**
* **testMergedRandoms.**

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

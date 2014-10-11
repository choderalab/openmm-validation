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


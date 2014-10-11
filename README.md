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

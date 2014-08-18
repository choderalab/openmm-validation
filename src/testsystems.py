"""
Module to generate Systems and positions for simple reference molecular systems for testing.

DESCRIPTION

This module provides functions for building a number of test systems of varying complexity,
useful for testing both OpenMM and various codes based on pyopenmm.

Note that the PYOPENMM_SOURCE_DIR must be set to point to where the PyOpenMM package is unpacked.

EXAMPLES

Create a 3D harmonic oscillator.

>>> import testsystems
>>> ho = testsystems.HarmonicOscillator()
>>> system, positions = ho.system, ho.positions

See list of methods for a complete list of provided test systems.

COPYRIGHT

@author Randall J. Radmer <radmer@stanford.edu>
@author John D. Chodera <jchodera@gmail.com>

All code in this repository is released under the GNU General Public License.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

TODO

* Add units checking code to check arguments.
* Change default arguments to Quantity objects, rather than None?

"""

import numpy as np
import numpy.random
import math
import copy
import scipy.special

import simtk.openmm as mm
import simtk.unit as unit
import simtk.openmm.app as app

from repex.utils import get_data_filename
from repex.thermodynamics import ThermodynamicState

kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA 


#=============================================================================================
# Abstract base class for test systems
#=============================================================================================

class TestSystem(object):
    """Abstract base class for test systems, demonstrating how to implement a test system.

    Parameters
    ----------
    
    Attributes
    ----------
    system : simtk.openmm.System
        Openmm system with the harmonic oscillator
    positions : list
        positions of harmonic oscillator

    Notes
    -----

    Unimplemented methods will default to the base class methods, which raise a NotImplementedException.

    Examples
    --------

    Create a test system.

    >>> testsystem = TestSystem()
    
    Retrieve System object.

    >>> system = testsystem.system

    Retrieve the positions.
    
    >>> positions = testsystem.positions

    Serialize system and positions to XML (to aid in debugging).

    >>> (system_xml, positions_xml) = testsystem.serialize()

    """
    def __init__(self, temperature=None, pressure=None):
        """Abstract base class for test system.

        Parameters
        ----------

        temperature : simtk.unit.Quantity, optional, units compatible with simtk.unit.kelvin
            The temperature of the system.

        pressure : simtk.unit.Quantity, optional, units compatible with simtk.unit.atmospheres
            The pressure of the system.

        """
        
        # Create an empty system object.
        self._system = mm.System()

        # Store positions.
        self._positions = unit.Quantity(np.zeros([0,3], np.float), unit.nanometers)

        # Store thermodynamic parameters.
        self._temperature = temperature
        self._pressure = pressure
        
        return

    @property
    def system(self):
        """The simtk.openmm.System object corresponding to the test system."""
        return copy.deepcopy(self._system)

    @system.setter
    def system(self, value):
        self._system = value

    @system.deleter
    def system(self):
        del self._system

    @property
    def positions(self):
        """The simtk.unit.Quantity object containing the particle positions, with units compatible with simtk.unit.nanometers."""
        return copy.deepcopy(self._positions)

    @positions.setter
    def positions(self, value):
        self._positions = value
    
    @positions.deleter
    def positions(self):
        del self._positions

    @property
    def analytical_properties(self):
        """A list of available analytical properties, accessible via 'get_propertyname(thermodynamic_state)' calls."""
        return [ method[4:] for method in dir(self) if (method[0:4]=='get_') ]

    def reduced_potential_expectation(self, state_sampled_from, state_evaluated_in):
        """Calculate the expected potential energy in state_sampled_from, divided by kB * T in state_evaluated_in.
        
        Notes
        -----
        
        This is not called get_reduced_potential_expectation because this function
        requires two, not one, inputs.
        """
        
        if hasattr(self, "get_potential_expectation"):
            U = self.get_potential_expectation(state_sampled_from)
            U_red = U / (kB * state_evaluated_in.temperature)
            return U_red
        else:
            raise AttributeError("Cannot return reduced potential energy because system lacks get_potential_expectation")

    def serialize(self):
        """Return the System and positions in serialized XML form.

        Returns
        -------
        
        system_xml : str
            Serialized XML form of System object.
            
        state_xml : str
            Serialized XML form of State object containing particle positions.

        """

        from simtk.openmm import XmlSerializer
        
        # Serialize System.
        system_xml = XmlSerializer.serialize(self._system)

        # Serialize positions via State.
        if self._system.getNumParticles() == 0:
        # Cannot serialize the State of a system with no particles.
            state_xml = None
        else:
            platform = mm.Platform.getPlatformByName('Reference')
            integrator = mm.VerletIntegrator(1.0 * unit.femtoseconds)
            context = mm.Context(self._system, integrator, platform)
            context.setPositions(self._positions)
            state = context.getState(getPositions=True)
            del context, integrator
            state_xml = XmlSerializer.serialize(state)

        return (system_xml, state_xml)

    @property
    def name(self):
        """The name of the test system."""
        return self.__class__.__name__

#=============================================================================================
# 3D harmonic oscillator
#=============================================================================================

class HarmonicOscillator(TestSystem):
    """Create a 3D harmonic oscillator, with a single particle confined in an isotropic harmonic well.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=90.0 * unit.kilocalories_per_mole/unit.angstrom**2
        harmonic restraining potential
    mass : simtk.unit.Quantity, optional, default=39.948 * unit.amu
        particle mass
    
    Attributes
    ----------
    system : simtk.openmm.System
        Openmm system with the harmonic oscillator
    positions : list
        positions of harmonic oscillator

    Notes
    -----

    The natural period of a harmonic oscillator is T = sqrt(m/K), so you will want to use an
    integration timestep smaller than ~ T/10.

    The standard deviation in position in each dimension is sigma = (kT / K)^(1/2)

    The expectation and standard deviation of the potential energy of a 3D harmonic oscillator is (3/2)kT.

    Examples
    --------

    Create a 3D harmonic oscillator with default parameters:

    >>> ho = HarmonicOscillator()
    >>> (system, positions) = ho.system, ho.positions

    Create a harmonic oscillator with specified mass and spring constant:

    >>> mass = 12.0 * unit.amu
    >>> K = 1.0 * unit.kilocalories_per_mole / unit.angstroms**2
    >>> ho = HarmonicOscillator(K=K, mass=mass)
    >>> (system, positions) = ho.system, ho.positions

    Get a list of the available analytically-computed properties.

    >>> print ho.analytical_properties
    ['potential_expectation', 'potential_standard_deviation']

    Compute the potential expectation and standard deviation

    >>> import simtk.unit as u
    >>> thermodynamic_state = ThermodynamicState(temperature=298.0*u.kelvin, system=system)
    >>> potential_mean = ho.get_potential_expectation(thermodynamic_state)
    >>> potential_stddev = ho.get_potential_standard_deviation(thermodynamic_state)
    
    """
    
    def __init__(self, K=100.0 * unit.kilocalories_per_mole / unit.angstroms**2, mass=39.948 * unit.amu, **kwargs):

        TestSystem.__init__(self, kwargs)

        # Create an empty system object.
        system = mm.System()

        # Add the particle to the system.
        system.addParticle(mass)

        # Set the positions.
        positions = unit.Quantity(np.zeros([1,3], np.float32), unit.angstroms)

        # Add a restrining potential centered at the origin.
        force = mm.CustomExternalForce('(K/2.0) * (x^2 + y^2 + z^2)')
        force.addGlobalParameter('K', K)
        force.addParticle(0, [])
        system.addForce(force)
        
        self.K, self.mass = K, mass
        self.system, self.positions = system, positions
        
        # Number of degrees of freedom.
        self.ndof = 3

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (3./2.) * kB * state.temperature
        
    def get_potential_standard_deviation(self, state):
        """Return the standard deviation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------
        
        potential_stddev : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            potential energy standard deviation if implemented, or else None
        
        """

        return (3./2.) * kB * state.temperature

class PowerOscillator(TestSystem):
    """Create a 3D Power oscillator, with a single particle confined in an isotropic x^b well.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=100.0
        harmonic restraining potential.  The units depend on the power, 
        so we accept unitless inputs and add units of the form 
        unit.kilocalories_per_mole / unit.angstrom ** b
    mass : simtk.unit.Quantity, optional, default=39.948 * unit.amu
        particle mass
    
    Attributes
    ----------
    system : simtk.openmm.System
        Openmm system with the harmonic oscillator
    positions : list
        positions of harmonic oscillator

    Notes
    -----

    Here we assume a potential energy of the form U(x) = k * x^b.  

    By the generalized equipartition theorem, the expectation of the 
    potential energy is 3 kT / b.
    
    """
    
    def __init__(self, K=100.0, b=2.0, mass=39.948 * unit.amu, **kwargs):

        TestSystem.__init__(self, kwargs)
        
        K = K * unit.kilocalories_per_mole / unit.angstroms ** b

        # Create an empty system object.
        system = mm.System()

        # Add the particle to the system.
        system.addParticle(mass)

        # Set the positions.
        positions = unit.Quantity(np.zeros([1,3], np.float32), unit.angstroms)

        # Add a restrining potential centered at the origin.
        force = mm.CustomExternalForce('(K) * (x^%d + y^%d + z^%d)' %(b, b, b))
        force.addGlobalParameter('K', K)
        force.addParticle(0, [])
        system.addForce(force)
        
        self.K, self.mass = K, mass
        self.b = b
        self.system, self.positions = system, positions
        
        # Number of degrees of freedom.
        self.ndof = 3

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (3.) * kB * state.temperature / self.b

    def _get_power_expectation(self, state, n):
        """Return the power of x^n.  Not currently used"""
        b = 1.0 * self.b
        beta = (1.0 * kB * state.temperature) ** -1.
        gamma = scipy.special.gamma
        return (self.K * beta) ** (-n / b) * gamma((n + 1.) / b) / gamma(1. / b)

    @classmethod
    def reduced_potential(cls, beta, a, b, a2, b2):
        gamma = scipy.special.gamma
        reduced_u = 3 * a2 * (a * beta) ** (-b2 / b) * gamma((b2 + 1.) / b) / gamma(1. / b) * beta
        return reduced_u

#=============================================================================================
# Diatomic molecule
#=============================================================================================

class Diatom(TestSystem):
    """Create a free diatomic molecule with a single harmonic bond between the two atoms.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=290.1 * unit.kilocalories_per_mole / unit.angstrom**2
        harmonic bond potential.  default is GAFF c-c bond
    r0 : simtk.unit.Quantity, optional, default=1.550 * unit.amu
        bond length.  Default is Amber GAFF c-c bond.
    constraint : bool, default=False
        if True, the bond length will be constrained
    m1 : simtk.unit.Quantity, optional, default=12.01 * unit.amu
        particle1 mass
    m2 : simtk.unit.Quantity, optional, default=12.01 * unit.amu
        particle2 mass
    use_central_potential : bool, optional, default=False
        if True, a soft central potential will also be added to keep the system from drifting away        
    

    Notes
    -----

    The natural period of a harmonic oscillator is T = sqrt(m/K), so you will want to use an
    integration timestep smaller than ~ T/10.

    Examples
    --------

    Create a Diatom:

    >>> diatom = Diatom()
    >>> system, positions = diatom.system, diatom.positions

    Create a Diatom with constraint in a central potential
    >>> diatom = Diatom(constraint=True, use_central_potential=True)
    >>> system, positions = diatom.system, diatom.positions

    """

    def __init__(self, 
        K=290.1 * unit.kilocalories_per_mole / unit.angstrom**2,
        r0=1.550 * unit.angstroms, 
        m1=39.948 * unit.amu,
        m2=39.948 * unit.amu,
        constraint=False,
        use_central_potential=False):

        # Create an empty system object.
        system = mm.System()

        # Add two particles to the system.
        system.addParticle(m1)
        system.addParticle(m2)

        # Add a harmonic bond.
        force = mm.HarmonicBondForce()
        force.addBond(0, 1, r0, K)
        system.addForce(force)

        if constraint:
            # Add constraint between particles.
            system.addConstraint(0, 1, r0)
        
        # Set the positions.
        positions = unit.Quantity(np.zeros([2,3], np.float32), unit.angstroms)
        positions[1,0] = r0

        if use_central_potential:
            # Add a central restraining potential.
            Kcentral = 1.0 * unit.kilocalories_per_mole / unit.nanometer**2
            force = mm.CustomExternalForce('(Kcentral/2.0) * (x^2 + y^2 + z^2)')
            force.addGlobalParameter('Kcentral', Kcentral)
            force.addParticle(0, [])
            force.addParticle(1, [])    
            system.addForce(force)

        self.system, self.positions = system, positions
        self.K, self.r0, self.m1, self.m2, self.constraint, self.use_central_potential = K, r0, m1, m2, constraint, use_central_potential
        
        # Store number of degrees of freedom.
        self.ndof = 6 - 1*constraint

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (self.ndof/2.) * kB * state.temperature

#=============================================================================================
# Diatomic fluid
#=============================================================================================

class DiatomicFluid(TestSystem):
    """Create a diatomic fluid.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=290.1 * unit.kilocalories_per_mole / unit.angstrom**2
        harmonic bond potential.  default is GAFF c-c bond
    r0 : simtk.unit.Quantity, optional, default=1.550 * unit.amu
        bond length.  Default is Amber GAFF c-c bond.
    constraint : bool, default=False
        if True, the bond length will be constrained
    m1 : simtk.unit.Quantity, optional, default=12.01 * unit.amu
        particle1 mass
    m2 : simtk.unit.Quantity, optional, default=12.01 * unit.amu
        particle2 mass
    epsilon : simtk.unit.Quantity, optional, default=0.1700 * unit.kilocalories_per_mole
        particle Lennard-Jones well depth
    sigma : simtk.unit.Quantity, optional, default=1.8240 * unit.angstroms
        particle Lennard-Jones sigma
    charge : simtk.unit.Quantity, optional, default=0.0 * unit.elementary_charge
        charge to place on atomic centers to create a dipole
    cutoff : simtk.unit.Quantity, optional, default=None
        if specified, the specified cutoff will be used; otherwise, half the box width will be used
    switch : bool, optional, default=True
        flag to use nonbonded switching function
    switch_width : simtk.unit.Quantity with units compatible with angstroms, optional, default=0.2*unit.angstroms
        switching function is turned on at cutoff - switch_width
    dispersion_correction : bool, optional, default=True
        if True, will use analytical dispersion correction (if not using switching function)
    nx, ny, nz : int, optional, default=6
        number of molecules in x, y, and z dimension


    Notes
    -----

    The natural period of a harmonic oscillator is T = sqrt(m/K), so you will want to use an
    integration timestep smaller than ~ T/10.

    Examples
    --------

    Create an uncharged Diatomic fluid.

    >>> diatom = DiatomicFluid()
    >>> system, positions = diatom.system, diatom.positions

    Create a dipolar fluid.

    >>> diatom = DiatomicFluid(charge=1.0*unit.elementary_charge)
    >>> system, positions = diatom.system, diatom.positions

    Create a Diatomic fluid with constraints instead of harmonic bonds

    >>> diatom = DiatomicFluid(constraint=True)
    >>> system, positions = diatom.system, diatom.positions

    Specify a different system size.

    >>> diatom = DiatomicFluid(constraint=True, nx=8, ny=8, nz=8)
    >>> system, positions = diatom.system, diatom.positions

    TODO
    ----
    
    Add subrandom selection of orientations.
    Add multiple electrostatics treatments.

    """

    def __init__(self, 
        K=424.0* unit.kilocalories_per_mole / unit.angstrom**2,
        r0=1.383 * unit.angstroms,  
        m1=14.01 * unit.amu,
        m2=14.01 * unit.amu,
        epsilon=0.1700 * unit.kilocalories_per_mole, 
        sigma=1.8240 * unit.angstroms,
        charge=0.0 * unit.elementary_charge,
        switch=True,
        switch_width=0.5*unit.angstroms,
        cutoff=None,
        constraint=False,
        dispersion_correction=True,
        nx=7, ny=7, nz=7):
        
        nmolecules = nx * ny * nz
        natoms = 2 * nmolecules

        # Create an empty system object.
        system = mm.System()

        # Add particles to the system.
        for molecule_index in range(nmolecules):
            system.addParticle(m1)
            system.addParticle(m2)

        if constraint:
            # Add constraint between particles.
            for molecule_index in range(nmolecules):
                system.addConstraint(2*molecule_index+0, 2*molecule_index+1, r0)
        else:
            # Add a harmonic bonds.
            force = mm.HarmonicBondForce()
            for molecule_index in range(nmolecules):
                force.addBond(2*molecule_index+0, 2*molecule_index+1, r0, K)
            system.addForce(force)

        # Set up nonbonded interactions.
        nb = mm.NonbondedForce()
            
        positions = unit.Quantity(np.zeros([natoms,3],np.float32), unit.angstrom)

        maxX = 0.0 * unit.angstrom
        maxY = 0.0 * unit.angstrom
        maxZ = 0.0 * unit.angstrom

        dx = r0/2
        dy = 0 * sigma
        dz = 0 * sigma

        scaleStepSizeX = 2.0
        scaleStepSizeY = 2.0
        scaleStepSizeZ = 2.0

        atom_index = 0
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    # Compute particle center.
                    x = sigma*scaleStepSizeX*ii
                    y = sigma*scaleStepSizeY*jj
                    z = sigma*scaleStepSizeZ*kk

                    # Add particle 1.
                    nb.addParticle(+charge, sigma, epsilon)

                    positions[atom_index,0] = x - dx
                    positions[atom_index,1] = y - dy
                    positions[atom_index,2] = z - dz
                    atom_index += 1

                    # Add particle 2.
                    nb.addParticle(-charge, sigma, epsilon)

                    positions[atom_index,0] = x + dx
                    positions[atom_index,1] = y + dy
                    positions[atom_index,2] = z + dz
                    atom_index += 1
                    
                    # Wrap center positions as needed.
                    if x>maxX: maxX = x
                    if y>maxY: maxY = y
                    if z>maxZ: maxZ = z

        # Add exceptions.
        for molecule_index in range(nmolecules):
            nb.addException(2*molecule_index+0, 2*molecule_index+1, 0.0*charge*charge, sigma, 0.0*epsilon)

        system.addForce(nb)
                    
        # Set periodic box vectors.
        x = maxX+2*sigma*scaleStepSizeX
        y = maxY+2*sigma*scaleStepSizeY
        z = maxZ+2*sigma*scaleStepSizeZ

        if not cutoff:
            min_width = min(maxX, maxY, maxZ)
            cutoff = min_width / 2.0 - 0.01 * unit.angstroms

        nb.setNonbondedMethod(mm.NonbondedForce.PME)
        nb.setUseDispersionCorrection(dispersion_correction)
        nb.setUseSwitchingFunction(switch)
        nb.setCutoffDistance(cutoff)
        nb.setSwitchingDistance(cutoff-switch_width)

        a = unit.Quantity((x,                0*unit.angstrom, 0*unit.angstrom))
        b = unit.Quantity((0*unit.angstrom,                y, 0*unit.angstrom))
        c = unit.Quantity((0*unit.angstrom, 0*unit.angstrom, z))
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Store number of degrees of freedom.
        self.ndof = 3*natoms - nmolecules*constraint

        # Store system and positions.
        self._system = system
        self._positions = positions

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (self.ndof/2.) * kB * state.temperature


class UnconstrainedDiatomicFluid(DiatomicFluid):
    """
    Examples
    --------

    Create an unconstrained diatomic fluid.

    >>> test = UnconstrainedDiatomicFluid()
    >>> system, positions = test.system, test.positions
    
    """
    def __init__(self, *args, **kwargs):
       super(UnconstrainedDiatomicFluid, self).__init__(constraint=False, *args, **kwargs)

class ConstrainedDiatomicFluid(DiatomicFluid):
    """
    Examples
    --------

    Create an constrained diatomic fluid.

    >>> test = ConstrainedDiatomicFluid()
    >>> system, positions = test.system, test.positions
    
    """
    def __init__(self, *args, **kwargs):
       super(ConstrainedDiatomicFluid, self).__init__(constraint=True, *args, **kwargs)

class DipolarFluid(DiatomicFluid):
    """
    Examples
    --------

    Create a dipolar fluid.

    >>> test = DipolarFluid()
    >>> system, positions = test.system, test.positions
    
    """
    def __init__(self, *args, **kwargs):
       super(DipolarFluid, self).__init__(charge=0.25*unit.elementary_charge, *args, **kwargs)

class UnconstrainedDipolarFluid(DipolarFluid):
    """
    Examples
    --------

    Create a dipolar fluid.

    >>> test = UnconstrainedDipolarFluid()
    >>> system, positions = test.system, test.positions

    """
    def __init__(self, *args, **kwargs):
       super(UnconstrainedDipolarFluid, self).__init__(constraint=False, *args, **kwargs)

class ConstrainedDipolarFluid(DipolarFluid):
    """
    Examples
    --------

    Create a dipolar fluid.

    >>> test = ConstrainedDipolarFluid()
    >>> system, positions = test.system, test.positions

    """
    def __init__(self, *args, **kwargs):
       super(ConstrainedDipolarFluid, self).__init__(constraint=True, *args, **kwargs)

#=============================================================================================
# Constraint-coupled harmonic oscillator
#=============================================================================================

class ConstraintCoupledHarmonicOscillator(TestSystem):
    """Create a pair of particles in 3D harmonic oscillator wells, coupled by a constraint.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=1.0 * unit.kilojoules_per_mole / unit.nanometer**2
        harmonic restraining potential
    d : simtk.unit.Quantity, optional, default=1.0 * unit.nanometer
        distance between harmonic oscillators.  Default is Amber GAFF c-c bond.
    mass : simtk.unit.Quantity, default=39.948 * unit.amu
        particle mass

    Attributes
    ----------
    system : simtk.openmm.System
    positions : list

    Notes
    -----

    The natural period of a harmonic oscillator is T = sqrt(m/K), so you will want to use an
    integration timestep smaller than ~ T/10.

    Examples
    --------

    Create a constraint-coupled harmonic oscillator with specified mass, distance, and spring constant.

    >>> ccho = ConstraintCoupledHarmonicOscillator()
    >>> mass = 12.0 * unit.amu
    >>> d = 5.0 * unit.angstroms
    >>> K = 1.0 * unit.kilocalories_per_mole / unit.angstroms**2
    >>> ccho = ConstraintCoupledHarmonicOscillator(K=K, d=d, mass=mass)
    >>> system, positions = ccho.system, ccho.positions
    """

    def __init__(self,
        K=1.0 * unit.kilojoules_per_mole/unit.nanometer**2,
        d=1.0 * unit.nanometer,
        mass=39.948 * unit.amu):

        # Create an empty system object.
        system = mm.System()

        # Add particles to the system.
        system.addParticle(mass)
        system.addParticle(mass)

        # Set the positions.
        positions = unit.Quantity(np.zeros([2,3], np.float32), unit.angstroms)
        positions[1,0] = d

        # Add a restrining potential centered at the origin.
        force = mm.CustomExternalForce('(K/2.0) * ((x-d)^2 + y^2 + z^2)')
        force.addGlobalParameter('K', K)
        force.addPerParticleParameter('d')
        force.addParticle(0, [0.0])
        force.addParticle(1, [d / unit.nanometers])
        system.addForce(force)

        # Add constraint between particles.
        system.addConstraint(0, 1, d)

        # Add a harmonic bond force as well so minimization will roughly satisfy constraints.
        force = mm.HarmonicBondForce()
        K = 10.0 * unit.kilocalories_per_mole / unit.angstrom**2 # force constant
        force.addBond(0, 1, d, K)
        system.addForce(force)

        self.system, self.positions = system, positions
        self.K, self.d, self.mass = K, d, mass

#=============================================================================================
# Harmonic oscillator array
#=============================================================================================

class HarmonicOscillatorArray(TestSystem):
    """Create a 1D array of noninteracting particles in 3D harmonic oscillator wells.

    Parameters
    ----------
    K : simtk.unit.Quantity, optional, default=90.0 * unit.kilocalories_per_mole/unit.angstroms**2
        harmonic restraining potential
    d : simtk.unit.Quantity, optional, default=1.0 * unit.nanometer
        distance between harmonic oscillators.  Default is Amber GAFF c-c bond.
    mass : simtk.unit.Quantity, default=39.948 * unit.amu
        particle mass
    N : int, optional, default=5
        Number of harmonic oscillators
    
    Attributes
    ----------
    system : simtk.openmm.System
    positions : list

    Notes
    -----

    The natural period of a harmonic oscillator is T = sqrt(m/K), so you will want to use an
    integration timestep smaller than ~ T/10.

    Examples
    --------

    Create a constraint-coupled 3D harmonic oscillator with default parameters.

    >>> ho_array = HarmonicOscillatorArray()
    >>> mass = 12.0 * unit.amu
    >>> d = 5.0 * unit.angstroms
    >>> K = 1.0 * unit.kilocalories_per_mole / unit.angstroms**2
    >>> ccho = HarmonicOscillatorArray(K=K, d=d, mass=mass)
    >>> system, positions = ccho.system, ccho.positions
    """

    def __init__(self, K=90.0 * unit.kilocalories_per_mole/unit.angstroms**2,
        d=1.0 * unit.nanometer,
        mass=39.948 * unit.amu ,
        N=5):

        # Create an empty system object.
        system = mm.System()

        # Add particles to the system.
        for n in range(N):
            system.addParticle(mass)

        # Set the positions for a 1D array of particles spaced d apart along the x-axis.
        positions = unit.Quantity(np.zeros([N,3], np.float32), unit.angstroms)
        for n in range(N):
            positions[n,0] = n*d

        # Add a restrining potential for each oscillator.
        force = mm.CustomExternalForce('(K/2.0) * ((x-x0)^2 + y^2 + z^2)')
        force.addGlobalParameter('K', K)
        force.addPerParticleParameter('x0')
        for n in range(N):
            parameters = (d*n / unit.nanometers, )
            force.addParticle(n, parameters)
        system.addForce(force)

        self.system, self.positions = system, positions
        self.K, self.d, self.mass, self.N = K, d, mass, N
        self.ndof = 3*N

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------

        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------

        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.

        """

        return (self.ndof/2.) * kB * state.temperature

    def get_potential_standard_deviation(self, state):
        """Return the standard deviation of the potential energy, computed analytically or numerically.

        Arguments
        ---------

        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------

        potential_stddev : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            potential energy standard deviation if implemented, or else None

        """

        return (self.ndof/2.) * kB * state.temperature

#=============================================================================================
# Sodium chloride FCC crystal.
#=============================================================================================

class SodiumChlorideCrystal(TestSystem):
    """Create an FCC crystal of sodium chloride.

    Each atom is represented by a charged Lennard-Jones sphere in an Ewald lattice.

    switch : bool, optional, default=True
        flag to use nonbonded switching function
    switch_width : simtk.unit.Quantity with units compatible with angstroms, optional, default=0.2*unit.angstroms
        switching function is turned on at cutoff - switch_width
    dispersion_correction : bool, optional, default=True
        if True, will use analytical dispersion correction (if not using switching function)

    Notes
    -----

    TODO

    * Lennard-Jones interactions aren't correctly being included now, due to LJ cutoff.  Fix this by hard-coding LJ interactions?
    * Add nx, ny, nz arguments to allow user to specify replication of crystal unit in x,y,z.
    * Choose more appropriate lattice parameters and lattice spacing.

    Examples
    --------

    Create sodium chloride crystal unit.

    >>> crystal = SodiumChlorideCrystal()
    >>> system, positions = crystal.system, crystal.positions
    """
    def __init__(self, switch=True, switch_width=0.2*unit.angstroms, dispersion_correction=True):
        # Set default parameters (from Tinker).
        mass_Na     = 22.990 * unit.amu
        mass_Cl     = 35.453 * unit.amu
        q_Na        = 1.0 * unit.elementary_charge
        q_Cl        =-1.0 * unit.elementary_charge
        sigma_Na    = 3.330445 * unit.angstrom
        sigma_Cl    = 4.41724 * unit.angstrom
        epsilon_Na  = 0.002772 * unit.kilocalorie_per_mole
        epsilon_Cl  = 0.118 * unit.kilocalorie_per_mole

        # Create system
        system = mm.System()

        # Set box vectors.
        box_size = 5.628 * unit.angstroms # box width
        a = unit.Quantity(np.zeros([3]), unit.nanometers); a[0] = box_size
        b = unit.Quantity(np.zeros([3]), unit.nanometers); b[1] = box_size
        c = unit.Quantity(np.zeros([3]), unit.nanometers); c[2] = box_size
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Create nonbonded force term.
        force = mm.NonbondedForce()

        # Set interactions to be periodic Ewald.
        force.setNonbondedMethod(mm.NonbondedForce.Ewald)

        # Set cutoff to be less than one half the box length.
        cutoff = box_size / 2.0 * 0.99
        force.setCutoffDistance(cutoff)

        # Set treatment.
        force.setUseDispersionCorrection(dispersion_correction)
        force.setUseSwitchingFunction(switch)
        force.setSwitchingDistance(cutoff-switch_width)

        # Allocate storage for positions.
        natoms = 2
        positions = unit.Quantity(np.zeros([natoms,3], np.float32), unit.angstroms)

        # Add sodium ion.
        system.addParticle(mass_Na)
        force.addParticle(q_Na, sigma_Na, epsilon_Na)
        positions[0,0] = 0.0 * unit.angstrom
        positions[0,1] = 0.0 * unit.angstrom
        positions[0,2] = 0.0 * unit.angstrom

        # Add chloride atom.
        system.addParticle(mass_Cl)
        force.addParticle(q_Cl, sigma_Cl, epsilon_Cl)
        positions[1,0] = 2.814 * unit.angstrom
        positions[1,1] = 2.814 * unit.angstrom
        positions[1,2] = 2.814 * unit.angstrom

        # Add nonbonded force term to the system.
        system.addForce(force)

        self.system, self.positions = system, positions

#=============================================================================================
# Lennard-Jones cluster
#=============================================================================================

class LennardJonesCluster(TestSystem):
    """Create a non-periodic rectilinear grid of Lennard-Jones particles in a harmonic restraining potential.

    Parameters
    ----------
    nx : int, optional, default=3
        number of particles in the x direction
    ny : int, optional, default=3
        number of particles in the y direction
    nz : int, optional, default=3
        number of particles in the z direction
    K : simtk.unit.Quantity, optional, default=1.0 * unit.kilojoules_per_mole/unit.nanometer**2
        harmonic restraining potential

    Examples
    --------

    Create Lennard-Jones cluster.

    >>> cluster = LennardJonesCluster()
    >>> system, positions = cluster.system, cluster.positions

    Create default 3x3x3 Lennard-Jones cluster in a harmonic restraining potential.

    >>> cluster = LennardJonesCluster(nx=10, ny=10, nz=10)
    >>> system, positions = cluster.system, cluster.positions
    """
    def __init__(self, nx=3, ny=3, nz=3, K=1.0 * unit.kilojoules_per_mole/unit.nanometer**2):

        # Default parameters
        mass_Ar     = 39.9 * unit.amu
        q_Ar        = 0.0 * unit.elementary_charge
        sigma_Ar    = 3.350 * unit.angstrom
        epsilon_Ar  = 0.001603 * unit.kilojoule_per_mole

        scaleStepSizeX = 1.0
        scaleStepSizeY = 1.0
        scaleStepSizeZ = 1.0

        # Determine total number of atoms.
        natoms = nx * ny * nz

        # Create an empty system object.
        system = mm.System()

        # Create a NonbondedForce object with no cutoff.
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

        positions = unit.Quantity(np.zeros([natoms,3],np.float32), unit.angstrom)

        atom_index = 0
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    system.addParticle(mass_Ar)
                    nb.addParticle(q_Ar, sigma_Ar, epsilon_Ar)
                    x = sigma_Ar*scaleStepSizeX*(ii - nx/2.0)
                    y = sigma_Ar*scaleStepSizeY*(jj - ny/2.0)
                    z = sigma_Ar*scaleStepSizeZ*(kk - nz/2.0)

                    positions[atom_index,0] = x
                    positions[atom_index,1] = y
                    positions[atom_index,2] = z
                    atom_index += 1

        # Add the nonbonded force.
        system.addForce(nb)

        # Add a restrining potential centered at the origin.
        energy_expression = '(K/2.0) * (x^2 + y^2 + z^2);'
        energy_expression += 'K = %f;' % (K / (unit.kilojoules_per_mole / unit.nanometers**2)) # in OpenMM units
        force = mm.CustomExternalForce(energy_expression)
        #force.addGlobalParameter('K', K)
        for particle_index in range(natoms):
            force.addParticle(particle_index, [])
        system.addForce(force)

        self.system, self.positions = system, positions

#=============================================================================================
# Lennard-Jones fluid
#=============================================================================================

class LennardJonesFluid(TestSystem):
    """Create a periodic rectilinear grid of Lennard-Jones particles.
    Parameters for argon are used by default. Cutoff is set to 3 sigma by default.

    Parameters
    ----------
    nx : int, optional, default=6
        number of particles in the x direction
    ny : int, optional, default=6
        number of particles in the y direction
    nz : int, optional, default=6
        number of particles in the z direction
    mass : simtk.unit.Quantity, optional, default=39.9 * unit.amu
        mass of each particle.
    sigma : simtk.unit.Quantity, optional, default=3.4 * unit.angstrom
        Lennard-Jones sigma parameter
    epsilon : simtk.unit.Quantity, optional, default=0.238 * unit.kilocalories_per_mole
        Lennard-Jones well depth
    cutoff : simtk.unit.Quantity, optional, default=None
        Cutoff for nonbonded interactions.  If None, defaults to 2.5 * sigma
    switch : simtk.unit.Quantity, optional, default=1.0 * unit.kilojoules_per_mole/unit.nanometer**2
        if specified, the switching function will be turned on at this distance (default: None)
    switch_width : simtk.unit.Quantity with units compatible with angstroms, optional, default=0.2*unit.angstroms
        switching function is turned on at cutoff - switch_width
    dispersion_correction : bool, optional, default=True
        if True, will use analytical dispersion correction (if not using switching function)

    Examples
    --------

    Create default-size Lennard-Jones fluid.

    >>> fluid = LennardJonesFluid()
    >>> system, positions = fluid.system, fluid.positions

    Create a larger 10x8x5 box of Lennard-Jones particles.

    >>> fluid = LennardJonesFluid(nx=10, ny=8, nz=5)
    >>> system, positions = fluid.system, fluid.positions

    Create Lennard-Jones fluid using switched particle interactions (switched off betwee 7 and 9 A) and more particles.

    >>> fluid = LennardJonesFluid(nx=10, ny=10, nz=10, switch=True, switch_width=7.0*unit.angstroms, cutoff=9.0*unit.angstroms)
    >>> system, positions = fluid.system, fluid.positions
    """

    def __init__(self, nx=6, ny=6, nz=6,
        mass=39.9 * unit.amu, # argon
        sigma=3.4 * unit.angstrom, # argon,
        epsilon=0.238 * unit.kilocalories_per_mole, # argon,
        cutoff=None,
        switch=False,
        switch_width=0.2*unit.angstrom,
        dispersion_correction=True):

        if cutoff is None:
            cutoff = 2.5 * sigma

        charge        = 0.0 * unit.elementary_charge

        scaleStepSizeX = 1.0
        scaleStepSizeY = 1.0
        scaleStepSizeZ = 1.0

        # Determine total number of atoms.
        natoms = nx * ny * nz

        # Create an empty system object.
        system = mm.System()

        # Set up periodic nonbonded interactions with a cutoff.
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        nb.setCutoffDistance(cutoff)
        nb.setUseDispersionCorrection(dispersion_correction)
        nb.setUseSwitchingFunction(switch)
        nb.setSwitchingDistance(cutoff-switch_width)

        positions = unit.Quantity(np.zeros([natoms,3],np.float32), unit.angstrom)

        maxX = 0.0 * unit.angstrom
        maxY = 0.0 * unit.angstrom
        maxZ = 0.0 * unit.angstrom

        atom_index = 0
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    system.addParticle(mass)
                    nb.addParticle(charge, sigma, epsilon)
                    x = sigma*scaleStepSizeX*ii
                    y = sigma*scaleStepSizeY*jj
                    z = sigma*scaleStepSizeZ*kk

                    positions[atom_index,0] = x
                    positions[atom_index,1] = y
                    positions[atom_index,2] = z
                    atom_index += 1

                    # Wrap positions as needed.
                    if x>maxX: maxX = x
                    if y>maxY: maxY = y
                    if z>maxZ: maxZ = z

        # Set periodic box vectors.
        x = maxX+2*sigma*scaleStepSizeX
        y = maxY+2*sigma*scaleStepSizeY
        z = maxZ+2*sigma*scaleStepSizeZ

        a = unit.Quantity((x,                0*unit.angstrom, 0*unit.angstrom))
        b = unit.Quantity((0*unit.angstrom,                y, 0*unit.angstrom))
        c = unit.Quantity((0*unit.angstrom, 0*unit.angstrom, z))
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Add the nonbonded force.
        system.addForce(nb)

        self.system, self.positions = system, positions

#=============================================================================================
# Custom Lennard-Jones fluid mixture of NonbondedForce and CustomNonbondedForce
#=============================================================================================

class CustomLennardJonesFluidMixture(TestSystem):
    """Create a periodic rectilinear grid of Lennard-Jones particled, but implemented via CustomBondForce and NonbondedForce.
    Parameters for argon are used by default. Cutoff is set to 3 sigma by default.

    Parameters
    ----------
    nx : int, optional, default=6
        number of particles in the x direction
    ny : int, optional, default=6
        number of particles in the y direction
    nz : int, optional, default=6
        number of particles in the z direction
    mass : simtk.unit.Quantity, optional, default=39.9 * unit.amu
        mass of each particle.
    sigma : simtk.unit.Quantity, optional, default=3.4 * unit.angstrom
        Lennard-Jones sigma parameter
    epsilon : simtk.unit.Quantity, optional, default=0.238 * unit.kilocalories_per_mole
        Lennard-Jones well depth
    cutoff : simtk.unit.Quantity, optional, default=None
        Cutoff for nonbonded interactions.  If None, defaults to 2.5 * sigma
    switch : simtk.unit.Quantity, optional, default=1.0 * unit.kilojoules_per_mole/unit.nanometer**2
        if specified, the switching function will be turned on at this distance (default: None)
    switch_width : simtk.unit.Quantity with units compatible with angstroms, optional, default=0.2*unit.angstroms
        switching function is turned on at cutoff - switch_width
    dispersion_correction : bool, optional, default=True
        if True, will use analytical dispersion correction (if not using switching function)

    Notes
    -----

    No analytical dispersion correction is included here.

    Examples
    --------

    Create default-size Lennard-Jones fluid.

    >>> fluid = CustomLennardJonesFluidMixture()
    >>> system, positions = fluid.system, fluid.positions

    Create a larger 10x8x5 box of Lennard-Jones particles.

    >>> fluid = CustomLennardJonesFluidMixture(nx=10, ny=8, nz=5)
    >>> system, positions = fluid.system, fluid.positions

    Create Lennard-Jones fluid using switched particle interactions (switched off betwee 7 and 9 A) and more particles.

    >>> fluid = CustomLennardJonesFluidMixture(nx=10, ny=10, nz=10, switch=True, switch_width=7.0*unit.angstroms, cutoff=9.0*unit.angstroms)
    >>> system, positions = fluid.system, fluid.positions
    """

    def __init__(self, nx=6, ny=6, nz=6,
        mass=39.9 * unit.amu, # argon
        sigma=3.4 * unit.angstrom, # argon,
        epsilon=0.238 * unit.kilocalories_per_mole, # argon,
        cutoff=None,
        switch=False,
        switch_width=0.2*unit.angstroms,
        dispersion_correction=True):

        if cutoff is None:
            cutoff = 2.5 * sigma

        charge        = 0.0 * unit.elementary_charge
        scaleStepSizeX = 1.0
        scaleStepSizeY = 1.0
        scaleStepSizeZ = 1.0

        charge = 0.0 * unit.elementary_charge

        # Determine total number of atoms.
        natoms = nx * ny * nz

        # determine number of atoms that will be treated by CustomNonbondedForce
        ncustom = int(natoms/2)

        # Create an empty system object.
        system = mm.System()

        # Set up periodic nonbonded interactions with a cutoff.
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        nb.setCutoffDistance(cutoff)
        nb.setUseDispersionCorrection(dispersion_correction)
        nb.setUseSwitchingFunction(switch)
        nb.setSwitchingDistance(cutoff-switch_width)

        # Set up periodic nonbonded interactions with a cutoff.
        if switch:
            energy_expression = "LJ * S;"
            energy_expression += "LJ = 4*epsilon*((sigma/r)^12 - (sigma/r)^6);"
            energy_expression += "sigma = 0.5*(sigma1+sigma2);"
            energy_expression += "epsilon = sqrt(epsilon1*epsilon2);"
            energy_expression += "S = (cutoff^2 - r^2)^2 * (cutoff^2 + 2*r^2 - 3*switch^2) / (cutoff^2 - switch^2)^3;"
            cnb = mm.CustomNonbondedForce(energy_expression)
            cnb.addGlobalParameter('switch', switch)
            cnb.addGlobalParameter('cutoff', cutoff)
            cnb.addPerParticleParameter('q')
            cnb.addPerParticleParameter('sigma')
            cnb.addPerParticleParameter('epsilon')
            cnb.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            cnb.setCutoffDistance(cutoff)
        else:
            energy_expression = "4*epsilon*((sigma/r)^12 - (sigma/r)^6);"
            cnb = mm.CustomNonbondedForce(energy_expression)
            cnb.addGlobalParameter('sigma', sigma)
            cnb.addGlobalParameter('epsilon', epsilon)
            cnb.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            cnb.setCutoffDistance(cutoff)
        # Only add interactions between first atom and rest.
        #atomset1 = range(0, 1)
        #atomset2 = range(1, natoms)
        #cnb.addInteractionGroup(atomset1, atomset2)

        positions = unit.Quantity(np.zeros([natoms,3],np.float32), unit.angstrom)

        maxX = 0.0 * unit.angstrom
        maxY = 0.0 * unit.angstrom
        maxZ = 0.0 * unit.angstrom

        atom_index = 0
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    system.addParticle(mass)
                    if (atom_index < ncustom):
                        cnb.addParticle([charge, sigma, epsilon])
                        nb.addParticle(0.0*charge, sigma, 0.0*epsilon)
                    else:
                        cnb.addParticle([0.0*charge, sigma, 0.0*epsilon])
                        nb.addParticle(charge, sigma, epsilon)
                    x = sigma*scaleStepSizeX*ii
                    y = sigma*scaleStepSizeY*jj
                    z = sigma*scaleStepSizeZ*kk

                    positions[atom_index,0] = x
                    positions[atom_index,1] = y
                    positions[atom_index,2] = z
                    atom_index += 1

                    # Wrap positions as needed.
                    if x>maxX: maxX = x
                    if y>maxY: maxY = y
                    if z>maxZ: maxZ = z

        # Set periodic box vectors.
        x = maxX+2*sigma*scaleStepSizeX
        y = maxY+2*sigma*scaleStepSizeY
        z = maxZ+2*sigma*scaleStepSizeZ

        a = unit.Quantity((x,                0*unit.angstrom, 0*unit.angstrom))
        b = unit.Quantity((0*unit.angstrom,                y, 0*unit.angstrom))
        c = unit.Quantity((0*unit.angstrom, 0*unit.angstrom, z))
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Add the nonbonded forces.
        system.addForce(nb)
        system.addForce(cnb)

        # Add long-range correction.
        if switch:
            # TODO
            pass
        else:
            volume = x*y*z
            density = natoms / volume
            per_particle_dispersion_energy = -(8./3.)*math.pi*epsilon*(sigma**6)/(cutoff**3)*density  # attraction
            per_particle_dispersion_energy += (8./9.)*math.pi*epsilon*(sigma**12)/(cutoff**9)*density  # repulsion
            energy_expression = "%f" % (per_particle_dispersion_energy / unit.kilojoules_per_mole)
            force = mm.CustomExternalForce(energy_expression)
            for i in range(natoms):
                force.addParticle(i, [])
            system.addForce(force)

        self.system, self.positions = system, positions


#=============================================================================================
# WCA Fluid
#=============================================================================================

class WCAFluid(TestSystem):

    def __init__(self, N=216, density=0.96, mass=39.9*unit.amu, epsilon=120.0*unit.kelvin*kB, sigma=3.4*unit.angstrom):
        """
        Create a Weeks-Chandler-Andersen system.

        Parameters:
        -----------
        N : int, optional, default = 216
            Number of particles.
        density : float, optional, default = 0.96
            Reduced density, N sigma^3 / V.
        mass : simtk.unit.Quantity with units compatible with angstrom, optional, default=39.9 amu
            Particle mass.
        epsilon : simtk.unit.Quantity with units compatible with kilocalories_per_mole, optional, default=120K*kB
            WCA well depth.
        sigma : simtk.unit.Quantity
            WCA sigma.

        """
        # Create system
        system = mm.System()

        # Compute total system volume.
        volume = N / density

        # Make system cubic in dimension.
        length = volume**(1.0/3.0)
        # TODO: Can we change this to use tuples or 3x3 array?
        a = unit.Quantity(numpy.array([1.0, 0.0, 0.0], numpy.float32), unit.nanometer) * length/unit.nanometer
        b = unit.Quantity(numpy.array([0.0, 1.0, 0.0], numpy.float32), unit.nanometer) * length/unit.nanometer
        c = unit.Quantity(numpy.array([0.0, 0.0, 1.0], numpy.float32), unit.nanometer) * length/unit.nanometer
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Add particles to system.
        for n in range(N):
            system.addParticle(mass)

        # Create nonbonded force term implementing Kob-Andersen two-component Lennard-Jones interaction.
        energy_expression = '4.0*epsilon*((sigma/r)^12 - (sigma/r)^6) + epsilon'

        # Create force.
        force = mm.CustomNonbondedForce(energy_expression)

        # Set epsilon and sigma global parameters.
        force.addGlobalParameter('epsilon', epsilon)
        force.addGlobalParameter('sigma', sigma)

        # Add particles
        for n in range(N):
            force.addParticle([])

        # Set periodic boundary conditions with cutoff.
        force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        force.setCutoffDistance(sigma)

        # Add nonbonded force term to the system.
        system.addForce(force)

        # Create initial coordinates using random positions.
        coordinates = unit.Quantity(numpy.random.rand(N,3), unit.nanometer) * (length / unit.nanometer)

        # Store system.
        self.system, self.positions = system, coordinates

#=============================================================================================
# Ideal gas
#=============================================================================================

class IdealGas(TestSystem):
    """Create an 'ideal gas' of noninteracting particles in a periodic box.

    Parameters
    ----------
    nparticles : int, optional, default=216
        number of particles
    mass : int, optional, default=39.9 * unit.amu
    temperature : int, optional, default=298.0 * unit.kelvin
    pressure : int, optional, default=1.0 * unit.atmosphere
    volume : None
        if None, defaults to (nparticles * temperature * unit.BOLTZMANN_CONSTANT_kB / pressure).in_units_of(unit.nanometers**3)

    Examples
    --------

    Create an ideal gas system.

    >>> gas = IdealGas()
    >>> system, positions = gas.system, gas.positions

    Create a smaller ideal gas system containing 64 particles.

    >>> gas = IdealGas(nparticles=64)
    >>> system, positions = gas.system, gas.positions

    """

    def __init__(self, nparticles=216, mass=39.9 * unit.amu, temperature=298.0 * unit.kelvin, pressure=1.0 * unit.atmosphere, volume=None):

        if volume is None: 
            volume = (nparticles * temperature * unit.BOLTZMANN_CONSTANT_kB / pressure).in_units_of(unit.nanometers**3)

        charge   = 0.0 * unit.elementary_charge
        sigma    = 3.350 * unit.angstrom # argon LJ 
        epsilon  = 0.0 * unit.kilojoule_per_mole # zero interaction

        # Create an empty system object.
        system = mm.System()

        # Compute box size.
        length = volume**(1.0/3.0)
        a = unit.Quantity((length,           0*unit.nanometer, 0*unit.nanometer))
        b = unit.Quantity((0*unit.nanometer,           length, 0*unit.nanometer))
        c = unit.Quantity((0*unit.nanometer, 0*unit.nanometer, length))
        system.setDefaultPeriodicBoxVectors(a, b, c)

        # Add particles.
        for index in range(nparticles):
            system.addParticle(mass)

        # Place particles at random positions within the box.
        # TODO: Use reproducible seed.
        # NOTE: This may not be thread-safe.

        state = np.random.get_state()
        np.random.seed(0)
        positions = unit.Quantity((length/unit.nanometer) * np.random.rand(nparticles,3), unit.nanometer)
        np.random.set_state(state)

        self.system, self.positions = system, positions
        self.ndof = 3 * nparticles

    def get_potential_expectation(self, state):
        """Return the expectation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return 0.0 * unit.kilojoules_per_mole
        
    def get_potential_standard_deviation(self, state):
        """Return the standard deviation of the potential energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------
        
        potential_stddev : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            potential energy standard deviation if implemented, or else None
        
        """

        return 0.0 * unit.kilojoules_per_mole

    def get_kinetic_expectation(self, state):
        """Return the expectation of the kinetic energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        potential_mean : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            The expectation of the potential energy.
        
        """

        return (3./2.) * kB * state.temperature 
        
    def get_kinetic_standard_deviation(self, state):
        """Return the standard deviation of the kinetic energy, computed analytically or numerically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------
        
        potential_stddev : simtk.unit.Quantity compatible with simtk.unit.kilojoules_per_mole
            potential energy standard deviation if implemented, or else None
        
        """

        return (3./2.) * kB * state.temperature 

    def get_volume_expectation(self, state):
        """Return the expectation of the volume, computed analytically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature and pressure defined
            The thermodynamic state at which the property is to be computed.
        
        Returns
        -------
        
        volume_mean : simtk.unit.Quantity compatible with simtk.unit.nanometers**3
            The expectation of the volume at equilibrium.
        
        Notes
        -----
        
        The true mean volume is used, rather than the large-N limit.

        """
        
        if not state.pressure:
            box_vectors = self.system.getDefaultPeriodicBoxVectors()
            volume = box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2] 
            return volume

        N = self._system.getNumParticles()
        return ((N+1) * unit.BOLTZMANN_CONSTANT_kB * state.temperature / state.pressure).in_units_of(unit.nanometers**3)
        
    def get_volume_standard_deviation(self, state):
        """Return the standard deviation of the volume, computed analytically.

        Arguments
        ---------
        
        state : ThermodynamicState with temperature and pressure defined
            The thermodynamic state at which the property is to be computed.

        Returns
        -------
        
        volume_stddev : simtk.unit.Quantity compatible with simtk.unit.nanometers**3
            The standard deviation of the volume at equilibrium.
        
        Notes
        -----
        
        The true mean volume is used, rather than the large-N limit.

        """
        
        if not state.pressure:
            return 0.0 * unit.nanometers**3

        N = self._system.getNumParticles()
        return (numpy.sqrt(N+1) * unit.BOLTZMANN_CONSTANT_kB * state.temperature / state.pressure).in_units_of(unit.nanometers**3)
    
#=============================================================================================
# Water box
#=============================================================================================

class WaterBox(TestSystem):
   """
   Create a water box test system.

   Examples
   --------
   
   Create a default (TIP3P) waterbox.

   >>> waterbox = WaterBox()

   Control the cutoff.
   
   >>> waterbox = WaterBox(box_edge=3.0*unit.nanometers, cutoff=1.0*unit.nanometers)

   Use a different water model.

   >>> waterbox = WaterBox(model='tip4pew')

   Don't use constraints.

   >>> waterbox = WaterBox(constrained=False)

   """

   def __init__(self, box_edge=2.5*unit.nanometers, cutoff=0.9*unit.nanometers, model='tip3p', switch=True, switch_width=0.5*unit.angstroms, constrained=True, dispersion_correction=True, nonbondedMethod=app.PME):
       """
       Create a water box test system.
       
       Parameters
       ----------
       
       box_edge : simtk.unit.Quantity with units compatible with nanometers, optional, default = 2.5 nm
          Edge length for cubic box [should be greater than 2*cutoff]
       cutoff : simtk.unit.Quantity with units compatible with nanometers, optional, default = 0.9 nm
          Nonbonded cutoff
       model : str, optional, default = 'tip3p'
          The name of the water model to use ['tip3p', 'tip4p', 'tip4pew', 'tip5p', 'spce']
       switch : bool, optional, default = True
          Turns the Lennard-Jones switching function on or off.
       switch_width : simtk.unit.Quantity with units compatible with nanometers, optional, default = 0.5 A
          Sets the width of the switch function for Lennard-Jones.
       constrained : bool, optional, default=True
          Sets whether water geometry should be constrained (rigid water implemented via SETTLE) or flexible.
       dispersion_correction : bool, optional, default=True
          Sets whether the long-range dispersion correction should be used.
       nonbondedMethod : simtk.openmm.app nonbonded method, optional, default=app.PME
          Sets the nonbonded method to use for the water box (one of app.CutoffPeriodic, app.Ewald, app.PME).

       Examples
       --------
       
       Create a default waterbox.
       
       >>> waterbox = WaterBox()
       >>> [system, positions] = [waterbox.system, waterbox.positions]
       
       Use reaction-field electrostatics instead.

       >>> waterbox = WaterBox(nonbondedMethod=app.CutoffPeriodic)

       Control the cutoff.
       
       >>> waterbox = WaterBox(box_edge=3.0*unit.nanometers, cutoff=1.0*unit.nanometers)
       
       Use a different water model.
       
       >>> waterbox = WaterBox(model='spce')

       Use a five-site water model.
       
       >>> waterbox = WaterBox(model='tip5p')

       Turn off the switch function.

       >>> waterbox = WaterBox(switch=False)

       Set the switch width.

       >>> waterbox = WaterBox(switch=True, switch_width=0.8*unit.angstroms)

       Turn of long-range dispersion correction.

       >>> waterbox = WaterBox(dispersion_correction=False)

       """

       import simtk.openmm.app as app
       
       supported_models = ['tip3p', 'tip4pew', 'tip5p', 'spce']
       if model not in supported_models:
           raise Exception("Specified water model '%s' is not in list of supported models: %s" % (model, str(supported_models)))

       # Load forcefield for solvent model.
       ff =  app.ForceField(model + '.xml')
       
       # Create empty topology and coordinates.
       top = app.Topology()
       pos = unit.Quantity((), unit.angstroms)
       
       # Create new Modeller instance.
       m = app.Modeller(top, pos)
       
       # Add solvent to specified box dimensions.
       boxSize = unit.Quantity(numpy.ones([3]) * box_edge/box_edge.unit, box_edge.unit)
       m.addSolvent(ff, boxSize=boxSize, model=model)

       # Get new topology and coordinates.
       newtop = m.getTopology()
       newpos = m.getPositions()
   
       # Convert positions to numpy.
       positions = unit.Quantity(numpy.array(newpos / newpos.unit), newpos.unit)
   
       # Create OpenMM System.
       system = ff.createSystem(newtop, nonbondedMethod=nonbondedMethod, nonbondedCutoff=cutoff, constraints=None, rigidWater=constrained, removeCMMotion=False)

       # Set switching function and dispersion correction.
       forces = { system.getForce(index).__class__.__name__ : system.getForce(index) for index in range(system.getNumForces()) }
       forces['NonbondedForce'].setUseSwitchingFunction(switch)
       forces['NonbondedForce'].setSwitchingDistance(cutoff - switch_width)
       forces['NonbondedForce'].setUseDispersionCorrection(dispersion_correction)

       self.ndof = 3*system.getNumParticles() - 3*constrained
       self.system, self.positions = system, positions

class FlexibleWaterBox(WaterBox):
   """
   Flexible water box.

   """

   def __init__(self, *args, **kwargs):
       """
       Create a flexible water box.
       
       Parameters are inherited from WaterBox (except for 'constrained').
              
       Examples
       --------
       
       Create a default flexible waterbox.
       
       >>> waterbox = FlexibleWaterBox()
       >>> [system, positions] = [waterbox.system, waterbox.positions]
       
       """
       super(FlexibleWaterBox, self).__init__(constrained=False, *args, **kwargs)

class FourSiteWaterBox(WaterBox):
   """
   Four-site water box (TIP4P-Ew).

   """

   def __init__(self, *args, **kwargs):
       """
       Create a water box test systemm using a four-site water model (TIP4P-Ew).
              
       Parameters are inherited from WaterBox (except for 'model').

       Examples
       --------
       
       Create a default waterbox.
       
       >>> waterbox = FourSiteWaterBox()
       >>> [system, positions] = [waterbox.system, waterbox.positions]
       
       Control the cutoff.
       
       >>> waterbox = FourSiteWaterBox(box_edge=3.0*unit.nanometers, cutoff=1.0*unit.nanometers)
       
       """
       super(FourSiteWaterBox, self).__init__(model='tip4pew', *args, **kwargs)

class FiveSiteWaterBox(WaterBox):
   """
   Five-site water box (TIP5P).

   """

   def __init__(self, *args, **kwargs):
       """
       Create a water box test systemm using a five-site water model (TIP5P).
       
       Parameters are inherited from WaterBox (except for 'model').
       
       Examples
       --------
       
       Create a default waterbox.
       
       >>> waterbox = FiveSiteWaterBox()
       >>> [system, positions] = [waterbox.system, waterbox.positions]
       
       Control the cutoff.
       
       >>> waterbox = FiveSiteWaterBox(box_edge=3.0*unit.nanometers, cutoff=1.0*unit.nanometers)
       
       """
       super(FiveSiteWaterBox, self).__init__(model='tip5p', *args, **kwargs)

class DischargedWaterBox(WaterBox):
   """
   Water box test system with zeroed charges.

   """

   def __init__(self, *args, **kwargs):
       """
       Create a water box test systemm using a four-site water model (TIP4P-Ew).
       
       Parameters are inherited from WaterBox.
       
       Examples
       --------
       
       Create a default waterbox.
       
       >>> waterbox = DischargedWaterBox()
       >>> [system, positions] = [waterbox.system, waterbox.positions]
       
       Control the cutoff.
       
       >>> waterbox = DischargedWaterBox(box_edge=3.0*unit.nanometers, cutoff=1.0*unit.nanometers)
       
       """
       super(DischargedWaterBox, self).__init__(*args, **kwargs)

       # Zero charges.
       system = self.system
       forces = { system.getForce(index).__class__.__name__ : system.getForce(index) for index in range(system.getNumForces()) }
       force = forces['NonbondedForce']
       for index in range(force.getNumParticles()):
           [charge, sigma, epsilon] = force.getParticleParameters(index)
           force.setParticleParameters(index, 0*charge, sigma, epsilon)
       for index in range(force.getNumExceptions()):
           [particle1, particle2, chargeProd, sigma, epsilon] = force.getExceptionParameters(index)
           force.setExceptionParameters(index, particle1, particle2, 0*chargeProd, sigma, epsilon)
           
       return

class DischargedWaterBoxHsites(WaterBox):
   """
   Water box test system with zeroed charges and Lennard-Jones sites on hydrogens.

   """

   def __init__(self, *args, **kwargs):
       """
       Create a water box with zeroed charges and Lennard-Jones sites on hydrogens.
       
       Parameters are inherited from WaterBox.
       
       Examples
       --------
       
       Create a default waterbox.
       
       >>> waterbox = DischargedWaterBox()
       >>> [system, positions] = [waterbox.system, waterbox.positions]
       
       Control the cutoff.
       
       >>> waterbox = DischargedWaterBox(box_edge=3.0*unit.nanometers, cutoff=1.0*unit.nanometers)
       
       """
       super(DischargedWaterBoxHsites, self).__init__(*args, **kwargs)

       # Zero charges.
       system = self.system
       forces = { system.getForce(index).__class__.__name__ : system.getForce(index) for index in range(system.getNumForces()) }
       force = forces['NonbondedForce']
       for index in range(force.getNumParticles()):
           [charge, sigma, epsilon] = force.getParticleParameters(index)
           charge *= 0
           if epsilon == 0.0 * unit.kilojoules_per_mole:
               # Add LJ site to hydrogens.
               epsilon = 0.0157 * unit.kilojoules_per_mole
               sigma = 0.06 * unit.angstroms
           force.setParticleParameters(index, charge, sigma, epsilon)
       for index in range(force.getNumExceptions()):
           [particle1, particle2, chargeProd, sigma, epsilon] = force.getExceptionParameters(index)
           chargeProd *= 0
           epsilon *= 0
           force.setExceptionParameters(index, particle1, particle2, chargeProd, sigma, epsilon)
           
       return

#=============================================================================================
# Alanine dipeptide in vacuum.
#=============================================================================================

class AlanineDipeptideVacuum(TestSystem):
    """Alanine dipeptide ff96 in vacuum.
    
    Parameters
    ----------
    constraints : optional, default=simtk.openmm.app.HBonds
    
    Examples
    --------
    
    Create alanine dipeptide with constraints on bonds to hydrogen
    >>> alanine = AlanineDipeptideVacuum()
    >>> (system, positions) = alanine.system, alanine.positions
    """

    def __init__(self, constraints=app.HBonds):

        prmtop_filename = get_data_filename("data/alanine-dipeptide-gbsa/alanine-dipeptide.prmtop")
        crd_filename = get_data_filename("data/alanine-dipeptide-gbsa/alanine-dipeptide.crd")

        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(implicitSolvent=None, constraints=constraints, nonbondedCutoff=None)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename)
        positions = inpcrd.getPositions(asNumpy=True)

        self.system, self.positions = system, positions

#=============================================================================================
# Alanine dipeptide in implicit solvent.
#=============================================================================================

class AlanineDipeptideImplicit(TestSystem):
    """Alanine dipeptide ff96 in OBC GBSA implicit solvent.
    
    Parameters
    ----------
    constraints : optional, default=simtk.openmm.app.HBonds
    
    Examples
    --------
    
    Create alanine dipeptide with constraints on bonds to hydrogen
    >>> alanine = AlanineDipeptideImplicit()
    >>> (system, positions) = alanine.system, alanine.positions
    """

    def __init__(self, constraints=app.HBonds):

        prmtop_filename = get_data_filename("data/alanine-dipeptide-gbsa/alanine-dipeptide.prmtop")
        crd_filename = get_data_filename("data/alanine-dipeptide-gbsa/alanine-dipeptide.crd")        

        # Initialize system.
        
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(implicitSolvent=app.OBC1, constraints=constraints, nonbondedCutoff=None)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename)
        positions = inpcrd.getPositions(asNumpy=True)

        self.system, self.positions = system, positions

#=============================================================================================
# Alanine dipeptide in explicit solvent
#=============================================================================================

class AlanineDipeptideExplicit(TestSystem):
    """Alanine dipeptide ff96 in TIP3P explicit solvent..

    Parameters
    ----------
    constraints : optional, default=simtk.openmm.app.HBonds
    rigid_water : bool, optional, default=True
    nonbondedCutoff : Quantity, optional, default=9.0 * unit.angstroms
    use_dispersion_correction : bool, optional, default=True
        If True, the long-range disperson correction will be used.
    nonbondedMethod : simtk.openmm.app nonbonded method, optional, default=app.PME
       Sets the nonbonded method to use for the water box (one of app.CutoffPeriodic, app.Ewald, app.PME).

    Examples
    --------

    >>> alanine = AlanineDipeptideExplicit()
    >>> (system, positions) = alanine.system, alanine.positions
    """

    def __init__(self, constraints=app.HBonds, rigid_water=True, nonbondedCutoff=9.0 * unit.angstroms, use_dispersion_correction=True, nonbondedMethod=app.PME):

        prmtop_filename = get_data_filename("data/alanine-dipeptide-explicit/alanine-dipeptide.prmtop")
        crd_filename = get_data_filename("data/alanine-dipeptide-explicit/alanine-dipeptide.crd")

        # Initialize system.
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(constraints=constraints, nonbondedMethod=nonbondedMethod, rigidWater=rigid_water, nonbondedCutoff=0.9*unit.nanometer)

        # Set dispersion correction use.
        forces = { system.getForce(index).__class__.__name__ : system.getForce(index) for index in range(system.getNumForces()) }
        forces['NonbondedForce'].setUseDispersionCorrection(use_dispersion_correction)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename, loadBoxVectors=True)
        positions = inpcrd.getPositions(asNumpy=True)

        # Set box vectors.
        box_vectors = inpcrd.getBoxVectors(asNumpy=True)
        system.setDefaultPeriodicBoxVectors(box_vectors[0], box_vectors[1], box_vectors[2])

        self.system, self.positions = system, positions

#=============================================================================================
# T4 lysozyme L99A mutant with p-xylene ligand.
#=============================================================================================

class LysozymeImplicit(TestSystem):
    """T4 lysozyme L99A (AMBER ff96) with p-xylene ligand (GAFF + AM1-BCC) in implicit OBC GBSA solvent.

    Parameters
    ----------
    flexibleConstraints : bool, optional, default=True
    constraints : simtk.openmm.app constraints (None, HBonds, HAngles, AllBonds)
       constraints to be imposed

    Examples
    --------

    >>> lysozyme = LysozymeImplicit()
    >>> (system, positions) = lysozyme.system, lysozyme.positions
    """

    def __init__(self, flexibleConstraints=True, constraints=app.HBonds, implicitSolvent=app.OBC1):

        prmtop_filename = get_data_filename("data/T4-lysozyme-L99A-implicit/complex.prmtop")
        crd_filename = get_data_filename("data/T4-lysozyme-L99A-implicit/complex.crd")

        # Initialize system.
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(implicitSolvent=app.OBC1, constraints=app.HBonds, nonbondedCutoff=None)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename)
        positions = inpcrd.getPositions(asNumpy=True)

        self.system, self.positions = system, positions


class SrcImplicit(TestSystem):
    """Src kinase in implicit AMBER 99sb-ildn with OBC GBSA solvent.

    Examples
    --------
    >>> src = SrcImplicit()
    >>> system, positions = src.system, src.positions
    """

    def __init__(self):

        pdb_filename = get_data_filename("data/src-implicit/implicit-refined.pdb")
        pdbfile = app.PDBFile(pdb_filename)

        # Construct system.
        forcefields_to_use = ['amber99sbildn.xml', 'amber99_obc.xml'] # list of forcefields to use in parameterization
        forcefield = app.ForceField(*forcefields_to_use)
        system = forcefield.createSystem(pdbfile.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)

        # Get positions.
        positions = pdbfile.getPositions()

        self.system, self.positions = system, positions

#=============================================================================================
# Src kinase in explicit solvent.
#=============================================================================================

class SrcExplicit(TestSystem):
    """Src kinase (AMBER 99sb-ildn) in explicit TIP3P solvent.

    Parameters
    ----------
    nonbondedMethod : simtk.openmm.app nonbonded method, optional, default=app.PME
       Sets the nonbonded method to use for the water box (CutoffPeriodic, app.Ewald, app.PME).

    Examples
    --------
    >>> src = SrcExplicit()
    >>> system, positions = src.system, src.positions

    """
    def __init__(self, nonbondedMethod=app.PME):

        system_xml_filename = get_data_filename("data/src-explicit/system.xml")
        state_xml_filename = get_data_filename("data/src-explicit/state.xml")

        # Read system.
        infile = open(system_xml_filename, 'r')
        system = mm.XmlSerializer.deserialize(infile.read())
        infile.close()

        # Read state.
        infile = open(state_xml_filename, 'r')
        serialized_state = mm.XmlSerializer.deserialize(infile.read())
        infile.close()

        # Select nonbonded method.
        forces = { system.getForce(index).__class__.__name__ : system.getForce(index) for index in range(system.getNumForces()) }
        from simtk.openmm import NonbondedForce
        methodMap = {app.NoCutoff:NonbondedForce.NoCutoff,
                     app.CutoffNonPeriodic:NonbondedForce.CutoffNonPeriodic,
                     app.CutoffPeriodic:NonbondedForce.CutoffPeriodic,
                     app.Ewald:NonbondedForce.Ewald,
                     app.PME:NonbondedForce.PME}
        forces['NonbondedForce'].setNonbondedMethod(methodMap[nonbondedMethod])

        # Get positions and set periodic box vectors.
        positions = serialized_state.getPositions()
        box_vectors = serialized_state.getPeriodicBoxVectors()
        system.setDefaultPeriodicBoxVectors(*box_vectors)

        self.system, self.positions = system, positions

#=============================================================================================
# Methanol box.
#=============================================================================================

class MethanolBox(TestSystem):
    """Methanol box.

    Parameters
    ----------
    flexibleConstraints : bool, optional, default=True
    shake : string, optional, default="h-bonds"
    nonbondedCutoff : Quantity, optional, default=7.0 * unit.angstroms
    nonbondedMethod : simtk.openmm.app nonbonded method, optional, default=app.PME
       Sets the nonbonded method to use for the water box (one of app.CutoffPeriodic, app.Ewald, app.PME).

    Examples
    --------

    >>> methanol_box = MethanolBox()
    >>> system, positions = methanol_box.system, methanol_box.positions
    """

    def __init__(self, flexibleConstraints=True, constraints=app.HBonds, nonbondedCutoff=7.0 * unit.angstroms, nonbondedMethod=app.CutoffPeriodic):

        system_name = 'methanol-box'
        prmtop_filename = get_data_filename("data/%s/%s.prmtop" % (system_name, system_name))
        crd_filename = get_data_filename("data/%s/%s.crd" % (system_name, system_name))

        # Initialize system.
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        system = prmtop.createSystem(constraints=constraints, nonbondedMethod=nonbondedMethod, rigidWater=True, nonbondedCutoff=0.9*unit.nanometer)

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename, loadBoxVectors=True)
        positions = inpcrd.getPositions(asNumpy=True)

        # Set box vectors.
        box_vectors = inpcrd.getBoxVectors(asNumpy=True)
        system.setDefaultPeriodicBoxVectors(box_vectors[0], box_vectors[1], box_vectors[2])

        self.system, self.positions = system, positions

#=============================================================================================
# Molecular ideal gas (methanol box).
#=============================================================================================

class MolecularIdealGas(TestSystem):
    """Molecular ideal gas (methanol box).

    Parameters
    ----------
    flexibleConstraints : bool, optional, default=True
    shake : string, optional, default=None
    nonbondedCutoff : Quantity, optional, default=7.0 * unit.angstroms
    nonbondedMethod : simtk.openmm.app nonbonded method, optional, default=app.PME
       Sets the nonbonded method to use for the water box (one of app.CutoffPeriodic, app.Ewald, app.PME).

    Examples
    --------

    >>> methanol_box = MolecularIdealGas()
    >>> system, positions = methanol_box.system, methanol_box.positions
    """

    def __init__(self, flexibleConstraints=True, shake=None, nonbondedCutoff=7.0 * unit.angstroms, nonbondedMethod=app.CutoffPeriodic):

        system_name = 'methanol-box'
        prmtop_filename = get_data_filename("data/%s/%s.prmtop" % (system_name, system_name))
        crd_filename = get_data_filename("data/%s/%s.crd" % (system_name, system_name))

        # Initialize system.
        prmtop = app.AmberPrmtopFile(prmtop_filename)
        reference_system = prmtop.createSystem(constraints=app.HBonds, nonbondedMethod=nonbondedMethod, rigidWater=True, nonbondedCutoff=0.9*unit.nanometer)

        # Make a new system that contains no intermolecular interactions.
        system = mm.System()

        # Add atoms.
        for atom_index in range(reference_system.getNumParticles()):
            mass = reference_system.getParticleMass(atom_index)
            system.addParticle(mass)

        # Add constraints
        for constraint_index in range(reference_system.getNumConstraints()):
            [iatom, jatom, r0] = reference_system.getConstraintParameters(constraint_index)
            system.addConstraint(iatom, jatom, r0)

        # Copy only intramolecular forces.
        nforces = reference_system.getNumForces()
        for force_index in range(nforces):
            reference_force = reference_system.getForce(force_index)
            if isinstance(reference_force, mm.HarmonicBondForce):
                # HarmonicBondForce
                force = mm.HarmonicBondForce()
                for bond_index in range(reference_force.getNumBonds()):
                    [iatom, jatom, r0, K] = reference_force.getBondParameters(bond_index)
                    force.addBond(iatom, jatom, r0, K)
                system.addForce(force)
            elif isinstance(reference_force, mm.HarmonicAngleForce):
                # HarmonicAngleForce
                force = mm.HarmonicAngleForce()
                for angle_index in range(reference_force.getNumAngles()):
                    [iatom, jatom, katom, theta0, Ktheta] = reference_force.getAngleParameters(angle_index)
                    force.addAngle(iatom, jatom, katom, theta0, Ktheta)
                system.addForce(force)
            elif isinstance(reference_force, mm.PeriodicTorsionForce):
                # PeriodicTorsionForce
                force = mm.PeriodicTorsionForce()
                for torsion_index in range(reference_force.getNumTorsions()):
                    [particle1, particle2, particle3, particle4, periodicity, phase, k] = reference_force.getTorsionParameters(torsion_index)
                    force.addTorsion(particle1, particle2, particle3, particle4, periodicity, phase, k)
                system.addForce(force)
            else:
                # Don't add any other forces.
                pass

        # Read positions.
        inpcrd = app.AmberInpcrdFile(crd_filename, loadBoxVectors=True)
        positions = inpcrd.getPositions(asNumpy=True)

        # Set box vectors.
        box_vectors = inpcrd.getBoxVectors(asNumpy=True)
        system.setDefaultPeriodicBoxVectors(box_vectors[0], box_vectors[1], box_vectors[2])
        
        self.system, self.positions = system, positions

#=============================================================================================
# System of particles with CustomGBForce
#=============================================================================================

class CustomGBForceSystem(TestSystem):
    """A system of particles with a CustomGBForce.

    Notes
    -----

    This example comes from TestReferenceCustomGBForce.cpp from the OpenMM distribution.
    
    Examples
    --------
    
    >>> gb_system = CustomGBForceSystem()
    >>> system, positions = gb_system.system, gb_system.positions
    """

    def __init__(self):

        numMolecules = 70
        numParticles = numMolecules*2
        boxSize = 10.0 * unit.nanometers

        # Default parameters
        mass     = 39.9 * unit.amu
        sigma    = 3.350 * unit.angstrom
        epsilon  = 0.001603 * unit.kilojoule_per_mole
        cutoff   = 2.0 * unit.nanometers
        
        system = mm.System()
        for i in range(numParticles):
            system.addParticle(mass)

        system.setDefaultPeriodicBoxVectors(mm.Vec3(boxSize, 0.0, 0.0), mm.Vec3(0.0, boxSize, 0.0), mm.Vec3(0.0, 0.0, boxSize))

        # Create NonbondedForce.
        nonbonded = mm.NonbondedForce()    
        nonbonded.setNonbondedMethod(mm.NonbondedForce.CutoffPeriodic)
        nonbonded.setCutoffDistance(cutoff)

        # Create CustomGBForce.
        custom = mm.CustomGBForce()
        custom.setNonbondedMethod(mm.CustomGBForce.CutoffPeriodic)
        custom.setCutoffDistance(cutoff)
        
        custom.addPerParticleParameter("q")
        custom.addPerParticleParameter("radius")
        custom.addPerParticleParameter("scale")

        custom.addGlobalParameter("solventDielectric", 80.0)
        custom.addGlobalParameter("soluteDielectric", 1.0)
        custom.addComputedValue("I", "step(r+sr2-or1)*0.5*(1/L-1/U+0.25*(1/U^2-1/L^2)*(r-sr2*sr2/r)+0.5*log(L/U)/r+C);"
                                      "U=r+sr2;"
                                      "C=2*(1/or1-1/L)*step(sr2-r-or1);"
                                      "L=max(or1, D);"
                                      "D=abs(r-sr2);"
                                      "sr2 = scale2*or2;"
                                      "or1 = radius1-0.009; or2 = radius2-0.009", mm.CustomGBForce.ParticlePairNoExclusions);
        custom.addComputedValue("B", "1/(1/or-tanh(1*psi-0.8*psi^2+4.85*psi^3)/radius);"
                                      "psi=I*or; or=radius-0.009", mm.CustomGBForce.SingleParticle);
        custom.addEnergyTerm("28.3919551*(radius+0.14)^2*(radius/B)^6-0.5*138.935485*(1/soluteDielectric-1/solventDielectric)*q^2/B", mm.CustomGBForce.SingleParticle);
        custom.addEnergyTerm("-138.935485*(1/soluteDielectric-1/solventDielectric)*q1*q2/f;"
                              "f=sqrt(r^2+B1*B2*exp(-r^2/(4*B1*B2)))", mm.CustomGBForce.ParticlePairNoExclusions);

        # Add particles.
        for i in range(numMolecules):
            if (i < numMolecules/2):
                charge = 1.0 * unit.elementary_charge
                radius = 0.2 * unit.nanometers
                scale = 0.5
                nonbonded.addParticle(charge, sigma, epsilon)
                custom.addParticle([charge, radius, scale])

                charge = -1.0 * unit.elementary_charge
                radius = 0.1 * unit.nanometers
                scale = 0.5
                nonbonded.addParticle(charge, sigma, epsilon)            
                custom.addParticle([charge, radius, scale]);
            else:
                charge = 1.0 * unit.elementary_charge
                radius = 0.2 * unit.nanometers
                scale = 0.8
                nonbonded.addParticle(charge, sigma, epsilon)
                custom.addParticle([charge, radius, scale])

                charge = -1.0 * unit.elementary_charge
                radius = 0.1 * unit.nanometers
                scale = 0.8
                nonbonded.addParticle(charge, sigma, epsilon)            
                custom.addParticle([charge, radius, scale]);

        system.addForce(nonbonded)
        system.addForce(custom)    

        # Place particles at random positions within the box.
        # TODO: Use reproducible random number seed.
        # NOTE: This may not be thread-safe.
        
        state = np.random.get_state()
        np.random.seed(0)
        positions = unit.Quantity((boxSize/unit.nanometer) * np.random.rand(numParticles,3), unit.nanometer)
        np.random.set_state(state)

        self.system, self.positions = system, positions

#=============================================================================================
# AMOEBA SYSTEMS
#=============================================================================================

class AMOEBAIonBox(TestSystem):
    """A single Ca2 ion in a water box.

    >>> testsystem = AMOEBAIonBox()
    >>> system, positions = testsystem.system, testsystem.positions
    
    """
    def __init__(self):
        pdb_filename = get_data_filename("data/amoeba/ion-in-water.pdb")
        pdbfile = app.PDBFile(pdb_filename)

        ff =  app.ForceField("amoeba2009.xml")
        # TODO: 7A is a hack
        system = ff.createSystem(pdbfile.topology, nonbondedMethod=app.PME, constraints=app.HBonds, useDispersionCorrection=True, nonbondedCutoff=7.0*unit.angstroms)

        positions = pdbfile.getPositions()
        
        self.system, self.positions = system, positions

class AMOEBAProteinBox(TestSystem):
    """PDB 1AP4 in water box.

    >>> testsystem = AMOEBAProteinBox()
    >>> system, positions = testsystem.system, testsystem.positions    

    """
    def __init__(self):
        pdb_filename = get_data_filename("data/amoeba/1AP4_14_wat.pdb")
        pdbfile = app.PDBFile(pdb_filename)

        ff =  app.ForceField("amoeba2009.xml")
        system = ff.createSystem(pdbfile.topology, nonbondedMethod=app.PME, constraints=app.HBonds, useDispersionCorrection=True)

        positions = pdbfile.getPositions()
        
        self.system, self.positions = system, positions

#=============================================================================================
# ALCHEMICALLY MODIFIED SYSTEMS
#=============================================================================================

class AlchemicalState(object):
    """
    Alchemical state description.
        
    These parameters describe the parameters that affect computation of the energy.

    Attributes
    ----------
    relativeRestraints : float
        Scaling factor for remaining receptor-ligand relative restraint terms (to help keep ligand near protein).     
    ligandElectrostatics : float
        Scaling factor for ligand charges, intrinsic Born radii, and surface area term.
    ligandSterics : float
        Scaling factor for ligand sterics (Lennard-Jones and Halgren) interactions.
    ligandTorsions : float
        Scaling factor for ligand non-ring torsions.
    annihilateElectrostatics : bool
        If True, electrostatics should be annihilated, rather than decoupled.
    annihilateSterics : bool
        If True, sterics (Lennard-Jones or Halgren potential) will be annihilated, rather than decoupled.

    TODO
    ----
    * Rework these structure members into something more general and flexible?
    * Add receptor modulation back in?
    """
        
    def __init__(self, relativeRestraints=0.0, ligandElectrostatics=1.0, ligandSterics=1.0, ligandTorsions=1.0, annihilateElectrostatics=True, annihilateSterics=False):
        """
        Create an Alchemical state.

        Parameters
        ----------
        relativeRestraints : float, optional, default = 0.0
            Scaling factor for remaining receptor-ligand relative restraint terms (to help keep ligand near protein).     
        ligandElectrostatics : float, optional, default = 1.0
            Scaling factor for ligand charges, intrinsic Born radii, and surface area term.
        ligandSterics : float, optional, default = 1.0
            Scaling factor for ligand sterics (Lennard-Jones or Halgren) interactions.
        ligandTorsions : float, optional, default = 1.0
            Scaling factor for ligand non-ring torsions.
        annihilateElectrostatics : bool, optional, default = True
            If True, electrostatics should be annihilated, rather than decoupled.
        annihilateSterics : bool, optional, default = False
            If True, sterics (Lennard-Jones or Halgren potential) will be annihilated, rather than decoupled.

        Examples
        --------

        Create a fully-interacting, unrestrained alchemical state.
        
        >>> alchemical_state = AlchemicalState(relativeRestraints=0.0, ligandElectrostatics=1.0, ligandSterics=1.0, ligandTorsions=1.0)
        >>> # This is equivalent to
        >>> alchemical_state = AlchemicalState()


        Annihilate electrostatics.

        >>> alchemical_state = AlchemicalState(annihilateElectrostatics=True, ligandElectrostatics=0.0)

        """

        self.relativeRestraints = relativeRestraints
        self.ligandElectrostatics = ligandElectrostatics
        self.ligandSterics = ligandSterics
        self.ligandTorsions = ligandTorsions
        self.annihilateElectrostatics = annihilateElectrostatics
        self.annihilateSterics = annihilateSterics

        return

class AlchemicalTestSystem(object):
    def _alchemicallyModifyLennardJones(cls, system, nonbonded_force, alchemical_atom_indices, alchemical_state, alpha=0.50, a=1, b=1, c=6):
        """
        Alchemically modify the Lennard-Jones force terms.

        This version uses the new group-based restriction capabilities of CustomNonbondedForce.


        Parameters
        ----------
        system : simtk.openmm.System
        System to modify.
        nonbonded_force : simtk.openmm.NonbondedForce
        The NonbondedForce to modify (will be changed).
        alchemical_atom_indices : list of int
        Atom indices to be alchemically modified.
        alchemical_state : AlchemicalState
        The alchemical state specification to be used in modifying Lennard-Jones terms.
        alpha : float, optional, default = 0.5
        Alchemical softcore parameter.
        a, b, c : float, optional, default a=1, b=1, c=6
        Parameters describing softcore force.
        
        """

        import simtk.openmm as openmm

        # Create CustomNonbondedForce to handle softcore interactions between alchemically-modified system and rest of system.

        energy_expression = "4*epsilon*(lambda^a)*x*(x-1.0);"
        energy_expression += "x = (1.0/(alpha*(1.0-lambda)^b + (r/sigma)^c))^(6/c);" 
        energy_expression += "epsilon = sqrt(epsilon1*epsilon2);" # mixing rule for epsilon
        energy_expression += "sigma = 0.5*(sigma1 + sigma2);" # mixing rule for sigma
        energy_expression += "lambda = lennard_jones_lambda;" # lambda

        # Create atom groups.
        natoms = system.getNumParticles()
        atomset1 = set(alchemical_atom_indices) # only alchemically-modified atoms
        atomset2 = set(range(system.getNumParticles())) - atomset1 # all atoms minus intra-alchemical region

        # Create alchemically modified nonbonded force.
        # TODO: Create a _createCustomNonbondedForce method to duplicate parameters?
        energy_expression += "alpha = %f;" % alpha
        energy_expression += "a = %f; b = %f; c = %f;" % (a,b,c)    
        custom_nonbonded_force = openmm.CustomNonbondedForce(energy_expression)            
        custom_nonbonded_force.setNonbondedMethod(nonbonded_force.getNonbondedMethod()) # TODO: Make sure these method indices are identical.
        custom_nonbonded_force.setUseSwitchingFunction(nonbonded_force.getUseSwitchingFunction()) 
        custom_nonbonded_force.setSwitchingDistance(nonbonded_force.getSwitchingDistance()) 
        custom_nonbonded_force.setUseLongRangeCorrection(nonbonded_force.getUseDispersionCorrection())
        custom_nonbonded_force.addGlobalParameter("lennard_jones_lambda", alchemical_state.ligandSterics);
        custom_nonbonded_force.addPerParticleParameter("sigma") # Lennard-Jones sigma
        custom_nonbonded_force.addPerParticleParameter("epsilon") # Lennard-Jones epsilon

        # Restrict interaction evaluation to be between alchemical atoms and rest of environment.
        # Only add custom nonbonded force if interacting groups are both nonzero in size.
        if (len(atomset1) != 0) and (len(atomset2) != 0):
            custom_nonbonded_force.addInteractionGroup(atomset1, atomset2)
            system.addForce(custom_nonbonded_force)

        # Create CustomBondedForce to handle softcore exceptions if alchemically annihilating ligand.
        if alchemical_state.annihilateSterics:
            energy_expression = "4*epsilon*(lambda^a)*x*(x-1.0);"
            energy_expression += "x = (1.0/(alpha*(1.0-lambda)^b + (r/sigma)^c))^(6/c);" 
            energy_expression += "alpha = %f;" % alpha
            energy_expression += "a = %f; b = %f; c = %f;" % (a,b,c)
            energy_expression += "lambda = lennard_jones_lambda;"
            custom_bond_force = openmm.CustomBondForce(energy_expression)            
            custom_bond_force.addGlobalParameter("lennard_jones_lambda", alchemical_state.ligandSterics);
            custom_bond_force.addPerBondParameter("sigma") # Lennard-Jones sigma
            custom_bond_force.addPerBondParameter("epsilon") # Lennard-Jones epsilon
            system.addForce(custom_bond_force)
        else:
            # Decoupling of sterics.
            # Add a second CustomNonbondedForce to restore "intra-alchemical" interactions to full strength.
            energy_expression = "4*epsilon*((sigma/r)^12 - (sigma/r)^6);" 
            energy_expression += "epsilon = sqrt(epsilon1*epsilon2);" # mixing rule for epsilon
            energy_expression += "sigma = 0.5*(sigma1 + sigma2);" # mixing rule for sigma
            custom_nonbonded_force2 = openmm.CustomNonbondedForce(energy_expression)            
            custom_nonbonded_force2.setUseSwitchingFunction(nonbonded_force.getUseSwitchingFunction()) 
            custom_nonbonded_force2.setSwitchingDistance(nonbonded_force.getSwitchingDistance()) 
            custom_nonbonded_force2.setUseLongRangeCorrection(nonbonded_force.getUseDispersionCorrection())
            custom_nonbonded_force2.addPerParticleParameter("sigma") # Lennard-Jones sigma
            custom_nonbonded_force2.addPerParticleParameter("epsilon") # Lennard-Jones epsilon
            system.addForce(custom_nonbonded_force2)
            # Restrict interaction evaluation to be between alchemical atoms and rest of environment.
            atomset1 = set(alchemical_atom_indices) # only alchemically-modified atoms
            atomset2 = set(alchemical_atom_indices) # only alchemically-modified atoms
            custom_nonbonded_force2.addInteractionGroup(atomset1, atomset2)

        # Copy Lennard-Jones particle parameters.
        for particle_index in range(nonbonded_force.getNumParticles()):
            # Retrieve parameters.
            [charge, sigma, epsilon] = nonbonded_force.getParticleParameters(particle_index)
            # Add corresponding particle to softcore interactions.
            if particle_index in alchemical_atom_indices:
                # Turn off Lennard-Jones contribution from alchemically-modified particles.
                nonbonded_force.setParticleParameters(particle_index, charge, sigma, epsilon*0.0) 
            # Add contribution back to custom force.
            custom_nonbonded_force.addParticle([sigma, epsilon])            
            if not alchemical_state.annihilateSterics:
                custom_nonbonded_force2.addParticle([sigma, epsilon])

        # Create an exclusion for each exception in the reference NonbondedForce, assuming that NonbondedForce will handle them.
        for exception_index in range(nonbonded_force.getNumExceptions()):
            # Retrieve parameters.
            [iatom, jatom, chargeprod, sigma, epsilon] = nonbonded_force.getExceptionParameters(exception_index)
            # Exclude this atom pair in CustomNonbondedForce.
            custom_nonbonded_force.addExclusion(iatom, jatom)
            if not alchemical_state.annihilateSterics:
                custom_nonbonded_force2.addExclusion(iatom, jatom)

            # If annihilating Lennard-Jones, move intramolecular interactions to custom_bond_force.
            if alchemical_state.annihilateSterics and (iatom in alchemical_atom_indices) and (jatom in alchemical_atom_indices):
                # Remove Lennard-Jones exception.
                nonbonded_force.setExceptionParameters(exception_index, iatom, jatom, chargeprod, sigma, epsilon * 0.0)
                # Add special CustomBondForce term to handle alchemically-modified Lennard-Jones exception.
                custom_bond_force.addBond(iatom, jatom, [sigma, epsilon])

        # Set periodicity and cutoff parameters corresponding to reference Force.
        if nonbonded_force.getNonbondedMethod() in [openmm.NonbondedForce.Ewald, openmm.NonbondedForce.PME]:
            # Convert Ewald and PME to CutoffPeriodic.
            custom_nonbonded_force.setNonbondedMethod( openmm.CustomNonbondedForce.CutoffPeriodic )
            if not alchemical_state.annihilateSterics:
                custom_nonbonded_force2.setNonbondedMethod( openmm.CustomNonbondedForce.CutoffPeriodic )
        else:
            custom_nonbonded_force.setNonbondedMethod( nonbonded_force.getNonbondedMethod() )
            if not alchemical_state.annihilateSterics:
                custom_nonbonded_force2.setNonbondedMethod(nonbonded_force.getNonbondedMethod() )

        custom_nonbonded_force.setCutoffDistance( nonbonded_force.getCutoffDistance() )
        if not alchemical_state.annihilateSterics:
            custom_nonbonded_force2.setCutoffDistance( nonbonded_force.getCutoffDistance() )

        return

class AlchemicalLennardJonesCluster(TestSystem,AlchemicalTestSystem):
    """Create an alchemically-perturbed version of LennardJonesCluster.


    Parameters
    ----------
    nx : int, optional, default=3
        number of particles in the x direction
    ny : int, optional, default=3
        number of particles in the y direction
    nz : int, optional, default=3
        number of particles in the z direction        
    K : simtk.unit.Quantity, optional, default=1.0 * unit.kilojoules_per_mole/unit.nanometer**2
        harmonic restraining potential

    Examples
    --------

    Create Lennard-Jones cluster.
    
    >>> cluster = AlchemicalLennardJonesCluster()
    >>> system, positions = cluster.system, cluster.positions

    Create default 3x3x3 Lennard-Jones cluster in a harmonic restraining potential.

    >>> cluster = AlchemicalLennardJonesCluster(nx=10, ny=10, nz=10)
    >>> system, positions = cluster.system, cluster.positions
    """
    def __init__(self, nx=3, ny=3, nz=3, K=1.0 * unit.kilojoules_per_mole/unit.nanometer**2):

        # Default parameters
        mass_Ar     = 39.9 * unit.amu
        q_Ar        = 0.0 * unit.elementary_charge
        sigma_Ar    = 3.350 * unit.angstrom
        epsilon_Ar  = 0.001603 * unit.kilojoule_per_mole

        scaleStepSizeX = 1.0
        scaleStepSizeY = 1.0
        scaleStepSizeZ = 1.0

        # Determine total number of atoms.
        natoms = nx * ny * nz

        # Create an empty system object.
        system = mm.System()

        # Create a NonbondedForce object with no cutoff.
        nb = mm.NonbondedForce()
        nb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

        positions = unit.Quantity(np.zeros([natoms,3],np.float32), unit.angstrom)

        atom_index = 0
        for ii in range(nx):
            for jj in range(ny):
                for kk in range(nz):
                    system.addParticle(mass_Ar)
                    nb.addParticle(q_Ar, sigma_Ar, epsilon_Ar)
                    x = sigma_Ar*scaleStepSizeX*(ii - nx/2.0)
                    y = sigma_Ar*scaleStepSizeY*(jj - ny/2.0)
                    z = sigma_Ar*scaleStepSizeZ*(kk - nz/2.0)

                    positions[atom_index,0] = x
                    positions[atom_index,1] = y
                    positions[atom_index,2] = z
                    atom_index += 1

        # Add the nonbonded force.
        system.addForce(nb)

        # Add a restrining potential centered at the origin.
        force = mm.CustomExternalForce('(K/2.0) * (x^2 + y^2 + z^2)')
        force.addGlobalParameter('K', K)
        for particle_index in range(natoms):
            force.addParticle(particle_index, [])
        system.addForce(force)

        # Alchemically modify system.
        alchemical_atom_indices = [0]
        delta = 1.0e-5
        alchemical_state = AlchemicalState(0, 0, 1-delta, 1, annihilateElectrostatics=True, annihilateSterics=False)
        self._alchemicallyModifyLennardJones(system, nb, alchemical_atom_indices, alchemical_state)

        self.system, self.positions = system, positions



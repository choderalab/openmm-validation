# TestEnergyForces.py by Lee Periolat, 8/6/2012
#
# This script parses a set of Mark Friedrichs' input files and runs an energy-force
# consistency check on them.

from SystemTestFunctions import *
from simtk.openmm.app import *
import simtk.openmm as mm
import simtk.unit as unit
import numpy
import math
import sys
import os
from multiprocessing import Process, Queue

# usage message
def printUsage():
    print "TestEnergyForces.py requires at least one argument of the form systemFile=filename"
    print "filename must contain the system(s) to be tested."  
    print "It also accepts up to three additional arguments of the form argumentName=value:"
    print "platform:  A platform (Reference, CPU, OpenCL, or CUDA); default is platform=OpenCL"
    print "eps:       A value for epsilon in nanometers; default is eps=0.1"
    print "tol:       A tolerance for the relative force difference; default it tol=1e-3"
    print "precision: the precision to use (single, mixed, or double); default is single"
    print
    print "Example:   python TestEnergyForces.py systemFile=test1 platform=cuda tol=1e-2"  

def runTest(path, systemName, platformName, eps, tolerance, precision, queue):
    """This function runs the test for a single system."""
    print "\nsystemName =", systemName

    # load the system from the xml file
    system = loadXMLFile(path, systemName)
    
    # set Ewald error tolerance to a sufficiently small value so the requested accuracy will be achievable
    for f in system.getForces():
        try:
            f.setEwaldErrorTolerance(tolerance/2)
        except:
            pass

    # read in the particle positions from the .pos file
    positions = loadPosFile(path, systemName)

    numParticles = len(positions)
    print "numParticles =", numParticles, "(from the .pos file)"

    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)

    print "mm.Platform.getNumPlatforms =", mm.Platform.getNumPlatforms() 
    platform = mm.Platform.getPlatformByName(platformName)
    print "platform.getName() =", platform.getName()

    print "Building \'context\' on \'platform\'."
    context = mm.Context(system, integrator, platform, properties)
    context.setPositions(positions)
    state0 = context.getState(getForces=True)
    forces0 = state0.getForces()


    # make sure the size of the force vector is equal to the number of particles
    assert (len(forces0)==numParticles)

    # calculate the norm of the forces
    force0NormSum = 0.0*unit.kilojoules**2/unit.mole**2/unit.nanometer**2
    for f in forces0:
        force0NormSum += unit.dot(f,f)
    force0Norm = unit.sqrt(force0NormSum)
    print "force0Norm =", force0Norm

    epsilon = eps*unit.nanometer
    step = epsilon/force0Norm
    print "step =", step

    # perturb the coordinates along the direction of forces0 and evaluate the energy
    context.setPositions([p-2*f*step for p,f in zip(positions, forces0)])
    pe1 = context.getState(getEnergy=True).getPotentialEnergy()
    context.setPositions([p-f*step for p,f in zip(positions, forces0)])
    pe2 = context.getState(getEnergy=True).getPotentialEnergy()
    context.setPositions([p+f*step for p,f in zip(positions, forces0)])
    pe3 = context.getState(getEnergy=True).getPotentialEnergy()
    context.setPositions([p+2*f*step for p,f in zip(positions, forces0)])
    pe4 = context.getState(getEnergy=True).getPotentialEnergy()

    # use a finite difference approximation to calculate the expected force0Norm
    expectedForceNorm = (-pe1+8*pe2-8*pe3+pe4)/(12*epsilon)
    relativeDifference = abs(expectedForceNorm-force0Norm)/force0Norm

    print "pe1 =", pe1
    print "pe2 =", pe2
    print "pe3 =", pe3
    print "pe4 =", pe4
    print "expectedForceNorm =", expectedForceNorm
    print "relativeDifference =", relativeDifference

    # check that the energy is within the desired range
    if relativeDifference > tolerance:
        print "*** ERROR EXCEEDS TOLERANCE ***"
        queue.put(False)
    else:
        queue.put(True)
    print "Test of", systemName, "complete."
    del(context)

if __name__ == '__main__':
    if len(sys.argv)<2 or len(sys.argv)>6:
        printUsage()
        sys.exit(1)

    # set default values
    platformName = 'OpenCL'
    eps = 1e-1
    tolerance = 1e-3
    precision = 'single'

    # parse the argument list
    argList = sys.argv[1:]
    argNames = []
    for arg in argList:
        argSplit = arg.split("=")
        argName = argSplit[0].lower()
        argNames.append(argName)
        argValue = argSplit[1]
        if argName=='systemfile':
            systemFile = argValue
        elif argName=='platform':
            if argValue.lower()=='opencl': platformName = 'OpenCL'
            elif argValue.lower()=='cuda': platformName = 'CUDA'
            elif argValue.lower()=='reference': platformName = 'Reference'
            elif argValue.lower()=='cpu': platformName = 'CPU'
            else:
                print "Error: Argument 'platform' must be OpenCL, CUDA, CPU, or Reference (case insensitive)"
                printUsage()
                sys.exit(1)
        elif argName=='eps':
            eps = float(argValue)
        elif argName=='tol':
            tolerance = float(argValue) 
        elif argName=='precision':
            precision = argValue

    # ensure that the argument systemFile is provided
    if not 'systemfile' in argNames:
        print "Error:  Argument list must contain a systemFile"
        printUsage()
        sys.exit(1) 

    # set the path to the xml and pos files
    if os.name=='nt': path = ".\systems\\"
    else: path = './systems/'

    print "\nsystemFile =", systemFile
    print "platformName =", platformName
    print "eps =", eps
    print "tolerance =", tolerance
    print "precision =", precision

    if platformName == 'OpenCL':
        properties = {'OpenCLPrecision':precision}
    elif platformName == 'CUDA':
        properties = {'CudaPrecision':precision}
    else:
        properties = dict()

    # Load in the system names from the systems file
    systemNames = loadSystemFile(systemFile)

    # Run each test in its own process to avoid CUDA bugs when creating too many contexts
    # in the same process.
    failures = []
    for systemName in systemNames:
        queue = Queue()
        p = Process(target=runTest, args=(path, systemName, platformName, eps, tolerance, precision, queue))
        p.start()
        p.join()
        try:
            if (queue.get_nowait() == False):
                failures.append(systemName)
        except:
            failures.append(systemName)

    # List any tests that failed

    print
    if len(failures) == 0:
        print "All tests passed"
    else:
        print len(failures), "tests failed:"
        for test in failures:
            print test

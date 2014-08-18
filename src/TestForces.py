# TestForces.py by Lee Periolat, 8/6/2012
#
# This script parses a set of Mark Friedrichs' input files, builds contexts on both
# the reference platform and a user chosen platform, and checks that the forces are
# consistent between the two.

from SystemTestFunctions import *
from simtk.openmm.app import *
import simtk.openmm as mm
import simtk.unit as unit
import numpy
import sys
import os
from multiprocessing import Process, Queue

def printUsage():
    print "TestForces.py requires at least one argument of the form systemFile=filename"
    print "filename must contain the system(s) to be tested."  
    print "It also accepts up to two additional arguments of the form argumentName=value:"
    print "platform:  A platform (CPU, OpenCL, or CUDA); default is platform=OpenCL"
    print "tol:       A tolerance for the ninetieth percentile of the relative force difference;"
    print "           default is tol=1e-4"
    print "precision: the precision to use (single, mixed, or double); default is single"
    print
    print "Example:   python TestForces.py systemFile=test1 platform=cuda tol=1e-2"  

def runTest(path, systemName, platformName, tolerance, precision, queue):
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

    # create the two integrators for comparison
    integrator0 = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    integrator1 = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)

    print "mm.Platform.getNumPlatforms =", mm.Platform.getNumPlatforms()
    platform0 = mm.Platform.getPlatformByName("Reference")
    platform1 = mm.Platform.getPlatformByName(platformName)
    print "platform0.getName() =", platform0.getName()
    print "platform1.getName() =", platform1.getName()

    print "Building context0 on platform0, the reference platform."
    context0 = mm.Context(system, integrator0, platform0)
    context0.setPositions(positions)
    state0 = context0.getState(getForces=True)
    forces0 = state0.getForces().value_in_unit(unit.kilojoules/unit.mole/unit.nanometer)

    print "Building context1 on platform1, the", platformName, "platform."
    context1 = mm.Context(system, integrator1, platform1, properties)
    context1.setPositions(positions)
    state1 = context1.getState(getForces=True)
    forces1 = state1.getForces(asNumpy=True).value_in_unit(unit.kilojoules/unit.mole/unit.nanometer)

    # make sure the two force vectors are the same length, and equal to the number of particles
    assert (len(forces0)==len(forces1))
    assert (len(forces0)==numParticles)

    # Calculate the relative difference for each particle in the system.  Along the way, look for
    # the maximum relative difference so that we can save the force norm at the maximum.
    sumRelativeDifferences = 0
    maxRelativeDifference1 = 0
    forceNormOfMax = 0
    relativeDifferences = [0]*numParticles
    forcesDiff = forces1 - forces0
    for i, f0 in enumerate(forces0):
        f1 = forces1[i]
        fdiff = forcesDiff[i]
        normf0 = numpy.linalg.norm(f0)
        normf1 = numpy.linalg.norm(f1)
        normfdiff = numpy.linalg.norm(fdiff)
        #if normf1 == 0: print "normf1 is zero"
        if normfdiff > 0:
            relativeDifferences[i] = 2.0*normfdiff/(normf1+normf0) 
        sumRelativeDifferences += relativeDifferences[i]
        if relativeDifferences[i] > maxRelativeDifference1:
            maxRelativeDifference1 = relativeDifferences[i]
            forceNormOfMax = normf1
    averageRelativeDifference = sumRelativeDifferences/float(numParticles)

    # calculate min, max, median, 10th percentile, and 90th percentile relative differences
    relativeDifferences.sort()
    maxRelativeDifference2 = relativeDifferences[numParticles-1]
    minRelativeDifference = relativeDifferences[0]
    medianRelativeDifference = relativeDifferences[int(round(0.5*numParticles+0.5))-1]
    tenthPercentileRelativeDifference = relativeDifferences[int(round(0.1*numParticles+0.5))-1]
    ninetiethPercentileRelativeDifference = relativeDifferences[int(round(0.9*numParticles+0.5))-1]

    print "averageRelativeDifference", averageRelativeDifference
    print "maxRelativeDifference1", maxRelativeDifference1
    print "maxRelativeDifference2", maxRelativeDifference2
    print "force norm of max relative difference", forceNormOfMax
    print "minRelativeDifference", minRelativeDifference
    print "medianRelativeDifference", medianRelativeDifference
    print "tenthPercentileRelativeDifference", tenthPercentileRelativeDifference
    print "ninetiethPercentileRelativeDifference", ninetiethPercentileRelativeDifference
    if ninetiethPercentileRelativeDifference > tolerance:
        print "*** ERROR EXCEEDS TOLERANCE ***"
        queue.put(False)
    else:
        queue.put(True)
    print "Test of", systemName, "complete."
    del(context0)
    del(context1)

if __name__ == '__main__':
    if len(sys.argv)<2 or len(sys.argv)>5:
        printUsage()
        sys.exit(1)

    # set default values
    platformName = 'OpenCL'
    tolerance = 1e-4
    precision = 'single'

    # parse the argument list
    argList = sys.argv[1:]
    argNames = []
    for arg in argList:
        argSplit = arg.split("=")
        argName = argSplit[0].lower();
        argNames.append(argName)
        argValue = argSplit[1]
        if argName=='systemfile':
            systemFile = argValue
        elif argName=='platform':
            if argValue.lower()=='opencl': platformName = 'OpenCL'
            elif argValue.lower()=='cuda': platformName = 'CUDA'
            elif argValue.lower()=='cpu': platformName = 'CPU'
            else:
                print "Error:  Argument 'platform' must be CPU, OpenCL, or CUDA (case insensitive)"
                printUsage()
                sys.exit(1)
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
        p = Process(target=runTest, args=(path, systemName, platformName, tolerance, precision, queue))
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

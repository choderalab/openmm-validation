# TestVerletEnergyConservation.py by Lee Periolat, 8/18/2012
#
# This script parses two of Mark Friedrichs' input files and runs an energy conservation test on them
# using either a standard Verlet integrator or a variable Verlet integrator.

from SystemTestFunctions import *
from simtk.openmm.app import *
import simtk.openmm as mm
import simtk.unit as unit
import numpy as np
import sys
import os
from multiprocessing import Process, Queue

def printUsage():
    print "TestVerletEnergyConservation.py requires at least one argument of the form systemFile=filename"
    print "filname must contain the system(s) to be tested."
    print "It also accepts up to six additional arguments of the form argumentName=value:"
    print "platform:     A platform (OpenCL, CUDA, CPU, or reference); default is platform=OpenCL"
    print "constTol:     A tolerance for the system constraints; default is constTol=1e-4"
    print "variableInt:  Whether to use a VariableVerletIntegrator or VerletIntegrator (True, t, False, or f);" 
    print "              default is variableInt=False, which means a VerletIntegrator is used by default"
    print "checks:       Number of times to check the energy and constraints; default is checks=1000"
    print "timeBtwChks:  Time between checks in picoseconds; default is timeBtwChks=1"
    print "errorTol:     The error tolerance for VariableLangevinIntegrator; default is errorTol=2e-6"
    print "stepSize:     The step size in picoseconds for VerletIntegrator; default is stepSize=1e-3"
    print "precision: the precision to use (single, mixed, or double); default is single"
    print
    print "Example:      python TestVerletEnergyConservation.py systemFile=test1 platform=opencl energyTol=1e-2"  

def runTest(path, systemName, platformName, constraintTolerance, variableIntegrator, checks, timeBetweenChecks,
            errorTolerance, stepSize, precision, queue):
    """This function runs the test for a single system."""
    print "\nsystemName =", systemName

    # load the system from the xml file
    system = loadXMLFile(path, systemName)

    # read in the particle positions from the .pos file
    positions = loadPosFile(path, systemName)

    numParticles = len(positions)
    print "numParticles =", numParticles, "(from the .pos file)"

    # calculate degrees of freedom in the system
    particles = system.getNumParticles()
    constraints = system.getNumConstraints()
    forces = [system.getForce(i) for i in range(system.getNumForces())]
    usesCMMR = any(isinstance(f, mm.CMMotionRemover) for f in forces)
    DOF = 3*particles-constraints 
    if usesCMMR: DOF -= 3

    print "particles =", particles, "(from the xml file)"
    print "constraints =", constraints
    print "usesCMMR =", usesCMMR
    print "DOF =", DOF

    # set simulation parameters
    BOLTZ = unit.MOLAR_GAS_CONSTANT_R
    equilibrationSteps = 30000
    temperatureEq = 300*unit.kelvin
    frictionCoeffEq = 1/unit.picosecond
    stepSizeEq = 0.001*unit.picoseconds
    energyTotal = 0.0*unit.kilojoules/unit.mole
    
    print "equilibrationSteps =", equilibrationSteps
    print "temperatureEq =", temperatureEq
    print "frictionCoeffEq =", frictionCoeffEq
    print "stepSizeEq =", stepSizeEq
    print "variableIntegrator =", variableIntegrator
    print "checks =", checks
    print "timeBetweenChecks =", timeBetweenChecks
    if variableIntegrator: print "errorTolerance =", errorTolerance
    else: 
        stepsBetweenChecks = int(timeBetweenChecks/stepSize)  
        print "stepSize =", stepSize
        print "stepsBetweenChecks =", stepsBetweenChecks
    print "constraintTolerance =", constraintTolerance
    #print "energyTolerance =", energyTolerance
    
    # create the equilibration integrator and context
    equilibrationIntegrator = mm.LangevinIntegrator(temperatureEq, frictionCoeffEq, stepSizeEq)
    equilibrationIntegrator.setConstraintTolerance(constraintTolerance)
    platform = mm.Platform.getPlatformByName(platformName)
    equilibrationContext = mm.Context(system, equilibrationIntegrator, platform)

    # set the positions according to the .pos file
    equilibrationContext.setPositions(positions)

    # equilibrate the system
    equilibrationContext.setVelocitiesToTemperature(temperatureEq)
    equilibrationContext.getIntegrator().step(equilibrationSteps)

    # get the positions and velocities in the equilibrium state
    equilibriumState = equilibrationContext.getState(getPositions=True, getVelocities=True)
    equilibriumPositions = equilibriumState.getPositions()
    equilibriumVelocities = equilibriumState.getVelocities()
    del(equilibrationContext)

    # create the simulation integrator, either variable or standard
    if variableIntegrator:
        integrator = mm.VariableVerletIntegrator(errorTolerance)
    else:
        integrator = mm.VerletIntegrator(stepSize)
    integrator.setConstraintTolerance(constraintTolerance)

    # create the context
    context = mm.Context(system, integrator, platform, properties)

    # transfer the positions and velocites from equilibration context to simulation context
    context.setPositions(equilibriumPositions)
    context.setVelocities(equilibriumVelocities)

    # run the simulation, accumulating the kinetic energy and checking constraints
    energiesOverDOF = [0]*(checks-1)
    times = [0]*(checks-1)
    passed = True
    for i in range(checks):
        state = context.getState(getEnergy=True, getPositions=True)
        if i>0: 
            energyOverDOF = computeEnergy(state, system, integrator.getStepSize())/DOF
            #print "i =", i, " energyOverDOF =", energyOverDOF
            energiesOverDOF[i-1] = energyOverDOF.value_in_unit(unit.kilojoules/unit.mole)
            times[i-1] = i*timeBetweenChecks.value_in_unit(unit.nanoseconds)
            #if i%10==0: print "energyOverDOF =", energyOverDOF
        # check the constraints to make sure they are still satisfied
        try:
            checkConstraints(state, system, constraintTolerance)
        except:
            print '*** CONSTRAINT VIOLATED ON CHECK', i+1, '***'
            passed = False
        # advance the simulation
        if variableIntegrator:
            integrator.stepTo((i+1)*timeBetweenChecks)
        else:
            integrator.step(stepsBetweenChecks)

    print "Simulation results:"
    if constraints>0: print "All constraints satisfied."
    
    # perform a linear regression to find the energy drift in units kT/DOF/ns
    x = np.array(times)
    y = np.array(energiesOverDOF)
    A = np.vstack([x, np.ones(len(x))]).T
    slope, intercept = np.linalg.lstsq(A, y)[0]
    slopeWithUnits = slope*unit.kilojoules/unit.moles/unit.nanoseconds
    print "slope =", slope
    print "slopeWithUnits = energyDrift = ", slopeWithUnits
    energyDrift = slopeWithUnits
    print "intercept =", intercept
    interceptWithUnits = intercept*unit.kilojoules/unit.mole
    print "interceptWithUnits =", interceptWithUnits
    print "Test of", systemName, "complete."
    queue.put(passed)
    del(context)

if __name__ == '__main__':
    if len(sys.argv)<2 or len(sys.argv)>9:
        printUsage()
        sys.exit(1)
    elif len(sys.argv)>=2:
        systemsFile = sys.argv[1]

    # set default values
    platformName = 'OpenCL'
    constraintTolerance = 1e-4
    variableIntegrator = False
    checks = 1000
    timeBetweenChecks = 1*unit.picoseconds
    errorTolerance = 2e-6
    stepSize = 1e-3*unit.picoseconds
    precision = 'single'

    # parse the argument list
    argList = sys.argv[1:]
    argNames = []
    for arg in argList:
        argSplit = arg.split("=")
        argName = argSplit[0]
        argNames.append(argName.lower())
        argValue = argSplit[1]
        if argName.lower()=='systemfile':
            systemFile = argValue
        elif argName.lower()=='platform':
            if argValue.lower()=='opencl': platformName = 'OpenCL'
            elif argValue.lower()=='cuda': platformName = 'CUDA'
            elif argValue.lower()=='cpu': platformName = 'CPU'
            else:
                print "Error:  Argument 'platform' must be CPU, OpenCL, or CUDA (case insensitive)"
                printUsage()
                sys.exit(1)
        elif argName.lower()=='consttol':
            constraintTolerance = float(argValue)
        elif argName.lower()=='variableint':
            if argValue.lower()=='t' or argValue.lower()=='true':
                variableIntegrator=True
            elif argValue.lower()=='f' or argValue.lower()=='false':
                variableIntegrator=False
            else:
                print "Error:  Argument 'variableInt' must be True, t, False, or f (case insensitive)"
                printUsage()
                sys.exit(1)
        elif argName.lower()=='checks':
            checks = int(argValue)
        elif argName.lower()=='timebtwchks':
            timeBetweenChecks = float(argValue)*unit.picoseconds
        elif argName.lower()=='errortol':
            errorTolerance = float(argValue)
        elif argName.lower()=='stepsize':
            stepSize = float(argValue)*unit.picoseconds
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

    # display input or default parameter values
    print "\nsystemFile =", systemFile
    print "platformName =", platformName
    print "constraintTolerance =", constraintTolerance
    print "variableIntegrator =", variableIntegrator
    print "checks =", checks
    print "timeBetweenChecks =", timeBetweenChecks
    if variableIntegrator:
        print "errorTolerance =", errorTolerance
    else:
        print "stepSize =", stepSize
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
        p = Process(target=runTest, args=(path, systemName, platformName, constraintTolerance, variableIntegrator,
                    checks, timeBetweenChecks, errorTolerance, stepSize, precision, queue))
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

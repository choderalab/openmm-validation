# TestLangevinThermostability.py by Lee Periolat, 8/6/2012
#
# This script parses a set of Mark Friedrichs' input files and runs a thermostability test on them
# using either a standard Langevin integrator or a variable Langevin integrator.

from SystemTestFunctions import *
from simtk.openmm.app import *
import simtk.openmm as mm
import simtk.unit as unit
import sys
import os
from multiprocessing import Process, Queue

def printUsage():
    print "TestLangevinThermostability.py requires at least one argument of the form systemFile=filename"
    print "filname must contain the system(s) to be tested."
    print "It also accepts up to nine additional arguments of the form argumentName=value:"
    print "platform:     A platform (CPU, OpenCL, CUDA, or reference); default is platform=OpenCL"
    print "constTol:     A tolerance for the system constraints; default is constTol=1e-4"
    print "energyTol:    A tolerance for the system energy; default is energyTol=0.03"
    print "variableInt:  Whether to use a VariableLangevinIntegrator or LangevinIntegrator (True, t, False, or f);" 
    print "              default is variableInt=False, which means a LangevinIntegrator is used by default"
    print "checks:       Number of times to check the energy and constraints; default is checks=1000"
    print "timeBtwChks:  Time between checks in picoseconds; default is timeBtwChks=1"
    print "temp:         The temperature in Kelvin of the simulation; default it temp=300"
    print "friction:     The friction coefficient in 1/picoseconds of the simulation; default is friction=1"
    print "errorTol:     The error tolerance for VariableLangevinIntegrator; default is errorTol=2e-6"
    print "stepSize:     The step size in picoseconds for LangevinIntegrator; default is stepSize=1e-3"
    print "precision: the precision to use (single, mixed, or double); default is single"
    print
    print "Example:      python TestLangevinThermostability.py systemFile=test1 platform=cuda energyTol=1e-2"  

def runTest(path, systemName, platformName, constraintTolerance, energyTolerance, variableIntegrator, checks,
            timeBetweenChecks, temperature, frictionCoeff, errorTolerance, stepSize, precision, queue):
    """This function runs the test for a single system."""
    print "\nsystemName =", systemName

    # load the system from the xml file
    system = loadXMLFile(path, systemName)

    # read in the particle positions from the .pos file
    positions = loadPosFile(path, systemName)

    numParticles = len(positions)
    print "numParticles =", numParticles, "(from .pos file)"

    # calculate degrees of freedom in the system
    particles = system.getNumParticles()
    constraints = system.getNumConstraints()
    forces = [system.getForce(i) for i in range(system.getNumForces())]
    usesCMMR = any(isinstance(f, mm.CMMotionRemover) for f in forces)
    DOF = 3*particles-constraints 
    if usesCMMR: DOF -= 3

    print "particles =", particles, "(from xml file)"
    print "constraints =", constraints
    print "usesCMMR =", usesCMMR
    print "DOF =", DOF

    # set simulation parameters
    BOLTZ = unit.MOLAR_GAS_CONSTANT_R
    equilibrationSteps = 30000
    temperatureEq = 300*unit.kelvin
    frictionCoeffEq = 1/unit.picosecond
    stepSizeEq = 0.001*unit.picoseconds
    kineticEnergyTotal = 0.0*unit.kilojoules/unit.mole
    print "temperatureEq =", temperatureEq
    print "frictionCoeffEq =", frictionCoeffEq
    print "stepSizeEq =", stepSizeEq    
    print "equilibrationSteps =", equilibrationSteps
    print "variableIntegrator =", variableIntegrator
    print "checks =", checks
    print "timeBetweenChecks =", timeBetweenChecks
    if variableIntegrator: print "errorTolerance =", errorTolerance
    else:
        stepsBetweenChecks = int(timeBetweenChecks/stepSize)
        print "stepSize =", stepSize
        print "stepsBetweenChecks =", stepsBetweenChecks
    print "constraintTolerance =", constraintTolerance
    print "energyTolerance =", energyTolerance

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
        integrator = mm.VariableLangevinIntegrator(temperature, frictionCoeff, errorTolerance)
    else:
        integrator = mm.LangevinIntegrator(temperature, frictionCoeff, stepSize)
    integrator.setConstraintTolerance(constraintTolerance)
    
    # create the context
    context = mm.Context(system, integrator, platform, properties)

    # transfer the positions and velocities from the equilibration context to simulation context
    context.setPositions(equilibriumPositions)
    context.setVelocities(equilibriumVelocities)

    # run the simulation, accumulating the kinetic energy and checking constraints
    passed = True
    for i in range(checks):
        state = context.getState(getEnergy=True, getPositions=True)
        kineticEnergyTotal += state.getKineticEnergy()
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

    # calculate the average kinetic energy, temperature, and expected kinetic energy
    averageKineticEnergy = kineticEnergyTotal/checks
    averageTemperature = averageKineticEnergy/(0.5*DOF*BOLTZ)
    expectedKE = 0.5*DOF*BOLTZ*temperature
    relativeDifference = abs((averageKineticEnergy-expectedKE)/averageKineticEnergy)
    
    # print results
    print "Simulation results:"
    if constraints>0: print "All constraints satisfied."
    print "averageKineticEnergy =", averageKineticEnergy
    print "expectedKE =", expectedKE
    print "relativeDifference =", relativeDifference
    print "averageTemperature =", averageTemperature

    # assert that the kinetic energy is close to the expected value
    if relativeDifference > energyTolerance:
        print "*** ERROR EXCEEDS TOLERANCE ***"
        passed = False
    print "Test of", systemName, "complete."
    queue.put(passed)
    del(context)

if __name__ == '__main__':
    if len(sys.argv)<2 or len(sys.argv)>12:
        printUsage()
        sys.exit(1)
    elif len(sys.argv)>=2:
        systemsFile = sys.argv[1]

    # set default values
    platformName = 'OpenCL'
    constraintTolerance = 1e-4
    energyTolerance = 0.03
    variableIntegrator = False
    checks = 1000
    timeBetweenChecks = 1*unit.picoseconds
    temperature = 300*unit.kelvin
    frictionCoeff = 1/unit.picoseconds
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
        elif argName.lower()=='energytol':
            energyTolerance = float(argValue) 
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
        elif argName.lower()=='temp':
            temperature = float(argValue)*unit.kelvin
        elif argName.lower()=='friction':
            frictionCoeff = float(argValue)/unit.picoseconds
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
    print "energyTolerance =", energyTolerance
    print "variableIntegrator =", variableIntegrator
    print "checks =", checks
    print "timeBetweenChecks =", timeBetweenChecks
    print "temperature =", temperature
    print "frictionCoeff =", frictionCoeff
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
        p = Process(target=runTest, args=(path, systemName, platformName, constraintTolerance, energyTolerance,
                    variableIntegrator, checks, timeBetweenChecks, temperature, frictionCoeff, errorTolerance,
                    stepSize, precision, queue))
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

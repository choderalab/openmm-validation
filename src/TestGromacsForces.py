from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
import os
import subprocess
import re

def printUsage():
    print "TestGromacsForces.py requires at least one argument of the form filename=(filename)"
    print "filename should be the name (without extension) of the top/gro/mdp files to load"
    print "It also accepts additional arguments of the form argumentName=value:"
    print "platform:  A platform (CPU, OpenCL, or CUDA); default is platform=OpenCL"
    print "tol:       A tolerance for the ninetieth percentile of the relative force difference;"
    print "           default is tol=1e-4"
    print "precision: The precision to use (single, mixed, or double); default is single"
    print "implicit:  Whether to use implicit solvent (true or false); default is false"
    print
    print "Example:   python TestGromacsForces.py filename=implicit platform=cuda tol=1e-2"  

# Set default values

filename = None
platformName = 'OpenCL'
tolerance = 1e-4
precision = 'single'
implicit = False

# Parse the argument list

for arg in sys.argv[1:]:
    argSplit = arg.split("=")
    argName = argSplit[0].lower();
    argValue = argSplit[1]
    if argName=='filename':
        filename = argValue
    elif argName=='platform':
        if argValue.lower()=='opencl': platformName = 'OpenCL'
        elif argValue.lower()=='cuda': platformName = 'CUDA'
        elif argValue.lower()=='cpu': platformName = 'CPU'
        else:
            print "Error:  Argument 'platform' must be OpenCL or CUDA (case insensitive)"
            printUsage()
            sys.exit(1)
    elif argName=='tol':
        tolerance = float(argValue)
    elif argName=='precision':
        precision = argValue
    elif argName=='implicit':
        if argValue.lower()=='true': implicit = True
        elif argValue.lower()=='false': implicit = False
        else:
            print "Error: Argument 'implicit' must be true or false (case insensitive)"
            printUsage()
            sys.exit(1)

# Ensure that the argument filename is provided

if filename is None:
    print "Error: Argument list must contain a filename"
    printUsage()
    sys.exit(1) 

print "filename =", filename
print "platformName =", platformName
print "tolerance =", tolerance
print "precision =", precision

if platformName == 'OpenCL':
    properties = {'OpenCLPrecision':precision}
elif platformName == 'CUDA':
    properties = {'CudaPrecision':precision}
else:
    properties = dict()

# Compute forces with Gromacs

gromacsBinDir = '/usr/local/gromacs/bin/'
gromacsTopDir = '/usr/local/gromacs/share/gromacs/top'
subprocess.call([gromacsBinDir+'grompp', '-f', filename+'.mdp', '-p', filename+'.top', '-c', filename+'.gro'])
subprocess.call([gromacsBinDir+'mdrun', '-nt', '1'])
process = subprocess.Popen([gromacsBinDir+'gmxdump', '-f', 'traj.trr'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
(output, errors) = process.communicate()
expr = re.compile('f\[.*?\]=\{(.*?),(.*?),(.*?)\}')
gromacsForces = []
for match in re.findall(expr, output):
    gromacsForces.append(Vec3(*[float(x) for x in match]))

# Compute forces with OpenMM

gro = GromacsGroFile(filename+'.gro')
top = GromacsTopFile(filename+'.top', unitCellDimensions=gro.getUnitCellDimensions(), includeDir=gromacsTopDir)
mdpOptions = {}
for line in open(filename+'.mdp'):
    if '=' in line:
        fields = [x.strip().lower() for x in line.split('=')]
    mdpOptions[fields[0]] = fields[1]
if mdpOptions.get('implicit_solvent') == 'gbsa':
    system = top.createSystem(nonbondedMethod=NoCutoff, implicitSolvent=OBC2, rigidWater=False)
elif mdpOptions.get('coulombtype') == 'pme':
    system = top.createSystem(nonbondedMethod=PME, ewaldErrorTolerance=1e-6, rigidWater=False)
else:
    system = top.createSystem(nonbondedMethod=NoCutoff, rigidWater=False)

integrator = VerletIntegrator(0.001*picoseconds)
simulation = Simulation(top.topology, system, integrator, Platform.getPlatformByName(platformName), properties)
simulation.context.setPositions(gro.positions)
state = simulation.context.getState(getForces=True)
forces = state.getForces().value_in_unit(kilojoules_per_mole/nanometer)
relativeDiff = [2*norm(f1-f2)/(norm(f1)+norm(f2)) for f1, f2 in zip(forces, gromacsForces)]
relativeDiff.sort()
print '10%:', relativeDiff[int(0.1*len(relativeDiff))]
print '50%:', relativeDiff[len(relativeDiff)/2]
print '90%:', relativeDiff[int(0.9*len(relativeDiff))]
print 'max:', relativeDiff[len(relativeDiff)-1]


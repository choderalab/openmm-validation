""" These functions are used in the system test scripts. """

import simtk.openmm as mm
import simtk.unit as unit
import os.path
import sys
import math

def loadSystemFile(systemFile):
    f = open(systemFile, 'r')
    lines = f.readlines()
    # remove trainling newlines, blank lines, and lines that start with #
    return [line.rstrip() for line in lines if len(line.rstrip())>0 and line[0]!="#"]

def loadXMLFile(path, systemName):
    filename = path + systemName + '.xml'
    f = open(filename, 'r')
    fxml = f.read()
    return mm.XmlSerializer.deserializeSystem(fxml)

def loadPosFile(path, systemName):
    filename = path + systemName + '.pos'
    while (not os.path.exists(filename)) and '_' in systemName:
        systemName = systemName[:systemName.rfind('_')]
        filename = path+systemName+'.pos'
    fpos = open(filename)

    # The first line of the .pos file should contain the number of positions.  Sometimes
    # the word 'Positions' is first and the number is second, sometimes the number is first
    # and the word 'Positions' is second.
    firstLine = fpos.readline()
    firstLineList = firstLine.split()
    if firstLineList[1] != 'Positions' and firstLineList[0] != 'Positions':
        print "Error; the first line of the .pos file should contain the number of particles."
        sys.exit(1)
    if firstLineList[1] == 'Positions':
        numParticles = int(firstLineList[0])
    else:
        numParticles = int(firstLineList[1])

    # read the remaining lines from the .pos file and parse the positions into the correct
    # data structure
    positions = [None]*numParticles
    remainingLines = fpos.readlines()
    for line in remainingLines:
        lineSplit = line.split()
        index = int(lineSplit[0])
        pos1 = float(lineSplit[1])
        pos2 = float(lineSplit[2])
        pos3 = float(lineSplit[3])
        positions[index] = mm.Vec3(pos1, pos2, pos3)*unit.nanometer
    return positions
# end of loadPosFile()

def computeEnergy(state, system, dt):
    """ Compute the energy of a state.
    """
    return state.getKineticEnergy()+state.getPotentialEnergy()
# end of computeEnergy()

def checkConstraints(state, system, constraintTolerance):
    positions = state.getPositions().value_in_unit(unit.nanometer)
    for j in range(system.getNumConstraints()):
        (particle1, particle2, expectedDistance) = system.getConstraintParameters(j)
        expectedDistance = expectedDistance.value_in_unit(unit.nanometer)
        p1 = positions[particle1]
        p2 = positions[particle2]
        distance = math.sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2+(p1[2]-p2[2])**2)
        relativeDifference = abs((distance-expectedDistance)/expectedDistance)
        assert (relativeDifference<constraintTolerance)
# end of checkConstraints()

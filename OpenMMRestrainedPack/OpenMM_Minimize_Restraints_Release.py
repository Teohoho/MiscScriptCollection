from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import datetime

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--files", type=str, help="Roots of files you want to Min,Heat,Equi; needs to be formatted like files.prmtop and files.inpcrd")
parser.add_argument("--OutputRoot", type=str, help="A string to be appended to the begining of the DCD,OUT,CHK and XML files generated by openMM")
parser.add_argument("--LoadState", type=str, help="XML file from previous simulation. Optional")
parser.add_argument("--RestrainedAtomsIn", type=str, help="A file containing the indices of the atoms to be restrained")
parser.add_argument("--Verbose", action="store_true", help="Print PDBs after every restrained minimization. If not used, only the positions obtained the final minimization will be saved.")
args = parser.parse_args()

## In order to reinitialize the restraints to the new positions
## of the particles, we need to generate the force objects as we 
## go. I will write a function that takes a set of positions and 
## a set of indices and returns a force object.

def GenerateRestraint(Positions, Indices):

	"""A function that generates an OpenMM Force object given
	   a set of positions and a list of atom indices that said 
	   force must affect.

	Parameters
   ----------

	Positions: simtk.unit.quantity.Quantity
		list of positions used in setting up the harmonic restraints

	Indices: list of int
		list of atom indices that will be affected by the force generated

	Returns
	------

	GeneratedRestraint: simtk.openmm.openmm.CustomExternalForce
		harmonic potential force object

	"""

	CustomHarmonicForce = CustomExternalForce("springConstant * (periodicdistance(x, y, z, x0, y0, z0)^2)")     #Keep periodicity in mind?
	CustomHarmonicForce.addGlobalParameter("springConstant", (10)*kilocalories_per_mole/nanometer**2)
	CustomHarmonicForce.addPerParticleParameter("x0")
	CustomHarmonicForce.addPerParticleParameter("y0")
	CustomHarmonicForce.addPerParticleParameter("z0")

	for PartIx in Indices:
		PartX = Positions[PartIx][0].value_in_unit(nanometer)
		PartY = Positions[PartIx][1].value_in_unit(nanometer)
		PartZ = Positions[PartIx][2].value_in_unit(nanometer)

		CustomHarmonicForce.addParticle(PartIx, [PartX, PartY, PartZ])


	print ("A new restraint was generated, containing {} atoms.".format(len(Indices)))

	return CustomHarmonicForce

# Read files
prmtop = AmberPrmtopFile(args.files + ".prmtop")
inpcrd = AmberInpcrdFile(args.files + ".inpcrd")

# Thermodynamics
T = 300*kelvin
p = 1*bar
timestep = 0.002*picoseconds

# First we create the system, context, simulation objects
#system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
#system.addForce(MonteCarloBarostat(p, T))
system = prmtop.createSystem(implicitSolvent=OBC2, soluteDielectric=1.0, solventDielectric=80.0, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1.2*nanometer, constraints=HBonds, implicitSolventSaltConc=0.15*moles/liter)

integrator = LangevinIntegrator(T, 1/picosecond, timestep)

simulation = Simulation(prmtop.topology, system, integrator, platformProperties={'DisablePmeStream':'true', 'DeviceIndex':'0,1'})
simulation.context.setPositions(inpcrd.positions) #We set the positions of the particles to be the same as in the inpcrd
if inpcrd.boxVectors is not None:
	simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# We iterate through all the sets of restraints, creating a list
# of all the index groups
if (args.RestrainedAtomsIn is not None):
	FlexList  = open(args.RestrainedAtomsIn) 
	LoopAtoms = FlexList.read().split("\n") 
	LoopAtoms = [x for x in LoopAtoms if x]  ##Remove empty elements 

	for RestraintIx in range(len(LoopAtoms)):
		LoopAtoms[RestraintIx] = LoopAtoms[RestraintIx].split()   ##Split into n groups
		LoopAtoms[RestraintIx] = [int(x) for x in LoopAtoms[RestraintIx] if x]  ##Turn to ints


	print ("There are {} restraint profiles:".format(len(LoopAtoms)))
	for RestraintIx in range(len(LoopAtoms)):
		print ("\t· Restraint Set {}: {} atoms restrained out of a total of {} atoms.".format(RestraintIx, len(LoopAtoms[RestraintIx]),system.getNumParticles()))

	##Minimize
	##We have to minimize while iterating through all the constraints set.

	TotalPotEnergy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
	print ("Before any minimization is done, the energy is {} kcal/mole".format(TotalPotEnergy))

	for RestraintIx in range(len(LoopAtoms)):

		CurrentPositions = simulation.context.getState(getPositions=True).getPositions()
		AtomIndices = LoopAtoms[RestraintIx]

		CurrentRestraint = GenerateRestraint(CurrentPositions, AtomIndices)

		system.addForce(CurrentRestraint)
		simulation.context.reinitialize(preserveState=True)

		simulation.minimizeEnergy()

		TotalPotEnergy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
		print ("After minimization number {}, the energy is {} kcal/mole".format(RestraintIx,TotalPotEnergy))

		if (args.Verbose):
			SaveState = simulation.context.getState(getPositions=True)
			pdbWrite = pdbreporter.PDBReporter("{}_RestraintProfile{}_min.pdb".format(args.OutputRoot, RestraintIx), reportInterval=0)
			pdbWrite.report(simulation, SaveState)

		##We remove the restraint, so a new one can take its place
		system.removeForce(system.getNumForces()-1)


##And at last we do a global minimization, or if no
##restrained atoms were specified, we do this from the start

CurrentPositions = simulation.context.getState(getPositions=True).getPositions()
simulation.context.reinitialize(preserveState=True)

simulation.minimizeEnergy()

TotalPotEnergy = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(kilocalories_per_mole)
print ("After global minimization, the energy is {} kcal/mole".format(TotalPotEnergy))

SaveState = simulation.context.getState(getPositions=True)
pdbWrite = pdbreporter.PDBReporter("{}_Global_min.pdb".format(args.OutputRoot), reportInterval=0)
pdbWrite.report(simulation, SaveState)

simulation.saveState("{}_FinalMinimization.xml".format(args.OutputRoot))

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy as np
import datetime

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--files", type=str, help="Roots of files you want to Min,Heat,Equi; needs to be formatted like files.prmtop and files.inpcrd")
parser.add_argument("--Solvent", type=str, help="What type of solvent to use. Possible values: explicit OR implicit.")
parser.add_argument("--OutputRoot", type=str, help="A string to be appended to the begining of the DCD,OUT,CHK and XML files generated by openMM")
parser.add_argument("--LoadState", type=str, help="XML file from previous simulation. Optional")
parser.add_argument("--RestrainedAtomsIn", type=str, help="A file containing the indices of the atoms to be restrained. The restraints that make up the restraint profile must be separated by a newline character.")
parser.add_argument("--RestrainedCoMAtomsIn", type=str, help="A file containing the lists of atoms which make up the two groups that will have their Centers of Mass restrained. The two groups have to be separated by a newline character.")
args = parser.parse_args()

if (version.short_version <= "7.2"):
	print ("Your version of OpenMM ({}) doesn't have the latest form of the 'reinitializeContext' method." \
"Please update to versions higher than 7.2".format(version.short_version))

if (args.OutputRoot is None):
	args.OutputRoot = args.files.split("/")[-1]

def GenerateCoMRestraint(Group1, Group2):

	"""A function that generates an OpenMM Force object, that applies
    an external force between the centers of groups of particles.
    This particular implementation was written for 2 groups, and with
    a flatbottom restraint in mind, but can be easily modified.
        
	Parameters
   ----------
	Group1: list of ints
        Particle indices belonging to the first group.

    Group2: list of ints
        Particle indices belonging to the second group.
	

	Returns
	------

	GeneratedCoMRestraint: simtk.openmm.openmm.CustomCentroidBondForce
		harmonic potential force object, between the centers of the two groups.
       
       
    Notes
    -----
    The center of each group is considered the center of mass of the particles
    of each group.
    
    The flatbottom restraint is defined as :
        { {[abs((d(x) - tol)) + (d(x) - tol)]/2}^2 } * k

	"""

	CustomCoMForce = CustomCentroidBondForce(2, "(((abs(abs(distance(g1,g2)) - tol) + (abs(distance(g1,g2)) - tol))/2)^2) * k")
	P1 = CustomCoMForce.addGlobalParameter("k", (100)*kilocalories_per_mole/nanometer**2)
	P2 = CustomCoMForce.addGlobalParameter("tol", 2*nanometer)

	G1 = CustomCoMForce.addGroup(Group1)
	G2 = CustomCoMForce.addGroup(Group2)
    
	CustomCoMForce.addBond((G1,G2))

	print ("A new Center of Mass restraint was generated, containing {} atoms.".format(len(Group1) + len(Group2)))

	return CustomCoMForce
 
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
if (args.Solvent.lower() == "explicit"):
	system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
	system.addForce(MonteCarloBarostat(p, T))
elif (args.Solvent.lower() == "implicit"):
	system = prmtop.createSystem(implicitSolvent=OBC2, soluteDielectric=1.0, solventDielectric=80.0, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1.2*nanometer, constraints=HBonds, implicitSolventSaltConc=0.15*moles/liter)
else:
	print ("No valid solvent type was selected. Exiting...")
	sys.exit()

integrator = LangevinIntegrator(T, 1/picosecond, timestep)
simulation = Simulation(prmtop.topology, system, integrator)

##We assume you want to keep the centers of mass restrained through
##the whole MD run, so we initialize that force
##before any other forces and keep it throughout.

if (args.RestrainedCoMAtomsIn is not None):
	FlexList  = open(args.RestrainedCoMAtomsIn) 
	LoopAtoms = FlexList.read().split("\n") 
	LoopAtoms = [x for x in LoopAtoms if x]  ##Remove empty elements 

	for RestraintIx in range(len(LoopAtoms)):
		LoopAtoms[RestraintIx] = LoopAtoms[RestraintIx].split()   ##Split into n groups
		LoopAtoms[RestraintIx] = [int(x) for x in LoopAtoms[RestraintIx] if x]  ##Turn to ints
    
	system.addForce(GenerateCoMRestraint(LoopAtoms[0], LoopAtoms[1]))
	print ("Centers of Mass have been restrained!")

	simulation.context.reinitialize(preserveState=True)

# Load previous positions and velocities
# First we need to instantiate the force object
# so the spring constant found in the constrained
# minimization done before can be assigned

system.addForce(GenerateRestraint([],[]))
simulation.context.reinitialize(preserveState=True)
simulation.loadState(args.LoadState)

# We then remove the empty force object
system.removeForce(system.getNumForces()-1)

##	HARMONICALLY RESTRAINED PRODUCTION + CoM ##

if (args.RestrainedAtomsIn is not None):

# We iterate through all the sets of restraints, creating a list
# of all the index groups


	nofSteps = 25000000 #(50 ns) 
	dcdFreq  = 25000 #(1000 frames) 
	outFreq  = 500 #(50000 lines)
	rstFreq  = 10000
	
	FlexList  = open(args.RestrainedAtomsIn) 
	LoopAtoms = FlexList.read().split("\n") 
	LoopAtoms = [x for x in LoopAtoms if x]  ##Remove empty elements 
	
	for RestraintIx in range(len(LoopAtoms)):
		LoopAtoms[RestraintIx] = LoopAtoms[RestraintIx].split()   ##Split into n groups
		LoopAtoms[RestraintIx] = [int(x) for x in LoopAtoms[RestraintIx] if x]  ##Turn to ints
	
	
	print ("There are {} restraint profiles:".format(len(LoopAtoms)))
	for RestraintIx in range(len(LoopAtoms)):
		print ("\t * Restraint Set {}: {} atoms restrained out of a total of {} atoms.".format(RestraintIx, len(LoopAtoms[RestraintIx]),system.getNumParticles()))
	
	##We add reporters
	for RestraintIx in range(len(LoopAtoms)):
		simulation.reporters = []
		simulation.reporters.append(StateDataReporter("{}_ProductionRestraintProfile{}.out".format(args.OutputRoot, RestraintIx), 1000, step=True,
												potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True,
												density=True, separator="\t"))
		simulation.reporters.append(DCDReporter("{}_ProductionRestraintProfile{}.dcd".format(args.OutputRoot, RestraintIx), 6000)) ##200 frames/restraint profile
		simulation.reporters.append(CheckpointReporter("{}_Prod_01.chk".format(args.OutputRoot), rstFreq))
	
		CurrentPositions = simulation.context.getState(getPositions=True).getPositions()
		AtomIndices = LoopAtoms[RestraintIx]
	
	##Add Restraint
		CurrentRestraint = GenerateRestraint(CurrentPositions, AtomIndices)
	
		system.addForce(CurrentRestraint)	
		simulation.context.reinitialize(preserveState=True)
	
		simulation.step(nofSteps)

	
		SaveState = simulation.context.getState(getPositions=True)
		pdbWrite = pdbreporter.PDBReporter("{}_PostProd_RestraintSet{}.pdb".format(args.OutputRoot, RestraintIx), reportInterval=0)
		pdbWrite.report(simulation, SaveState)

	##Remove Restraint
		system.removeForce(system.getNumForces()-1)

##	END HARMONIX ##

##	CoM-ONLY PRODUCTION ##

nofSteps = 250000000 #(500 ns)
dcdFreq  = 25000 #(10000 frames) 
outFreq  = 500 #(500000 lines)
rst1Freq  = 500000
rst2Freq  = 750000
print("Production run for " + str(nofSteps * timestep))
print("Simulation started at " + str(datetime.datetime.now()))

simulation.reporters.append(StateDataReporter("{}_Prod.out".format(args.OutputRoot), outFreq, step=True,
                              potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True,
                              density=True, totalSteps=nofSteps, remainingTime=True, speed=True, progress=True, separator="\t"))
simulation.reporters.append(DCDReporter("{}_Prod.dcd".format(args.OutputRoot), dcdFreq)) 
simulation.reporters.append(CheckpointReporter("{}_Prod_01.chk".format(args.OutputRoot), rst1Freq))
simulation.reporters.append(CheckpointReporter("{}_Prod_02.chk".format(args.OutputRoot), rst2Freq))

simulation.step(nofSteps)
simulation.saveState(args.OutputRoot + 'Final.xml') #We save it, since we may want to continue this simulation.
simulation.saveCheckpoint(args.OutputRoot + "Final.chk") #We also save a checkpoint, if we want to continue this simulation later.
print("Done at " + str(datetime.datetime.now()))

## END CoM-ONLY PRODUCTION ##

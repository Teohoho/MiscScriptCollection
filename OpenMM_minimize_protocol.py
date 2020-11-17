#this script runs minimizations in the following workflow:
#first it minimizes the ALL sidechains, then sidechains
#and loops (backbone) and then the whole system
#The selections are done using atom indices, obtained using
#MDTraj

#TODO: add a flag for explicit solvent

from simtk.openmm.app import *
from simtk.openmm import *
import simtk.unit
import mdtraj as md
import sys
import datetime


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--files", type=str, nargs='+', help="Roots of files you want to Min,Heat,Equi; needs to be formatted like files.prmtop and files.inpcrd")
parser.add_argument("--loopIDs", type=str, nargs='+', help="What ResIDs are found in Loops. If not given, only global minimization wil be done. IF you want to run succesive minimization rounds on different sets of loops, add them in the order you want to minimize them. Must be MDTraj DSL language.")
args = parser.parse_args()

if (args.loopIDs is None):
    print ("Loop IDs not supplied. Global minimization will be done.")

#Generate the list of files, given that you want to run this script on files that have the same root
lof = args.files
print("Running protocol on: ", end='')
print(lof)

# Read files
for Molecule in lof:
    prmtop = AmberPrmtopFile(Molecule + ".prmtop")
    inpcrd = AmberInpcrdFile(Molecule + ".inpcrd")
    args.preffix = str(Molecule).split("/")[-1] #This is for naming purposes.

# Read into MDTraj as well

    MDTrajTrajectoryObject = md.load(Molecule + ".inpcrd", top = Molecule + ".prmtop")

# Thermodynamics
    T = 0*simtk.unit.kelvin
    p = 1*simtk.unit.bar
    timestep = 0.000*simtk.unit.picoseconds

# Setup the system (Explicit Solvent)
#    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*simtk.unit.nanometer, constraints=HBonds)
#    system.addForce(MonteCarloBarostat(p, T))

# Setup the system (Implicit Solvent)
    system = prmtop.createSystem(implicitSolvent=OBC2, soluteDielectric=1.0, solventDielectric=80.0, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1.2*simtk.unit.nanometer, constraints=None, implicitSolventSaltConc=0.15*simtk.unit.moles/simtk.unit.liter) 

    integrator = LangevinIntegrator(T, 1/simtk.unit.picosecond, timestep)
    simulation = Simulation(prmtop.topology, system, integrator, platformProperties={'DisablePmeStream':'true'})
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    if (args.loopIDs is not None):
        ##We extract atom indices using mdtraj
        BackboneAtoms = MDTrajTrajectoryObject.topology.select("backbone")
        LoopAtoms = len(args.loopIDs) * [None]
        for LoopIx in range (len(LoopAtoms)):
            LoopAtoms[LoopIx] = MDTrajTrajectoryObject.topology.select("not " + args.loopIDs[LoopIx])   #Note: we use "not" because we want these to be mobile, not constrained

    energy_PM = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(simtk.unit.kilocalories_per_mole)
    print ("Energy before minimization is {}  Kcal/mole".format(energy_PM), end='\n\n')

# We need to store the "real" mass of the atoms so we can
# restore it later

    ActualMasses = []

    for PartIx in range(system.getNumParticles()):
        ActualMasses.append(system.getParticleMass(PartIx).value_in_unit(simtk.unit.dalton))

    if (args.loopIDs is not None):

# Minimization of sidechains

        for PartIx in BackboneAtoms:
            system.setParticleMass(int(PartIx), 0)
   
        print("Minimizing {} sidechains…".format(Molecule))
        simulation.minimizeEnergy()
        print("Minimization of {} sidechains done.".format(Molecule), end='\n\n')

        SaveState = simulation.context.getState(getPositions=True)
        pdbWrite = pdbreporter.PDBReporter("{}_SC_min.pdb".format(Molecule), reportInterval=0)
        pdbWrite.report(simulation, SaveState)

    
        energy_PM = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(simtk.unit.kilocalories_per_mole)
        print ("Energy after sidechain minimization is {} Kcal/mole".format(energy_PM))
    
        #Restore masses:
        for PartIx in range(system.getNumParticles()):
            system.setParticleMass(PartIx, ActualMasses[PartIx])

# Minimization of all Loops given

        for LoopIx in range(len(LoopAtoms)):
            for PartIx in LoopAtoms[LoopIx]:
                system.setParticleMass(int(PartIx), 0)
    
            print("Minimizing {} (Loop number {})…".format(Molecule, LoopIx))
            simulation.minimizeEnergy()
            energy_PM = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(simtk.unit.kilocalories_per_mole)
            print ("Energy after loop #{} minimization is {} Kcal/mole".format(LoopIx,energy_PM))
   
            SaveState = simulation.context.getState(getPositions=True)
            pdbWrite = pdbreporter.PDBReporter("{}_Loop{}_min.pdb".format(Molecule,LoopIx), reportInterval=0)
            pdbWrite.report(simulation, SaveState)

    
            ##Remember to restore the masses:
            for PartIx in range(system.getNumParticles()):
                system.setParticleMass(PartIx, ActualMasses[PartIx])


# Finally, we run a global minimization.

    print("Minimizing {} globally…".format(Molecule))
    simulation.minimizeEnergy()

    SaveState = simulation.context.getState(getPositions=True)
    pdbWrite = pdbreporter.PDBReporter("{}_All_min.pdb".format(Molecule), reportInterval=0)
    pdbWrite.report(simulation, SaveState)

    print("Minimized {} globally…".format(Molecule))
    energy_PM = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(simtk.unit.kilocalories_per_mole)
    print ("Energy after global minimization is {} Kcal/mole".format(energy_PM))

#this script runs minimizations in the following workflow:
#first it minimizes the ALL sidechains, then sidechains
#and loops (backbone) and then the whole system
#The selections are done using atom indices, obtained using
#MDTraj

from simtk.openmm.app import *
from simtk.openmm import *
import simtk.unit
import mdtraj as md
import sys
import datetime


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--files", type=str, nargs='+', help="Roots of files you want to Min,Heat,Equi; needs to be formatted like files.prmtop and files.inpcrd")
parser.add_argument("--FixedIDs", type=str, nargs='+', help="What atoms to keep constrained. If not given, only global minimization wil be done. IF you want to run succesive minimization rounds on different sets of loops, add them in the order you want to minimize them. Must be MDTraj DSL language.")
parser.add_argument("--solvent", type=str, help="What kind of solvent to use? Possible arguments: explicit | implicit")
args = parser.parse_args()

if (args.FixedIDs is None):
    print ("Fixed IDs not supplied. Global minimization will be done.")

#Generate the list of files, given that you want to run this script on files that have the same root
lof = args.files
print("Running protocol on: ", end='')
print(lof)

# Read files
for Molecule in lof:
    prmtop = AmberPrmtopFile(Molecule + ".prmtop")
    inpcrd = AmberInpcrdFile(Molecule + ".inpcrd")

# Read into MDTraj as well

    MDTrajTrajectoryObject = md.load(Molecule + ".inpcrd", top = Molecule + ".prmtop")

# Thermodynamics
    T = 0*simtk.unit.kelvin
    p = 1*simtk.unit.bar
    timestep = 0.000*simtk.unit.picoseconds

# Setup the system (Implicit solvent)
    if (args.solvent == "implicit"):
        system = prmtop.createSystem(implicitSolvent=OBC2, soluteDielectric=1.0, solventDielectric=80.0, nonbondedMethod=CutoffNonPeriodic, nonbondedCutoff=1.2*nanometer, constraints=HBonds, implicitSolventSaltConc=0.15*moles/liter)

#Setup the system (Explicit Solvent)

    elif (args.solvent == "explicit"):
        system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
        system.addForce(MonteCarloBarostat(p, T))

#in case the user selects an invalid solvent type

    else:
        print ("Solvent type not recognized. Please choose 'explicit' or 'implicit' and try again.")
        sys.exit(0)


    integrator = LangevinIntegrator(T, 1/simtk.unit.picosecond, timestep)
    simulation = Simulation(prmtop.topology, system, integrator, platformProperties={'DisablePmeStream':'true'})
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    if (args.FixedIDs is not None):
        ##We extract atom indices using mdtraj
        BackboneAtoms = MDTrajTrajectoryObject.topology.select("backbone")
        FixedAtoms = len(args.FixedIDs) * [None]
        for FixIx in range (len(FixedAtoms)):
            FixedAtoms[FixedIx] = MDTrajTrajectoryObject.topology.select("not ({})".format(args.FixedIDs[FixedIx]))   #Note: we use "not" because we want these to be mobile, not constrained

    energy_PM = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(simtk.unit.kilocalories_per_mole)
    print ("Energy before minimization is {}  Kcal/mole".format(energy_PM), end='\n\n')

# We need to store the "real" mass of the atoms so we can
# restore it later

    ActualMasses = []

    for PartIx in range(system.getNumParticles()):
        ActualMasses.append(system.getParticleMass(PartIx).value_in_unit(simtk.unit.dalton))

    if (args.FixedIDs is not None):

# Minimization of all Fixed Domains given

        for FixedIx in range(len(FixedAtoms)):
            for PartIx in FixedAtoms[FixedIx]:
                system.setParticleMass(int(PartIx), 0)
    
            print("Minimizing {} (Fixed number {})…".format(Molecule, FixedIx))
            simulation.minimizeEnergy()
            energy_PM = simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(simtk.unit.kilocalories_per_mole)
            print ("Energy after loop #{} minimization is {} Kcal/mole".format(FixedIx,energy_PM))
   
            SaveState = simulation.context.getState(getPositions=True)
            pdbWrite = pdbreporter.PDBReporter("{}_FixedDomain{}_min.pdb".format(Molecule,FixedIx), reportInterval=0)
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

import mdtraj as md
import numpy as np
import scipy.stats

import datetime

import argparse
parser = argparse.ArgumentParser()

parser.add_argument("--TrajIn", type=str, help="Trajectory File In")
parser.add_argument("--TopIn" , type=str, help="Topology File In")
parser.add_argument("--ScoreFile" , type=str, help="Output scores for each frame to this file.")
parser.add_argument("--ReceptorSele" , type=str, help="A DSL selection for your receptor")
parser.add_argument("--LigandSele" , type=str, help="A DSL selection for your ligand")

args = parser.parse_args()

## Load Trajectory
MDTrajTrajectoryObject = md.load(args.TrajIn, top=args.TopIn)
MDTrajTrajectoryObject = MDTrajTrajectoryObject[0:50]
print ("Trajectory loaded: {} frames".format(MDTrajTrajectoryObject.n_frames))

## Define Ligand and Receptor selections, in DSL selection language

LigandSele = MDTrajTrajectoryObject.topology.select(args.LigandSele)
ReceptorSele = MDTrajTrajectoryObject.topology.select(args.ReceptorSele)

###################
## CONTACT SCORE ##
###################

StartTime = datetime.datetime.now()

## Define indices corresponding to charged residues

AcidSele  = MDTrajTrajectoryObject.topology.select("(resname GLU and name CD) or (resname ASP and name CG)")
BasicSele = MDTrajTrajectoryObject.topology.select("(resname ARG and name CZ) or (resname LYS and name CE)")

## Generate lists for neighbors

AcidReceptor  = np.intersect1d(ReceptorSele, AcidSele)
BasicReceptor = np.intersect1d(ReceptorSele, BasicSele)
AcidLigand    = np.intersect1d(LigandSele, AcidSele)
BasicLigand   = np.intersect1d(LigandSele, BasicSele)

## Define neighbor cutoff, in nm

cutoff = 0.7

## Compute neighbors for each pair of charges; 
## Nomenclature explanation: LARA = Ligand Acid Receptor Acid

LARA = md.compute_neighbors(MDTrajTrajectoryObject,cutoff=cutoff,query_indices=AcidLigand,haystack_indices=AcidReceptor)
LARB = md.compute_neighbors(MDTrajTrajectoryObject,cutoff=cutoff,query_indices=AcidLigand,haystack_indices=BasicReceptor)
LBRA = md.compute_neighbors(MDTrajTrajectoryObject,cutoff=cutoff,query_indices=BasicLigand,haystack_indices=AcidReceptor)
LBRB = md.compute_neighbors(MDTrajTrajectoryObject,cutoff=cutoff,query_indices=BasicLigand,haystack_indices=BasicReceptor)
    
## We generate a numpy array to store the charge scores

ContactScore = np.empty(MDTrajTrajectoryObject.n_frames, dtype="int64")
    
## We compute and add the counts from each frame in the array

for FrameIx in range(MDTrajTrajectoryObject.n_frames):
    FrameCount = (len(LARB[FrameIx]) + len(LBRA[FrameIx])) - (len(LBRB[FrameIx]) + (len(LARA[FrameIx])))
    ContactScore[FrameIx] = FrameCount
    
elapsedTime = datetime.datetime.now() - StartTime

print ("Contacts score computed in {} seconds for {} frames".format(elapsedTime.total_seconds(), MDTrajTrajectoryObject.n_frames))


###################
##   SASA Score  ##
###################

StartTime = datetime.datetime.now()


## Generate separate systems

LigandMDTObject   = MDTrajTrajectoryObject.atom_slice(LigandSele)
ReceptorMDTObject = MDTrajTrajectoryObject.atom_slice(ReceptorSele)

## Compute SASAs

ComplexSASA  = md.shrake_rupley(MDTrajTrajectoryObject)
LigandSASA   = md.shrake_rupley(LigandMDTObject)
ReceptorSASA = md.shrake_rupley(ReceptorMDTObject)

## Generate SASA array

SASAScore = np.full(MDTrajTrajectoryObject.n_frames, -1.0)

## Compute SASA difference, and assign it in the array

for FrameIx in range(MDTrajTrajectoryObject.n_frames):
    FrameScore = np.sum(ReceptorSASA[FrameIx]) + np.sum(LigandSASA[FrameIx]) - np.sum(ComplexSASA[FrameIx])
    SASAScore[FrameIx] = FrameScore
    
#print (SASAScore)

elapsedTime = datetime.datetime.now() - StartTime

print ("SASA computed in {} seconds for {} frames".format(elapsedTime.total_seconds(), MDTrajTrajectoryObject.n_frames))


###################
##   Score File  ##
###################

ScoreFileOut = open(args.ScoreFile, "w")

## General Statistics
ScoreFileOut.write("Number of frames: {}\n".format(MDTrajTrajectoryObject.n_frames))
ScoreFileOut.write("Contact score edges: {}\t{}\n\n".format(np.min(ContactScore), np.max(ContactScore)))
ScoreFileOut.write("SASA score edges: {}\t{}\n".format(np.min(SASAScore), np.max(SASAScore)))

## Header
ScoreFileOut.write("{}\t{}\t{}\t{}\n".format("Frame", "Contact ", "SASA Score", "SASA Percentile"))

##Print actual data
for FrameIx in range(MDTrajTrajectoryObject.n_frames):
    #ScoreFileOut.write("{:>6}\t{:>7}\t{:>15.6}\t{:>17}\n".format(FrameIx, ContactScore[FrameIx], SASAScore[FrameIx], scipy.stats.percentileofscore(SASAScore,SASAScore[FrameIx])))
    ScoreFileOut.write("{}\t{:>3}\t\t{:.6}\t\t\t{}\n".format(FrameIx, ContactScore[FrameIx], SASAScore[FrameIx], scipy.stats.percentileofscore(SASAScore,SASAScore[FrameIx])))

ScoreFileOut.close()

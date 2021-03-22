from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sys

class GenerateRestraint:
    """
    A class that makes it easy to generate a wide array
    of restraining force objects, given an OpenMM simulation
    object

    Parameters
    ----------

    OpenMMSim:  OpenMMSimulationObject
                OpenMM Simulation Object where we get the particle
                positions from.

    HarmonicIndices:     str
                Path to a file containing atom indices that will be used
                in the harmonic restraint.
    
    FlatbottomIndices:     str
                Path to a file containing atom indices that will be later
                used in the following force definitions.
    """
    def __init__(self, OpenMMSim, HarmonicIndices, FlatbottomIndices):
        OpenMMSim = self.OpenMMSim

        ## Read HarmonicIndices file, and generate a list of indices.
        HarmonicList  = open(args.RestrainedCoMAtomsIn).readlines().split("\n")
	HarmonicList = [x for x in HarmonicList if x]  ##Remove empty elements 

	for RestraintIx in range(len(HarmonicList)):
	    HarmonicList[RestraintIx] = HarmonicList[RestraintIx].split()   ##Split into n groups
	    HarmonicList[RestraintIx] = [int(x) for x in HarmonicList[RestraintIx] if x]  ##Turn to ints

        for Ix in HarmonicList:
            print ("Harmonic restraint containing {} atoms detected.".format(len(Ix)))
            

        FlatbottomList  = open(args.RestrainedCoMAtomsIn).readlines().split("\n")
	FlatbottomList = [x for x in FlatbottomList if x]  ##Remove empty elements 

	for RestraintIx in range(len(FlatbottomList)):
	    FlatbottomList[RestraintIx] = FlatbottomList[RestraintIx].split()   ##Split into n groups
	    FlatbottomList[RestraintIx] = [int(x) for x in FlatbottomList[RestraintIx] if x]  ##Turn to ints

        if (len(FlatbottomList) != 2):
            print ("Warning! The Flatbottom restraint file needs to have exactly 2 lines. (got {} lines)!".format(len(FlatbottomList)))
            sys.exit()

    def HarmonicRestraint(self, ParticleIx, K=10, dummy=False):
        """
        Generate a Harmonic, atom position-based restraint.

        Parameters
        ----------

        ParticleIx:     list of Ints
                        Indices of particles affected by the harmonic restraint.

        K:              float
                        Value of the spring constant for the harmonic restraint. Given in 
                        units of kcal/mol/nm^2

        dummy:          bool
                        When getting state info from a xml file, OpenMM needs
                        to have all the parameters from said XML defined in the simulation
                        object in which the XML is loaded. Use this to generate an empty force 
                        object, so you can run LoadState on an XML where this force was used.

        Returns
        -------

        HarmonicForceObject:    CustomExternalForce Object
                                harmonic potential force object
        """

        ## Define the equation for the Harmonic Force
        CustomHarmonicForce = CustomExternalForce("springConstant * (periodicdistance(x, y, z, x0, y0, z0)^2)")
        CustomHarmonicForce.addGlobalParameter("springConstant", k*kilocalories_per_mole/nanometer**2)
        CustomHarmonicForce.addPerParticleParameter("x0")
        CustomHarmonicForce.addPerParticleParameter("y0")
        CustomHarmonicForce.addPerParticleParameter("z0")

        if (dummy==False):
        ## Get positions from the OpenMMSim
            Positions = self.OpenMMSim.context.getState(getPositions=True).getPositions()
    
            ## Iterate through the given ParticleIx and add the
            ## proper atoms to the force object
            for PartIx in ParticleIx:
                PartX = Positions[PartIx][0].value_in_unit(nanometer)
                PartY = Positions[PartIx][1].value_in_unit(nanometer)
                PartZ = Positions[PartIx][2].value_in_unit(nanometer)
                
                CustomHarmonicForce.addParticle(PartIx, [PartX, PartY, PartZ])
    
            ## Print some info about the object
            print ("A new harmonic force has been created:\nIs periodic: {}\nNumber of Particles: {}".format(CustomHarmonicForce.usesPeriodicBoundaryConditions(), CustomHarmonicForce.getNumParticles()))
    
        return CustomHarmonicForce

    def FlatBottomRestraint(Group1,Group2, springConstant, tol, dummy=False):

        """A function that generates an OpenMM Force object, that applies
        an external force between the centers of groups of particles.
        This particular implementation was written for 2 groups, and with
        a flatbottom restraint in mind, but can be easily modified.
        
	Parameters
        ----------
	Group1:     list of ints
                    Particle indices belonging to the first group.

        Group2:     list of ints
                    Particle indices belonging to the second group.

        springConstant:     float
                    Value of spring constant for the potential

        tol:        float
                    Width of the actual flat bottom.

        dummy:      bool
                    When getting state info from a xml file, OpenMM needs
                    to have all the parameters from said XML defined in the simulation
                    object in which the XML is loaded. Use this to generate an empty force 
                    object, so you can run LoadState on an XML where this force was used.
	
	Returns
	------
	GeneratedCoMRestraint: simtk.openmm.openmm.CustomCentroidBondForce
		Flatbottom potential force object, between the centers of the two groups.
       
       
        Notes
        -----
        The center of each group is considered the center of mass of the particles
        of each group.
        
        The flatbottom restraint is defined as :
            { {[abs((d(x) - tol)) + (d(x) - tol)]/2}^2 } * k

        """
        
        ## Define the equation for the flatbottom potential
	CustomCoMForce = CustomCentroidBondForce(2, "(((abs(abs(distance(g1,g2)) - tol) + (abs(distance(g1,g2)) - tol))/2)^2) * k")
	P1 = CustomCoMForce.addGlobalParameter("k", (100)*kilocalories_per_mole/nanometer**2)
	P2 = CustomCoMForce.addGlobalParameter("tol", 2*nanometer)

        if (dummy == False):

            ## Define the two groups
            G1 = CustomCoMForce.addGroup(Group1)
            G2 = CustomCoMForce.addGroup(Group2)
            
            CustomCoMForce.addBond((G1,G2))
        
            print ("A new Center of Mass restraint was generated, containing {} atoms.".format(len(Group1) + len(Group2)))

	return CustomCoMForce 

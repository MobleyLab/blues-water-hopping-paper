import parmed
from simtk import unit
import mdtraj
import numpy as np
import sys, traceback
import math
import copy
import random
import os
from openeye.oechem import *

class Move(object):

    """Move provides methods for calculating properties on the
    object 'move' (i.e ligand) being perturbed in the NCMC simulation.
    This is the base Move class.
    Ex.
        from blues.ncmc import Model
        ligand = Model(structure, 'LIG')
        ligand.calculateProperties()
    Attributes
    ----------
    """

    def __init__(self):
        """Initialize the Move object
        Currently empy.
        """

    def initializeSystem(self, system, integrator):
        """If the system or integrator needs to be modified to perform the move
        (ex. adding a force) this method is called during the start
        of the simulation to change the system or integrator to accomodate that.
        Parameters
        ----------
        system : simtk.openmm.System object
            System to be modified.
        integrator : simtk.openmm.Integrator object
            Integrator to be modified.
        Returns
        -------
        system : simtk.openmm.System object
            The modified System object.
        integrator : simtk.openmm.Integrator object
            The modified Integrator object.
        """
        new_sys = system
        new_int = integrator
        return new_sys, new_int

    def beforeMove(self, context):
        """This method is called at the start of the NCMC portion if the
        context needs to be checked or modified before performing the move
        at the halfway point.
        Parameters
        ----------
        context: simtk.openmm.Context object
            Context containing the positions to be moved.
        Returns
        -------
        context: simtk.openmm.Context object
            The same input context, but whose context were changed by this function.
        """
        return context

    def afterMove(self, context):
        """This method is called at the end of the NCMC portion if the
        context needs to be checked or modified before performing the move
        at the halfway point.
        Parameters
        ----------
        context: simtk.openmm.Context object
            Context containing the positions to be moved.
        Returns
        -------
        context: simtk.openmm.Context object
            The same input context, but whose context were changed by this function.
        """
        return context

    def _error(self, context):
        """This method is called if running during NCMC portion results
        in an error. This allows portions of the context, such as the
        context parameters that would not be fixed by just reverting the
        positions/velocities of the context.
        Parameters
        ----------
        context: simtk.openmm.Context object
            Context containing the positions to be moved.
        Returns
        -------
        context: simtk.openmm.Context object
            The same input context, but whose context were changed by this function.
        """

        return context

    def move(self, context):
        return context


class WaterTranslationRotationMove(Move):
    """ Move that translates a random water within a specified radius of the protein's
    center of mass to another point within that radius
    Parameters
    ----------
    structure:
        topology: parmed.Topology
            ParmEd topology object containing atoms of the system.
        water_name: str, optional, default='WAT'
            Residue name of the waters in the system.
        radius: float*unit compatible with simtk.unit.nanometers, optional, default=2.0*unit.nanometers
            Defines the radius within the protein center of mass to choose a water
            and the radius in which to randomly translate that water.
    """

    def __init__(self, structure, water_name='WAT', radius=2.0*unit.nanometers):
        print("START OF INIT")
        #initialize self attributes
        self.radius = radius #
        self.water_name = water_name
        self.water_residues = []
        protein_atoms = []
        self.before_ncmc_check = True
        self.traj = mdtraj.load('1znk_out.pdb')

        residues = list(structure.topology.residues())
        for res in residues:
            if res.name == water_name:
                water_mol = []
                for atom in res.atoms():
                    water_mol.append(atom.index)
                self.water_residues.append(water_mol)
            if res.name == 'F09':
                lig = []
                for atom in res.atoms():
                    lig.append(atom.index)
                protein_atoms.append(lig)
        self.atom_indices = self.water_residues[0]
        self.protein_atoms = protein_atoms[0][9]#9

        print("self.atom_indices",self.atom_indices)
        print("TYPE self.atom_indices", type(self.atom_indices))

        print("structure", structure)
        print("TYPE structure", type(structure))

        print("water_name",water_name)
        print("TYPE water_name",type(water_name))

        print("self.radius",self.radius)
        print("TYPE self.radius",type(self.radius))

        print("self.traj",self.traj)
        print("TYPE self.traj",type(self.traj))

        print("self.water_residues",len(self.water_residues))
        print("TYPE self.water_residues",type(self.water_residues))

        print("self.protein_atoms", self.protein_atoms)
        print("TYPE self.protein_atoms",type(self.protein_atoms))


    def _random_sphere_point(self, radius, origin):
        """function to generate a uniform random point
        in a sphere of a specified radius.
        Used to randomly translate the water molecule

        Parameters
        ----------
        radius: float
            Defines the radius of the sphere in which a point
            will be uniformly randomly generated.
        """
        print("Radius of _random_sphere_point", radius)
        r = radius * ( np.random.random()**(1./3.) )  #r (radius) = specified radius * cubed root of a random number between 0.00 and 0.99999
        phi = np.random.uniform(0,2*np.pi) #restriction of phi (or azimuth angle) is set from 0 to 2pi. random.uniform allows the values to be chosen w/ an equal probability
        costheta = np.random.uniform(-1,1) #restriction set from -1 to 1
        theta = np.arccos(costheta) #calculate theta, the angle between r and Z axis
        x = np.sin(theta) * np.cos(phi) #x,y,and z are cartesian coordinates
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        sphere_point = np.array([x, y, z]) * r + origin
        return sphere_point

    def getMasses(self, topology):
        """Returns a list of masses of the specified ligand atoms.
        Parameters
        ----------
        topology: parmed.Topology
            ParmEd topology object containing atoms of the system.
        """
        print('This is the start of the getMasses function....')
        masses = unit.Quantity(np.zeros([int(topology.getNumAtoms()),1],np.float32), unit.dalton)
        for idx,atom in enumerate(topology.atoms()):
            masses[idx] = atom.element._mass #gets the mass of the atom, adds to list (along with index)
        print('This is the end of the getMasses function....')
        return masses

    def getCenterOfMass(self, positions, masses):
        """Returns the calculated center of mass of the ligand as a np.array
        Parameters
        ----------
        positions: parmed.Structure
            ParmEd positions of the atoms to be moved.
        masses : numpy.array
            np.array of particle masses
        """
        coordinates = np.asarray(positions._value, np.float32) #gives the value of atomic positions as an array
        center_of_mass = parmed.geometry.center_of_mass(coordinates, masses) * positions.unit
        print('This is the end of the getCOM function....')
        return center_of_mass

    def beforeMove(self, nc_context):
        """
        Temporary fterHop/unction (until multiple alchemical regions are supported),
        which is performed at the beginning of a ncmc iteration. Selects
        a random water within self.radius of the protein's center of mass
        and switches the positions and velocities with the alchemical water
        defined by self.atom_indices, effecitvely duplicating mulitple
        alchemical region support.
        Parameters
        ----------
        nc_context: simtk.openmm Context object
            The context which corresponds to the NCMC simulation.
        """
        print('This is the start of the beforeMove function....')
        start_state = nc_context.getState(getPositions=True, getVelocities=True)
        start_pos = start_state.getPositions(asNumpy=True) #gets starting positions
        #print('start_pos', start_pos[self.atom_indices[0]]) #prints starting position of the first water atom
        start_vel = start_state.getVelocities(asNumpy=True) #gets starting velocities
        print("start_vel",start_vel)
        switch_pos = np.copy(start_pos)*start_pos.unit
        switch_vel = np.copy(start_vel)*start_vel.unit
        print('ALCH before_switch', switch_pos[self.atom_indices[0]])
        print("Looking for a random water inside the radius...")
        is_inside_sphere = False
        #TODO use random.shuffle to pick random particles (limits upper bound)
        while not is_inside_sphere:
            #water_index = np.random.choice(range(1, 4))
            water_index = np.random.randint(2,8679)

            print("water_index", water_index)

            water_choice = self.water_residues[water_index]
            print("water_choice", water_choice)
            print("water_choice type", type(water_choice))
            print("start_pos[water_choice[:]]", start_pos[water_choice[:]])

            print("start_pos",start_pos)
            print("start_pos", np.shape(start_pos))
            print("self.traj.xyz[0,:,:] BEFORE", self.traj.xyz[0,:,:])
            print("self.traj.xyz", np.shape(self.traj.xyz[0,:,:]))

            self.traj.xyz[0,:,:] = start_pos;
            #print("self.traj.xyz[0,:,:] after", self.traj.xyz[0,:,:])
            start_pos2 = start_state.getPositions(asNumpy=True)
            print("current POS (should match ones above)", start_pos2)
            print("np.array(water_choice[0]).flatten()",np.array(water_choice[0]).flatten())

            #print("np.array(self.protein_atoms[0]).flatten()",np.array(self.protein_atoms[0]).flatten())
            print("np.array(self.protein_atoms).flatten()",np.array(self.protein_atoms).flatten())

            #pairs = self.traj.topology.select_pairs(np.array(water_choice[0]).flatten(), np.array(self.protein_atoms[0]).flatten())
            pairs = self.traj.topology.select_pairs(np.array(water_choice[0]).flatten(), np.array(self.protein_atoms).flatten())

            water_distance = mdtraj.compute_distances(self.traj, pairs, periodic=True)
            print("water_distance",water_distance)
            water_dist = np.linalg.norm(water_distance)
            print("water_dist",water_dist)

            if water_dist <= (self.radius.value_in_unit(unit.nanometers)):
                is_inside_sphere = True
                print("water is within radius")
            else:
                print("water was outside of the radius")
        print("....a random water has been found inside the radius")
        #replace chosen water's positions/velocities with alchemical water
        for i in range(3):
            print("Now replacing the chosen water's positions & velocities w/ alchemical water")
            #set indices of the alchemical waters atoms equal to the indices of the starting positions of the random waters atoms
            switch_pos[self.atom_indices[i]] = start_pos[water_choice[i]]
            #do the same for velocity
            switch_vel[self.atom_indices[i]] = start_vel[water_choice[i]]
            print("Now replacing the alch. water's positions & velocities w/ the random water")
            #set indices of the randomly chosen waters atom equal to alchemical waters atom indices. Same w/ velocity
            switch_pos[water_choice[i]] = start_pos[self.atom_indices[i]]
            switch_vel[water_choice[i]] = start_vel[self.atom_indices[i]]
            #print("Velocities and positions have been switched")
        print('ALCH after switch', switch_pos[self.atom_indices[0]]) #prints the new indices of the alchemical water
        nc_context.setPositions(switch_pos)
        nc_context.setVelocities(switch_vel)
        print("This is the end of the beforeMove function...")
        return nc_context

    def move(self, context):
        """
        This function is called by the blues.MoveEngine object during a simulation.
        Translates the alchemical water randomly within a sphere of self.radius.
        """
        print("START OF MOVE FUNCTION")
        #get the position of the system from the context
        #print("THIS IS THE MOVE FUNCTION")
        before_move_pos = context.getState(getPositions=True).getPositions(asNumpy=True)
        #protein_pos = before_move_pos[self.protein_atoms[0]] #gets the positions from the indices of the atoms in the protein residues in relation to the system
        protein_pos = before_move_pos[self.protein_atoms] #gets the positions from the indices of the atoms in the protein residues in relation to the system

        #print("protein_pos",protein_pos)#find the center of mass and the displacement
        #origin = before_move_pos[self.protein_atoms[0]]
        origin = before_move_pos[self.protein_atoms]

        #Generate uniform random point in a sphere of a specified radius
        point_within_sphere = self._random_sphere_point(self.radius, origin)
        movePos = np.copy(before_move_pos)*before_move_pos.unit #makes a copy of the position of the system from the context
        self.traj.xyz[0,:,:] = movePos;
        #pairs = self.traj.topology.select_pairs(np.array(self.atom_indices[0]).flatten(), np.array(self.protein_atoms[0]).flatten())
        pairs = self.traj.topology.select_pairs(np.array(self.atom_indices[0]).flatten(), np.array(self.protein_atoms).flatten())

        water_distance = mdtraj.compute_distances(self.traj, pairs, periodic=True)

        print("point_within_sphere", point_within_sphere)

        for index, resnum in enumerate(self.atom_indices):
            print('before', before_move_pos[resnum])

        print("self.atom_indices", self.atom_indices)
        # Replace the oxygens coordinates to the new location
        movePos[self.atom_indices[0]] = point_within_sphere
        # Get vector distance of the two hydrogens from the oxygen
        H1V = movePos[self.atom_indices[1]] - before_move_pos[self.atom_indices[0]]
        H2V = movePos[self.atom_indices[2]] - before_move_pos[self.atom_indices[0]]
        # Move the hydrogens to new location
        # self.atom_indices[1] = hydrogen1 and self.atom_indices[0] = oxygen
        movePos[self.atom_indices[1]] = H1V + point_within_sphere
        movePos[self.atom_indices[2]] = H2V + point_within_sphere
        for index, resnum in enumerate(self.atom_indices):
            print('AFTER', movePos[resnum])
        context.setPositions(movePos)
        return context

    def afterMove(self, nca_context):
        """This method is called at the end of the NCMC portion if the
        context needs to be checked or modified before performing the move
        at the halfway point.
        Parameters
        ----------
        context: simtk.openmm.Context object
            Context containing the positions to be moved.
        Returns
        -------
        context: simtk.openmm.Context object
            The same input context, but whose context were changed by this function.
        """
        print("Beginning of the afterMove")

        before_final_move_pos = nca_context.getState(getPositions=True).getPositions(asNumpy=True)

        #Update positions for distance calculation
        print("before_final_move_pos", before_final_move_pos)
        movePos_a = np.copy(before_final_move_pos)*before_final_move_pos.unit
        print("movePos_a",movePos_a)
        self.traj.xyz[0,:,:] = movePos_a;
        #print("self.traj.xyz[0,:,:]",self.traj.xyz[0,:,:])

        #pairs = self.traj.topology.select_pairs(np.array(self.atom_indices[0]).flatten(), np.array(self.protein_atoms[0]).flatten())
        pairs = self.traj.topology.select_pairs(np.array(self.atom_indices[0]).flatten(), np.array(self.protein_atoms).flatten())

        water_distance = mdtraj.compute_distances(self.traj, pairs, periodic=True)
        print("water_distance",water_distance)
        print("np.linalg.norm(self.radius._value)",np.linalg.norm(self.radius._value))
        if water_distance > np.linalg.norm(self.radius._value):
        #if np.linalg.norm(water_distance) > 1.7*unit.nanometers._value:
            nca_context._integrator.setGlobalVariableByName("protocol_work", 999999)
        print("end of afterMove")
        return nca_context

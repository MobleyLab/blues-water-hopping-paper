from mup_WaterTranslationRotation import WaterTranslationRotationMove
from blues.moves import MoveEngine
from blues.simulation import *
import json
from blues.settings import Settings

opt = Settings('mup.yaml').asDict()
structure = opt['Structure']
restart = parmed.amber.Rst7('/oasis/tscc/scratch/bergazin/water/mup_NPT_MD.rst7')#change to dir where rst is
#restart = parmed.amber.Rst7('/oasis/tscc/scratch/bergazin/water/MUP-1/NCMCsteps_variation_trials/1000_NCMCsteps_trials/t1/t1.rst7')

structure.positions = restart.positions
#structure.velocities = restart.velocities
structure.box = restart.box
print("structure.box",structure.box)
print(json.dumps(opt, sort_keys=True, indent=2, skipkeys=True, default=str))

#Select move type
ligand = WaterTranslationRotationMove(structure, water_name='WAT')

#Iniitialize object that selects movestep
ligand_mover = MoveEngine(ligand)

#Generate the openmm.Systems outside SimulationFactory to allow modifications
systems = SystemFactory(structure, ligand.atom_indices, opt['system'])

#Freeze atoms in the alchemical system
#systems.md = systems.freeze_atoms(structure, systems.md, **opt['freeze'])
#systems.alch = systems.freeze_atoms(structure, systems.alch, **opt['freeze'])

# Restrain atoms in the MD and alchemical system
#systems.md = systems.restrain_positions(structure, systems.md, **opt['restraints'])
#systems.alch = systems.restrain_positions(structure, systems.alch, **opt['restraints'])

#Generate the OpenMM Simulations
simulations = SimulationFactory(systems, ligand_mover, opt['simulation'], opt['md_reporters'], opt['ncmc_reporters'])

#Energy minimize system
#simulations.md.minimizeEnergy(maxIterations=0)
#simulations.md.step(50000)
#state = simulations.md.context.getState(getPositions=True, getEnergy=True)
#print('Minimized energy = {}'.format(state.getPotentialEnergy().in_units_of(unit.kilocalorie_per_mole)))

# MD simulation
#simulations.md.step(opt['simulation']['nstepsMD'])

#BLUES Simulation
blues = BLUESSimulation(simulations)
blues.run(**opt['simulation'])
accept_before = 0
accept_after = 0
nIter = 10000
for i in range(nIter):
    blues.run(**opt['simulation'])  #specify in YAML to do 1 iteration
    accept_after = blues.accept
    if accept_after > accept_before:
        sim = blues._ncmc_sim
        topology = sim.topology
        system = sim.context.getSystem()
        state = sim.context.getState(
            getPositions=True,
            getVelocities=True,
            getParameters=True,
            getForces=True,
            getParameterDerivatives=True,
            getEnergy=True,
            enforcePeriodicBox=True)

        # Generate the ParmEd Structure
        structure = parmed.openmm.load_topology(topology, system, xyz=state.getPositions())
        outfname = 't2_%s-acc_mup_unrestrained_15000ncmc.pdb'%(i)
        structure.save(outfname, overwrite=True)
    accept_before = accept_after

# MC Simulation
#mc = MonteCarloSimulation(simulations, opt['simulation'])
#mc.run()

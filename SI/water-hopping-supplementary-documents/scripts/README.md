# Example scripts and used in the calculations

The systems were minimized until forces were below a tolerance of 10 kJ/mol. Long-range electrostatics were calculated using Particle Mesh Ewald. Simulations were run using a 4fs timestep and hydrogen mass repartitioning. The alpha carbons and/or ligands in the systems were restrained with a 5 kcal/(mol A^2) force constant. The temperature was set to 300 K in all cases except the water box with graphene sheets, which was set to 500 K.

Molecular dynamics simulation details: For the Buckyball system, equilibration consisted of 250 ps of NVT MD and 10 ns of NPT MD of equilibration. The MD production run was for 10 ns in the NPT ensemble. For the water box with dividing graphene sheets, equilibration consisted of 5 ns NVT MD. The MUP-1 protein-ligand system was equilibrated for 1 ns of NVT MD and 10ns NPT MD. The MD production run was for 40 ns in the NPT ensemble. The HSP90 protein-ligand system was equilibrated for 1 ns of NVT MD and 80ns NPT MD. The MD production run was for 285 ns in the NPT ensemble.

BLUES simulation details: For the water box system with dividing graphene sheets, BLUES with translational water moves was executed for 240000 iterations, using 2500 NCMC steps and 1000 MD steps. The buckyball system was simulated for a total of 1000 iterations, using 2500 NCMC steps and 1000 MD steps. Additional amounts of NCMC steps were tested for the MUP-1 system, and varied from 1250 to 30000 steps. Both of the solvated MUP-1 and HSP90 (PBD: 5j64) systems were simulated for a total of 10,000 BLUES iterations. For the MUP-1 system, various amounts of NCMC steps were tested, varying from 1250 to 30000 NCMC steps, but the number of MD steps in all cases was 1000 MD steps per iteration. For the HSP90 system, each BLUES iteration consisted of 2500 NCMC steps and 1000 MD steps.

For documentation on the BLUES modules see https://github.com/MobleyLab/blues/blob/master/README.md.
For a tutorial on how to use BLUES see https://mobleylab-blues.readthedocs.io/en/latest/tutorial.html.

### Manifest

- example.py: Example python script for running equilibration.
  - Has options for restraining atoms in the MD and alchemical system.
- example.pbs: Example bash script to execute `example.py` on a TSCC cluster.
- example.yaml: Basic YAML script specifying the parameters for the simulations.
 - `simulation: pressure` - sets the simulation to run NPT. When not used, the default is NVT.
 - `nstepsMD`, `nstepsNC` and `nIter` - control the number of MD steps, NCMC steps and the number of iterations used, respectively.
 - `restraints: selection` - AmberMask selection to apply positional restraints and `weight` is the restraint weight for xyz atom restraints in kcal/(mol A^2). Used to restrain all of the alpha carbons, carbons, and/or ligands. The specific restraint selection was '@CA&!(:158)' for MUP-1, ':209|@CA' for HSP90 and '@C' for the graphene sheets in the water box and the C60 buckyball systems.

- WaterTranslation.py: Example python script for the water hopping moves.
  - The move requires the user to indicate an atom and a radius defining a sphere encompassing an area of interest around the atoms position. This script requires user selection of an atom, the residue name of water, a PDB of the system (used to update positions for distance calculations), and a radius.
  The HSP90 system used atom C10 of the ligand (atom number 3282), roughly the center of the ligand.
  The MUP-1 system used atom C9 of the ligand (atom number 2505)

# Example scripts used in the calculations

For documentation on the BLUES modules see https://github.com/MobleyLab/blues/blob/master/README.md.
For a tutorial on how to use BLUES see https://mobleylab-blues.readthedocs.io/en/latest/tutorial.html.
Additional details of individual system setups are available in the paper.

### Manifest

- `example.py`: Example python script for running simulations. Requires `example.yaml` and `WaterTranslation.py`.
  - Equilibration: First edit `example.yaml` to provide the names of input files (under `structure`) and other details (ie specify the number of MD steps with `nstepsMD`). A `.rst7` file will be generated once equilibration finishes.
  - MD production: In `example.py`, uncomment lines 9-11 and use the equilibration `.rst7` file as input and comment out lines 38-40 (which pertain to energy minimization). In `example.yaml`, update `nstepsMD` for production as needed.
  - BLUES simulation: In `example.py`, uncomment lines 9-11 and use the equilibration `.rst7` file as input and comment out lines 38-40 (which pertain to energy minimization). In `example.yaml`, update the number of `nstepsMD` (number of MD steps to use per BLUES cycle), `nstepsNC` (number of NCMC steps to use per BLUES cycle) and `nIter` (total number of blues cycles).

- `example.pbs`: Example bash script to execute `example.py` on a TSCC cluster.

- `example.yaml`: Basic YAML script specifying the parameters for the simulations.
  - `simulation: pressure` - sets the simulation to run NPT. When not used, the default is NVT.
  - `nstepsMD`, `nstepsNC` and `nIter` - control the number of MD steps, NCMC steps and the number of iterations used, respectively.
  - `restraints: selection` - AmberMask selection to apply positional restraints and `weight` is the restraint weight for xyz atom restraints in kcal/(mol A^2). Used to restrain all of the alpha carbons, carbons, and/or ligands.
    - The specific restraint selection was '@CA&!(:158)' for MUP-1, ':209|@CA' for HSP90 and '@C' for the graphene sheets in the water box and the C60 buckyball systems.

- `WaterTranslation.py`: Example python script to generate the OpenMM Simulations and to perform the water hopping moves. Required by `example.py`.
  - The move requires the user to indicate an atom and a radius defining a sphere encompassing an area of interest around the atoms position. This script requires user selection of an atom (line 142), the residue name of the water in the system (line 116), a PDB of the system (line 124) which is used to update positions for distance calculations, and a radius (line 117). The first water in the system is used as a bookkeeping device (line 140).
    - System atom selection and radius details (more details can be found in the paper):
      - The HSP90 system used atom C10 of the ligand (atom number 3282), roughly the center of the ligand, and the radius was set to 1.5 nanometers.
      - The MUP-1 system used atom C9 of the ligand (atom number 2505) and the radius was set to 2.0 nanometers.
      - The two wall system with dividing graphene sheets used a carbon atom in the middle of one of the sheets and the radius was set to 1.5 nanometers.
      - The buckyball used on of the carbon atoms and the radius was set to 1.2 nanometers.

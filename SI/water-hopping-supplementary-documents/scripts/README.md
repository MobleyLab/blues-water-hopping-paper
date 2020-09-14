# Example scripts used in the calculations

For documentation on the BLUES modules see https://github.com/MobleyLab/blues/blob/master/README.md.

For a tutorial on how to use BLUES see https://mobleylab-blues.readthedocs.io/en/latest/tutorial.html.

Additional details of individual system setups are available in the paper.


## Usage
Steps for reproducing a calculation:
- Equilibration: update [`example.yaml`](example.yaml) with input files and specify the number of MD steps to run in `nstepsMD`. In [`example.py`](example.py), uncomment lines pertaining to minimization and MD. Use [`example.pbs`](example.pbs) to execute [`example.py`](example.py) on a cluster (TSCC in this case).
- BLUES simulation: In [`example.py`](example.py) use the `.rst7` output file from equilibration as input (lines 9-11) and comment out lines which pertain to energy minimization and plain MD. If restraints are used, make sure to uncomment lines 31 and 32 in [`example.py`](example.py) to restrain atoms in both the MD and alchemical system. In [`example.yaml`](example.yaml), update the number of  MD steps to use per BLUES cycle in `nstepsMD`, the number of NCMC steps to use per BLUES cycle in `nstepsNC` and the total number of blues cycles in `nIter`. In [`example.yaml`](example.yaml), use `restraints` to apply positional restraints using AmberMask syntax (some examples can be found [here](https://amber-md.github.io/pytraj/latest/atom_mask_selection.html)). In [`WaterTranslation.py`](WaterTranslation.py) on line 117, input desired radius for defining a sphere and the residue name of the water in the system. Input a PDB of the system on line 124. Input the protein residue name in line 135 and indicate the specific atom that the radius will span from in line 143. Use [`example.pbs`](example.pbs) to execute [`example.py`](example.py) again.


### Manifest

- [`example.py`](example.py): Example python script for running simulations. Requires [`example.yaml`](example.yaml) and [`WaterTranslation.py`](WaterTranslation.py).
  - Equilibration: First edit [`example.yaml`](example.yaml) to provide the names of input files (under `structure`) and other details (ie specify the number of MD steps with `nstepsMD`). A `.rst7` file will be generated once equilibration finishes.
  - MD production: In [`example.py`](example.py), uncomment lines 9-11 and use the equilibration `.rst7` file as input and comment out lines 38-40 (which pertain to energy minimization). In [`example.yaml`](example.yaml`), update `nstepsMD` for production as needed.
  - BLUES simulation: In [`example.py`](example.py), uncomment lines 9-11 and use the equilibration `.rst7` file as input and comment out lines 38-40 (which pertain to energy minimization). In [`example.yaml`](example.yaml), update the number of `nstepsMD` (number of MD steps to use per BLUES cycle), `nstepsNC` (number of NCMC steps to use per BLUES cycle) and `nIter` (total number of blues cycles).

- [`example.pbs`](example.pbs): Example bash script to execute [`example.py`](example.py) on a TSCC cluster.

- [`example.yaml`](example.yaml): Basic YAML script specifying the parameters for the simulations.
  - `simulation: pressure` - sets the simulation to run NPT. When not used, the default is NVT.
  - `nstepsMD`, `nstepsNC` and `nIter` - control the number of MD steps, NCMC steps and the number of iterations used, respectively.
  - `restraints: selection` - AmberMask selection to apply positional restraints and `weight` is the restraint weight for xyz atom restraints in kcal/(mol A^2). Used to restrain all of the alpha carbons, carbons, and/or ligands.
    - The specific restraint selection was '@CA&!(:158)' for MUP-1, ':209|@CA' for HSP90 and '@C' for the graphene sheets in the water box and the C60 buckyball systems.

- [`WaterTranslation.py`](WaterTranslation.py): Example python script to generate the OpenMM Simulations and to perform the water hopping moves. Required by [`example.py`](example.py).
  - The move requires the user to indicate an atom and a radius defining a sphere encompassing an area of interest around the atoms position. This script requires user selection of an atom (line 142), the residue name of the water in the system (line 116), a PDB of the system (line 124) which is used to update positions for distance calculations, and a radius (line 117). The first water in the system is used as a bookkeeping device (line 140).
    - System atom selection and radius details (more details can be found in the paper):
      - The HSP90 system used atom C10 of the ligand (atom number 3282), roughly the center of the ligand, and the radius was set to 1.5 nanometers.
      - The MUP-1 system used atom C9 of the ligand (atom number 2505) and the radius was set to 2.0 nanometers.
      - The two wall system with dividing graphene sheets used a carbon atom in the middle of one of the sheets and the radius was set to 1.5 nanometers.
      - The buckyball used on of the carbon atoms and the radius was set to 1.2 nanometers.

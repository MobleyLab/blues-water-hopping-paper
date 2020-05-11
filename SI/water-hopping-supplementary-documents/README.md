### Manifest
- `SI-tables/`: Folder containing tables as reported in the water darting SI. Tables include acceptance ratios of the replicate simulations at different NCMC step amounts for the MUP-1 system, the average acceptance rate of all moves and average number of force evaluations across replicates for the buried cavity in MUP-1 system to become hydrated, and acceptance ratios of all attempted moves for each replicate simulation of the HSP90 protein-ligand system.
- `scripts/`: Folder containing example scripts used for MD and BLUES simulations.
  - Contains:
    - `example.py`: Example python script for running simulations.
    - `example.pbs`: Example bash script to execute `example.py` on a TSCC cluster.
    - `example.yaml`: Basic YAML script specifying the parameters for the simulations.
    - `WaterTranslation.py`: Example python script for the water hopping moves, also used to generate the OpenMM Simulations.
- `input-files/`: Folder containing the system input files. Contains AMBER files for Major urinary protein 1 (MUP-1), Heat shock protein 90 (HSP90), Buckminsterfullerene C60, and the water box system with dividing graphene walls.

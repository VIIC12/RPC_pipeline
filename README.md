## [[RFD->[PMPNN->Consensus->ColabFold]]] - RPC-Pipeline 
## Logic
```bash
 >>ARRAY
┌──────────────────┐
│>>ROUNDS          │
│    RFdiffusion   │-->RMSD
│ ┌───────▼──────┐ │
│ │>>CYCLES      │ │
│ │  ┌─────┐     │ │
│ │  │     ▼     │ │
│ │  │   PMPNN   │ │
│ │  │     │     │ │
│ │  │     ▼     │ │
│ │  │ Consensus │ │
│ │  │     │     │ │
│ │  │     ▼     │ │
│ │  │ ColabFold │ │-->RMSD
│ │  └─────┤     │ │
│ └────────┼─────┘ │
└──────────┼───────┘
           ▼
       PROJECT_NAME
       |-rfd_op -> RFdiffusion output
         |-PROJECT_ARRAY
       |-pmpnn_ip -> PMPNN input
         |-PROJECT_ARRAY
       |-pmpnn_op -> PMPNN output
         |-PROJECT_ARRAY
       |-cf_op -> ColabFold output
       |-cf_weights -> ColabFold weights
       |-consensus -> Consensus output
       |-PROJECT_NAME_best_cf_model -> Best ColabFold model
       |-PROJECT_NAME_best_af_model -> Best AlphaFold model
       PROJECT_ARRAY_rmsd_results.txt -> RMSD results
```

## _Acknowledgments_
rmsd.py & sec_structure.py (c) Max Beining
calculate_consensus.py (c) T. Schiffner

## _Settings_
#### Number of arrays
`#SBATCH --array=1-2`

#### Environemnt settings
`PMPNN_PATH=/work/ta905ttoo-rfd/ProteinMPNN`<br>
`PMPNN_PYTHON=/home/sc.uni-leipzig.de/ta905ttoo/.conda/envs/mlfold/bin/python`<br>
`PMPNN_ENV=/home/sc.uni-leipzig.de/ta905ttoo/.conda/envs/mlfold/`<br>

`RDIFF_PATH=/home/sc.uni-leipzig.de/mj334kfyi/software/RFdiffusion`<br>
`RFDIFF_ENV=/home/sc.uni-leipzig.de/mj334kfyi/.conda/envs/rfdiff`<br>

#### RFdiffusion settings
`input_pdb="/work/ta905ttoo-rfd/1_rfd/input/input.pdb"`<br>
`contigs="[20-70/B1-15/20-70/0 A1-243/0]"`<br>
`inpaint_seq="[B1-3,B5-7,B10,B13,B15]"`<br>
`length="80-95"`

#### Fixed positions in the binder
#*E.g. KSFINKVPSKR -> K=1, S=2 -> Same numbering as in the PDB file, if the chain starts with 1<br>
`FIX_POS=(4 8 9 11 12 14)`

#### ChainID of the binder (e.g. Antibody) in the input PDB file
`BINDER_CHAIN=A`<br>

#### ChainID of the target (e.g. Antigen) in the input PDB file
`TARGET_CHAIN=B`<br>

#### RMSD treshold
#TODO -> if RMSD > rmsd_treshold, abort<br>
`#rmsd_treshold=5`

#### Sequence of the binder (e.g. Antibody)
`target_seq="ABYAC"`

#### Which strcucture prediction software should be used? AlphaFold or ColabFold
`FOLD=ColabFold`<br>
#! AlphaFold does not work with the current setup<br>

#### Project name 
`PROJECT_NAME="es"`<br>

#### How many RFdiffusion designs should be generated?
`ROUNDS=2`<br>

#### How often should the PMPNN->Consensus->ColabFold cycle be repeated?
`CYLCES=2`<br>

#### PMPNN settings
`chains_to_design="A"`<br>
`model="vanilla"`<br>

v0.9.0
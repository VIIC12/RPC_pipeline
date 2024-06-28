import os
import argparse
import numpy as np
from numpy import array, dot
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.PDB.PDBParser import PDBParser


def parse_args():
    parser = argparse.ArgumentParser(description="Calculate the RMSD between two motifs in a structure, by aligning " + 
                                     "a reference structure first and calculating the RMSD for the motif. Toms want this " +
                                     "for his protein design project. Generally the script assumes that the native " +
                                     "structure contains a chain A to be aligned with an RFdiffusion model chain B.\n\n" +
                                     "The scripts aligns only the backbone atoms (CA, N, C, O) of the specified residues.",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--native", type=str, help="The native pdb file.", required=True)
    parser.add_argument("--native_chain", type=str, help="The chain in the native pdb file to be aligned. Default: A", required=False, default="A")
    parser.add_argument("--residue_chain_native", type=str, help="The chain in the native pdb file to extract residues for RMSD. Default: B", required=False, default="B")
    parser.add_argument("--residue_list_native", type=int, nargs="+", help="The residues to calculate the RMSD in native model.", required=True)
    parser.add_argument("--model", type=str, help="The designed model pdb file.", required=True)
    parser.add_argument("--model_chain", type=str, help="The chain in the model pdb file to be aligned. Default: B", required=False, default="B")
    parser.add_argument("--residue_chain_model", type=str, help="The chain in the designed model pdb file to extract residues for RMSD. Default: A", required=False, default="A")
    parser.add_argument("--residue_list_model", type=int, nargs="+", help="The residues to calculate the RMSD in design model.", required=True)
    
    return parser.parse_args()
    

def align(native, native_chain:str, model, model_chain:str, atom_types = ["CA", "N", "C", "O"]) -> SVDSuperimposer:
    """
    Aligns a model structure onto a native structure. Returns the SVDSuperimposer object to use
    for transforming the model structure.
    """
    native = PDBParser(QUIET=True).get_structure("native", native)
    model = PDBParser(QUIET=True).get_structure("model", model)
    
    for chains in native.get_chains():
        if chains.id == native_chain:
            native_coords = [ a.coord for a in chains.get_atoms() if a.name in atom_types ]
            break
    for chains in model.get_chains():
        if chains.id == model_chain:
            model_coords = [ a.coord for a in chains.get_atoms() if a.name in atom_types ]
            break
        
    if len(native_coords) != len(model_coords):
        raise ValueError("The number of atoms in the native and model structures do not match.")
    
    si = SVDSuperimposer()
    si.set(np.array(native_coords), np.array(model_coords))
    si.run() # Run the SVD alignment
    
    return si


def get_residue_coords(pdbfile:str, chain_id:str, residues:list[int], atom_types:list[str] = ["CA", "N", "C", "O"]) -> np.ndarray:
    """
    Get the coordinates of the specified atom types for the specified residues in the specified chain.
    """
    structure = PDBParser(QUIET=True).get_structure("structure", pdbfile)
    coords = []
    for chain in structure.get_chains():
        if chain.id == chain_id:
            for resi in chain.get_residues():
                if resi.id[1] in residues:
                    coords += [ a.coord for a in resi.get_atoms() if a.name in atom_types ]
            break
        
    if not coords:
        raise ValueError(f"No coordinates found for the specified residues in chain {chain_id}.")
    
    return np.array(coords)


def rmsd_no_align(native_coords:np.array, model_coords:np.array) -> float:
    """
    Calculate the RMSD between two structures without aligning them.
    """
    if len(native_coords) != len(model_coords):
        raise ValueError("The number of atoms in the native and model structures do not match.")
    
    return np.sqrt(np.sum((native_coords - model_coords) ** 2) / native_coords.shape[0])


def main(args):
    superimposer = align(args.native, args.native_chain, args.model, args.model_chain)
    
    native_coords = get_residue_coords(args.native, args.residue_chain_native, args.residue_list_native)
    model_coords = get_residue_coords(args.model, args.residue_chain_model, args.residue_list_model)
    
    # Transform the model coordinates
    rot, tran = superimposer.get_rotran()
    new_model_coords = dot(model_coords, rot) + tran
    print(f"RMSD: {rmsd_no_align(native_coords, new_model_coords)}")
    
if __name__ == '__main__':
    args = parse_args()
    if not os.path.isfile(args.native):
        raise FileNotFoundError(f"File {args.native} not found.")
    if not os.path.isfile(args.model):
        raise FileNotFoundError(f"File {args.model} not found.")
    
    main(args)
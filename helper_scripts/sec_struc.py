#!/usr/bin/env python3

import sys
import pyrosetta
import os
import contextlib

@contextlib.contextmanager
def suppress_stdout():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr

def calc_secstructure(file: str, chain: str = 'A'):
    '''
    Calculate the secondary structure composition of a chain from a pdbfile.
    Returns the secondary structure sequence and the fraction of helix, sheet and loops.
    '''
    with suppress_stdout():
        pyrosetta.init('-out:level 100')  # Set logging level to suppress most of the output
    
    pose = pyrosetta.pose_from_pdb(file)
    Dssp = pyrosetta.rosetta.protocols.moves.DsspMover()
    Dssp.apply(pose)
    
    ss = pose.secstruct()
    ss_seq = [ss[i] for i in range(len(ss)) if pose.pdb_info().pose2pdb(i+1).split()[1] == chain]
    ss_seq = ''.join(ss_seq)
    total = len(ss_seq)
    frac_helix = ss_seq.count('H') / total
    frac_sheet = ss_seq.count('E') / total
    frac_loops = ss_seq.count('L') / total

    return ss_seq, frac_helix, frac_sheet, frac_loops

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: script.py <pdb_file> <chain>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    chain = sys.argv[2]
    
    ss_seq, frac_helix, frac_sheet, frac_loops = calc_secstructure(pdb_file, chain)
    
    print(f"Secondary Structure Sequence: {ss_seq}")
    print(f"Fraction of Helix: {frac_helix:.2f}")
    print(f"Fraction of Sheet: {frac_sheet:.2f}")
    print(f"Fraction of Loops: {frac_loops:.2f}")
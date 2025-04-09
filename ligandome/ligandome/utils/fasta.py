"""Module for working with and producing fasta files."""
from __future__ import annotations

import os
import requests
from pathlib import Path

def fetch_human_proteome_fasta(fasta_outfile: Path) -> None:
    """Fetch a fasta file of the entire human reference proteome.

    Args:
        fasta_outfile (Path): Path to save the fasta file.
    """
    canonical_proteome = requests.get('https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28reviewed%3Atrue%20AND%20proteome%3Aup000005640%29')

    if not os.path.isfile(fasta_outfile):
        with open(fasta_outfile, 'a+') as fasta:
            fasta.write(canonical_proteome.text)

def read_fasta_file(fasta_file: Path) -> dict[str, str]:
    """Read a fasta file as a dictionary.

    Args:
        fasta_file (Path): Path of file to read.

    Returns:
        dict[str, str]: Dictionary of name, sequence.
    """    
    sequences = {}

    filestring = open(fasta_file, 'r+').read().split('>')[1:]
    for entry in filestring:
        lines = entry.split('\n')
        sequences[lines[0]] = ''.join(lines[1:])
    
    return sequences

def peptide_list_to_fasta_file(peptides: list[str], fasta_outfile: Path) -> None:
    """Save a list of peptides as a fasta file with numbered sequence annotations.

    Args:
        peptides (list[str]): Peptides to save.
        fasta_outfile (Path): Fasta file to create and write to.
    """
    
    if os.path.isfile(fasta_outfile):
        os.remove(fasta_outfile)
    
    with open(fasta_outfile, 'a+') as fasta:
        for i, peptide in enumerate(peptides):
            fasta.write(f'>{i}\n')
            fasta.write(f'{peptide}\n')
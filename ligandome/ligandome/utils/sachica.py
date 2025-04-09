"""Module for python wrappers to run sachica algorithm."""
from __future__ import annotations

import os
from tqdm import tqdm
from pathlib import Path
from scipy.sparse import lil_matrix
from ligandome.utils.constants import DATABASE_EXPORTS_PATH

def peptide_list_to_sachica_input(peptides: list[str], output_path: Path) -> None:
    """Write a peptide list to a sachica input file.

    Args:
        peptides (list[str]): List of peptides.
        output_path (Path): Output fasta file to write peptides to.
    """
    with open(output_path, 'a+') as f:
        for p in peptides:
            f.write(p)
            f.write('\n\n')

def run_sachica(peptides_file: Path, edit_distance: int, temp_dir: Path) -> list[tuple[int, int]]:
    """Run sachica algorithm on list of peptides.

    Args:
        peptides_file (Path): File with peptides to compare.
        edit_distance (int): Distance threshold to compare peptides.
        temp_dir (Path): Temporary directory to use for files.

    Returns:
        tuple[int, int]: List of tuples within given edit distance.
    """ 
    if os.path.isfile(f'{str(temp_dir)}/sachica_out_{edit_distance}_{peptides_file.stem}.sachica'):
        os.remove(f'{str(temp_dir)}/sachica_out_{edit_distance}_{peptides_file.stem}.sachica')
    command = f'{str(DATABASE_EXPORTS_PATH.parent)}/third_party_tools/SACHICA/sachica %pdsSCH {peptides_file} 9 {edit_distance} {str(temp_dir)}/sachica_out_{edit_distance}_{peptides_file.stem}.sachica'
    os.system(command)
    print('Reading results...')
    results = open(f'{str(temp_dir)}/sachica_out_{edit_distance}_{peptides_file.stem}.sachica','r').read().split('\n')
    print('Results opened!')
    return results

def sachica_scores_to_sparse_matrix(peptide_edges: list[str], peptide_index_table: dict[str, int]) -> lil_matrix:
    """Parse sachica scores text file into sparse distance matrix.

    Args:
        peptide_edges (list[str]): Sachica output split into lines.
        peptide_index_table (dict[str, int]): Peptide hash table for indexing.

    Returns:
        lil_matrix: Sparse distance matrix of peptide edit distances.
    """
    distance_matrix: lil_matrix = lil_matrix(
        (len(peptide_index_table), len(peptide_index_table)), 
        dtype='int8'
    )
    for edge in tqdm(peptide_edges):
        metrics = edge.split()
        if len (metrics) == 7:
            peptide_one = metrics[-2]
            peptide_two = metrics[-1]
            if len(peptide_one) == 9 and len(peptide_two) == 9:
                distance_matrix[peptide_index_table.get(peptide_one), 
                                peptide_index_table.get(peptide_two)] = int(metrics[4])
                distance_matrix[peptide_index_table.get(peptide_two), 
                                peptide_index_table.get(peptide_one)] = int(metrics[4])
    
    return distance_matrix
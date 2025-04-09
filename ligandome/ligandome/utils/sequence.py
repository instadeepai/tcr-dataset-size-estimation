"""Module for working with peptide sequence data."""
from __future__ import annotations

import pandas as pd
from tqdm import tqdm
from itertools import combinations, product, chain

def extract_kmers(k: int, sequence: str, id: str) -> dict[str, str]:
    """Extract all k-mers using a sliding window method.

    Args:
        k (int): Length of k-mers to extract.
        sequence (str): Sequence to extract k-mers for.
        id (str): ID of sequence being split.

    Returns:
        dict[str, str]: Dictionary of k-mers and source ids.
    """    
    kmers = {}
    # Calculate how many kmers of length k there are
    num_kmers = len(sequence) - k + 1
    # Loop over the kmer start positions
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = sequence[i:i+k]
        kmers[kmer] = id
    # Return the final counts
    return kmers

def get_kmer_union(fasta_dict: dict[str, str], k: int) -> dict[str, str]:
    """Get a union of all k-mers for a batch of sequences.

    Args:
        fasta_dict (dict[str, str]): Dictionary of sequences, id.
        k (int): Length of k-mers to extract.

    Returns:
        dict[str, str]: Union of unique k-mers with their source sequence ids.
    """
    kmer_dict = {}
    for id, sequence in tqdm(fasta_dict.items()):
        
        kmers = extract_kmers(k, sequence, id)
        for kmer, id in kmers.items():
            kmer_dict[kmer] = id
    return kmer_dict

def mutant_generator(peptide: str, edit_distance: int, alphabet: list) -> str:
    """Generate all peptides within given edit distance of input peptide.

    Args:
        peptide (str): Peptide to generate edit distance combinations from.
        edit_distance (int): Distance of generated peptides from original.
        alphabet (list): Restricted alphabet to use.
        
    Yields:
        Iterator[str]: edit distance n peptides.
    """    
    for positions in combinations(range(len(peptide)), edit_distance):
        for replacements in product(range(len(alphabet) - 1), repeat=edit_distance):
            cousin = list(peptide)
            p = 0
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin), p

def synthesise_all_mutants(peptide: str, edit_distance: int, alphabet: list[str]) -> list[str]:
    """Generate mutants whose edit distance from peptide equal to edit_distance.
    
    Args:
        peptide (str): Peptide to generate edit distance combinations from.
        edit_distance (int): Distance of generated peptides from original.
        alphabet (list): Restricted alphabet to use.
    
    Returns:
        list[str]: List of all mutants within given edit distance.
    """
    return list(chain.from_iterable(mutant_generator(peptide, i, alphabet)
                               for i in range(edit_distance + 1)))

def get_non_overlapping_mutant_subset(peptides: list[str], 
                                      edit_distance: int, 
                                      alphabet: list[str], 
                                      existing_subset: list[str]=None) -> dict[str]:
    """Get subset of mutants at given edit distance which do not appear in the original set.
    
    For example, we can supply the entire human proteome to this function as a list of 9mers, 
    and it will return all edit distance N peptides which do not appear at other locations in
    the proteome. 

    Args:
        peptides (list[str]): Peptides to mutate.
        edit_distance (int): Edit distance threshold.
        alphabet (list[str]): Mutation alphabet options. 
        existing_subset (list[str]): Set that is checked for overlaps. If none, defaults to input.

    Returns:
        list[str]: Unique mutants not present in initial set.
    """
    input_peptides = set(peptides)
    if existing_subset is None:
        existing_subset = peptides
    existing_subset = set(existing_subset)
    mutanome = {}

    for peptide in tqdm(input_peptides):
        peptide_mutanome = synthesise_all_mutants(peptide, edit_distance, alphabet)
        for mutant_peptide, _ in peptide_mutanome:
            if mutant_peptide not in existing_subset:
                mutanome[mutant_peptide] = peptide
    
    return mutanome

def sliding_window_get_all_kmers(peptides: list[str], ids: list[str], k: int) -> pd.DataFrame:
    """Break down all peptides into 9-mers using sliding window.

    Args:
        peptides (list[str]): List of unique peptides of varying lengths.
        ids (list[str]): List of ids or sources for each peptide in list.
        k (int): Length of k-mers to extract.

    Returns:
        pd.DataFrame: Cleaned enumerated dataframe.
    """
    # for non 9-mers, break down into 9-mers
    union_neoantigens = {}

    for peptide, mutation_id in tqdm(zip(peptides, ids), total=len(peptides)):
        if len(peptide) < k:
            pass
        elif len(peptide) == k:
            if union_neoantigens.get(peptide) is None:
                    union_neoantigens[peptide] = [mutation_id]
            else:
                union_neoantigens[peptide].append(mutation_id)
        else:
            ninemer_permutations = extract_kmers(k, peptide, mutation_id)
            for kmer, id in ninemer_permutations.items():
                if union_neoantigens.get(kmer) is None:
                    union_neoantigens[kmer] = [id]
                else:
                    union_neoantigens[kmer].append(id)
    
    return pd.DataFrame({'peptide':list(union_neoantigens.keys()),
                         'peptide_source': list(union_neoantigens.values())})
    
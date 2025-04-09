"""Module for NetMHCPan utuls and python wrappers."""
from __future__ import annotations

import os
import math
import pandas as pd
from pathlib import Path
from subprocess import Popen, PIPE
from joblib import Parallel, delayed
from ligandome.utils.constants import DATABASE_EXPORTS_PATH
from ligandome.utils.cleaning import divide_list_into_chunks


def run_netmhcpan(
    prefix: str, allele: str, peptide_list: list[str], job_idx: int, temp_dir: Path
) -> dict[str, float]:
    """Run netmhcpan predictions on a list of peptides.

    Args:
        prefix(str): Prefix for filenames.
        allele (str): Which allele to run predictions against.
        peptide_list (list[str]): List of peptides to run predictions for.
        job_idx (int): Index of job.
        temp_dir (Path): Temporary directory to use for files.

    Returns:
        list[float]: List of NetMHCPan binding affinity predictions.
    """
    with open(f"{str(temp_dir)}/{prefix}_netmhcpan_{allele}_{job_idx}.txt", "a+") as pep_file:
        for peptide in peptide_list:
            pep_file.write(f"{peptide}\n")

    command = f'{str(DATABASE_EXPORTS_PATH.parent)}/third_party_tools/netMHCpan-4.1/netMHCpan -p \
                {str(temp_dir)}/{prefix}_netmhcpan_{allele}_{job_idx}.txt \
                -a "HLA-{allele[:3]}:{allele[3:]}" -BA | grep "1 HLA"'
    with Popen(command, stdout=PIPE, shell=True) as netmhcpan_run:
        results = netmhcpan_run.communicate()[0].decode("utf-8").split("\n")
        float_score_results: list[float] = [
            float(r.split()[12]) for r in results if len(r.split()) > 13
        ]
    os.remove(f"{str(temp_dir)}/{prefix}_netmhcpan_{allele}_{job_idx}.txt")
    assert len(float_score_results) == len(peptide_list), f'NetMHCPan scores have length {len(float_score_results)} whilst {len(peptide_list)} peptides were supplied.'
    return float_score_results

def get_netmhcpan_hits(df: pd.DataFrame, peptide_col: str, allele: str, el_rank_threshold: float, threads: int, job_tag: str, temp_dir: Path) -> pd.DataFrame:
    """Get subset of dataframe which are predicted binders by NetMHCPan4.1b.

    Args:
        df (pd.DataFrame): Dataframe of peptides.
        peptide_col (str): Name of column with peptide sequences.
        allele (str): Allele to predict against. 
        el_rank_threshold (float): Threshold for keeping peptides as binders.
        threads (int): Threads to parallelise scoring across.
        job_tag (str): Tag for temporary files.
        temp_dir (Path): Temporary directory to use for files.

    Returns:
        pd.DataFrame: Input dataframe subsetted to NetMHCPan predicted binders.
    """    
    peptide_list_chunks = divide_list_into_chunks(df[peptide_col], chunksize=math.ceil(len(df[peptide_col])/threads))
    nested_scores = Parallel(n_jobs=threads)(delayed(run_netmhcpan)(job_tag, allele, chunk, i, temp_dir) for i, chunk in enumerate(peptide_list_chunks))
    flattened_scores = [val for sublist in nested_scores for val in sublist]
    df['netmhcpan_rank_el'] = flattened_scores
    return df.loc[df['netmhcpan_rank_el'] < el_rank_threshold].reset_index(drop=True)

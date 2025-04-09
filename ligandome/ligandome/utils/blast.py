"""Module for blastp utils and python wrappers."""
from __future__ import annotations

import numpy as np
import pandas as pd
from io import StringIO
from pathlib import Path
from subprocess import Popen, PIPE
from ligandome.utils.constants import DATABASE_EXPORTS_PATH
from ligandome.utils.fasta import peptide_list_to_fasta_file

def make_blast_protein_db(fasta_file: Path, db_outdir: Path, blast_path: Path=Path(f'{str(DATABASE_EXPORTS_PATH.parent)}/third_party_tools/ncbi-blast-2.14.1+/bin/makeblastdb')) -> None:
    """Create blast protein database from a fasta file.

    Args:
        fasta_file (Path): Fasta file of protein sequence.
        db_outdir (Path): Where to save the output protein blastdb
        blast_path (Path, optional): Where the blast executable is located. Defaults to Path(f'{str(DATABASE_EXPORTS_PATH.parent)}/third_party_tools/ncbi-blast-2.14.1+/bin/makeblastdb').
    """    
    command = f'{str(blast_path)} -in {str(fasta_file)} -dbtype prot -out {str(db_outdir)}'
    Popen(command, stdout=PIPE, shell=True).communicate()

    
def blast_fasta_against_db(fasta_file: Path, 
                           db_path: Path=Path(f'{str(DATABASE_EXPORTS_PATH)}/blast_dbs/UniProtCanonicalProteome'), 
                           blast_path: Path=Path(f'{str(DATABASE_EXPORTS_PATH.parent)}/third_party_tools/ncbi-blast-2.14.1+/bin/blastp')) -> pd.DataFrame:
    """Run BLASTP of a fasta file of sequences against a given blast database.

    Args:
        fasta_file (Path): Filepath of fasta ro run blast for.
        db_path (Path, optional): Database to blast against. Defaults to Path('../Data/Airlock/BlastDBs/UniProtCanonicalProteome.fasta').
        blast_path (Path, optional): Path to blast binary to use. Defaults to Path('../Software/ncbi-blast-2.14.1+/bin/blastp').

    Returns:
        pd.DataFrame: Dataframe of hits.
    """
    command = f'{str(blast_path)} -query {str(fasta_file)} -db {str(db_path)} -outfmt 6'
    blast_run = Popen(command, stdout=PIPE, shell=True)
    results = StringIO(blast_run.communicate()[0].decode('utf-8'))
    return pd.read_csv(results, sep="\t", names=['qseqid',
                                                 'sseqid',
                                                 'pident',
                                                 'length',
                                                 'mismatch',
                                                 'gapopen',
                                                 'qstart',
                                                 'qend',
                                                 'sstart',
                                                 'send',
                                                 'evalue',
                                                 'bitscore'])

def human_proteome_blast_filter(df: pd.DataFrame, tmp_dir: Path, alignment_lengths: list[int], identity_threshold: int, job_tag: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Filter dataframe of peptides into matches and non-matches to human proteome.

    Args:
        df (pd.DataFrame): Dataframe of peptides to blast.
        tmp_dir (Path): Temporary directory for working files.
        alignment_lengths (list[int]): Lengths of alignments accepted as matches.
        identity_threshold (int): Identity percentage threshold accepted as matches.
        job_tag (str): Tag for temporary job files.

    Returns:
        tuple[pd.DataFrame, pd.DataFrame]: Dataframe of matched peptides, dataframe of unmatched peptides.
    """    
    peptide_list_to_fasta_file(df['peptide'], tmp_dir / f'{job_tag}Peptides.fasta')
    blast_results = blast_fasta_against_db(tmp_dir / f'{job_tag}Peptides.fasta')
    blast_results['wt_map_bool'] = np.where(((blast_results['pident'] > identity_threshold) & (blast_results['length'].isin(alignment_lengths))), 1, 0)
    blast_results = blast_results.sort_values(by=['qseqid','wt_map_bool','pident']).drop_duplicates(keep='first')
    df['uniprot_wt'] = df.index.map(dict(zip(blast_results['qseqid'], blast_results['sseqid'])))
    df['wt_alignment_start'] = df.index.map(dict(zip(blast_results['qseqid'], blast_results['sstart'])))
    df['wt_alignment_end'] = df.index.map(dict(zip(blast_results['qseqid'], blast_results['send'])))
    df['wt_map_bool'] = df.index.map(dict(zip(blast_results['qseqid'], blast_results['wt_map_bool']))).fillna(0)
    unmapped_to_proteome = df.loc[df['wt_map_bool'] == 0]
    mapped_to_proteome = df.loc[df['wt_map_bool'] == 1]
    return mapped_to_proteome, unmapped_to_proteome

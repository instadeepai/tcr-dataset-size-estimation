"""Script for fetching the non-self ligandome for a given allele from public databases."""
from __future__ import annotations

import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from ligandome.utils.fasta import read_fasta_file
from ligandome.utils.cleaning import aggregate_and_clean_peptide_df
from ligandome.utils.constants import (
    DATABASE_EXPORTS_PATH, 
    AMINO_ACIDS, 
    KMER_LENGTH
)
from ligandome.utils.netmhcpan import get_netmhcpan_hits
from ligandome.utils.blast import human_proteome_blast_filter, make_blast_protein_db
from ligandome.utils.sequence import (
    get_kmer_union,
    get_non_overlapping_mutant_subset,
)
from ligandome.utils.databases import query_all_databases
from ligandome.utils.graphs import calculate_dominating_set_members

HLA_LIGAND_ATLAS_DATA = DATABASE_EXPORTS_PATH / 'HLA_ligand_atlas_aggregated.csv'
IEDB_HEALTHY_DATA = DATABASE_EXPORTS_PATH / 'IEDB_healthy_data.csv'
IEDB_TUMOUR_DATA = DATABASE_EXPORTS_PATH / 'IEDB_tumour_data.csv'
IEDB_VIRAL_DATA = DATABASE_EXPORTS_PATH / 'IEDB_viral_data.csv'
IEDB_MHC_MAPPING_DATA = DATABASE_EXPORTS_PATH / 'IEDB_ligand_full.csv'
VDJDB_VIRAL_DATA = DATABASE_EXPORTS_PATH / 'VDJDB_viral.csv'
NETMHCPAN_DATA = DATABASE_EXPORTS_PATH / 'NetMHCPan_monoallelic_training_data.csv'
UNIPROT_HUMAN_PROTEOME_FASTA = DATABASE_EXPORTS_PATH / 'UniProtCanonicalProteome.fasta'

def parse_arguments() -> argparse.ArgumentParser:
    """Parse command line arguments.

    Returns:
        argparse.ArgumentParser: Parsed user arguments.
    """    

    parser = argparse.ArgumentParser()

    parser.add_argument('--allele', '-hla', required=True, type=str, help='Which HLA to calculate ligandome for in simplified format (e.g. A0201).')
    parser.add_argument('--output','-o', required=True, type=str, help='Folder to store the output results; will be created if it does not already exist.')
    parser.add_argument('--dominating_set','-d', required=False, action='store_true', default=False, help='Whether to calculate the dominating sets at edit distances 1, 2 and 3.')
    parser.add_argument('--edit_distance_thresholds', '-e', required=False, type=str, default='1,2,3', help='Threshold(s) to calculate dominating sets for.')
    parser.add_argument('--threads', '-t', required=False, type=int, default=1, help='Number of threads to parallelise workflow over.')

    args = parser.parse_args()
    args.edit_distance_thresholds = [int(threshold) for threshold in args.edit_distance_thresholds.split(',')]

    return args

def setup_temp_dirs(args: argparse.ArgumentParser) -> tuple[argparse.ArgumentParser, Path]:
    """Setup timestamped working dirs for ligandome runs.

    Args:
        args (argparse.ArgumentParser): Parsed user arguments.

    Returns:
        tuple[argparse.ArgumentParser, Path]: Parsed arguments and temporary directory path.
    """    
    args.output = Path(args.output) / args.allele
    args.output.mkdir(exist_ok=True, parents=True)
    tag = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
    tmp_dir = args.output / tag
    tmp_dir.mkdir(exist_ok=True)
    (tmp_dir / 'for_netmhcpan_inference').mkdir(exist_ok=True)
    (tmp_dir / 'for_sachica_inference').mkdir(exist_ok=True)
    return args, tmp_dir

def main():
    """
    Main workflow call
    """
    args = parse_arguments()
    args, tmp_dir = setup_temp_dirs(args)

    # Make sure our blastdb exists
    if not (DATABASE_EXPORTS_PATH / 'blast_dbs' / 'UniProtCanonicalProteome').exists():
        make_blast_protein_db(UNIPROT_HUMAN_PROTEOME_FASTA,
                              DATABASE_EXPORTS_PATH / 'blast_dbs' / 'UniProtCanonicalProteome')

    # Query all our databases.
    experimental_full = query_all_databases(hla_ligand_atlas_data=HLA_LIGAND_ATLAS_DATA,
                                            iedb_healthy_data=IEDB_HEALTHY_DATA,
                                            iedb_tumour_data=IEDB_TUMOUR_DATA,
                                            iedb_viral_data=IEDB_VIRAL_DATA,
                                            iedb_mhc_mapping_data=IEDB_MHC_MAPPING_DATA,
                                            vdjdb_viral_data=VDJDB_VIRAL_DATA,
                                            netmhcpan_data=NETMHCPAN_DATA,
                                            allele=args.allele)

    # Split proteome into 9-mers for comparison
    human_proteome = read_fasta_file(UNIPROT_HUMAN_PROTEOME_FASTA)
    human_proteome_ninemers = get_kmer_union(human_proteome, k=KMER_LENGTH)
    
    # Separate self and non-self via a blast search
    endogenous_experimental = experimental_full.loc[experimental_full['peptide_origin'] == 'Healthy Tissues']
    exogenous_experimental = experimental_full.loc[experimental_full['peptide_origin'] != 'Healthy Tissues']
    endogenous_experimental = aggregate_and_clean_peptide_df(endogenous_experimental, agg_col='peptide')
    endogenous_experimental, unmapped_experimental = human_proteome_blast_filter(df=endogenous_experimental, 
                                                                                 tmp_dir=tmp_dir,
                                                                                 alignment_lengths=[8,9,10],
                                                                                 identity_threshold=88,
                                                                                 job_tag='healthy')
    
    # Create a mutanome for verified self peptides avoiding re-creating self peptides
    experimental_mutanome = get_non_overlapping_mutant_subset(endogenous_experimental['peptide'], 
                                                              edit_distance=1, 
                                                              alphabet=AMINO_ACIDS, 
                                                              existing_subset=list(human_proteome_ninemers.keys()))
    experimental_mutanome = pd.DataFrame({'peptide':experimental_mutanome.keys(),
                                        'mhc_allele':args.allele,
                                        'source_peptide':experimental_mutanome.values(),
                                        'peptide_origin':'Synthetic Healthy Mutanome'})
    
    # Create a mutanome for database tumour origin peptides avoiding re-creating self peptides
    tumour_origin_data = exogenous_experimental.loc[exogenous_experimental['peptide_origin'] == 'Tumour Sample']
    tumour_mutanome = get_non_overlapping_mutant_subset(tumour_origin_data['peptide'], 
                                                        edit_distance=1, 
                                                        alphabet=AMINO_ACIDS, 
                                                        existing_subset=list(human_proteome_ninemers.keys()))
    tumour_mutanome = pd.DataFrame({'peptide':tumour_mutanome.keys(),
                                    'mhc_allele':args.allele,
                                    'source_peptide':tumour_mutanome.values(),
                                    'peptide_origin':'Synthetic Tumour Mutanome'})
    
    # Keep only predicted NetMHCPan4.1 hits from our created mutanomes
    tumour_mutanome_hits = get_netmhcpan_hits(df=tumour_mutanome,
                                              peptide_col='peptide',
                                              allele=args.allele,
                                              el_rank_threshold=0.5,
                                              threads=args.threads,
                                              job_tag='tumour',
                                              temp_dir=tmp_dir / 'for_netmhcpan_inference')
    
    experimental_mutanome_hits = get_netmhcpan_hits(df=experimental_mutanome,
                                                    peptide_col='peptide',
                                                    allele=args.allele,
                                                    el_rank_threshold=0.5,
                                                    threads=args.threads,
                                                    job_tag='experimental',
                                                    temp_dir=tmp_dir / 'for_netmhcpan_inference')
    
    # Filter experimental non-self peptides which match human proteome exactly
    exogenous_experimental = aggregate_and_clean_peptide_df(exogenous_experimental, agg_col='peptide')
    proteome_matches, exogenous_experimental = human_proteome_blast_filter(df=exogenous_experimental,
                                                                           tmp_dir=tmp_dir,
                                                                           alignment_lengths=[9],
                                                                           identity_threshold=99,
                                                                           job_tag='exogenous')
    
    # Concatenate predicted mutanome binders and experimental non-self peptides
    final_exp_ligandome = pd.concat(
        [experimental_mutanome_hits, tumour_mutanome_hits, exogenous_experimental]
    )

    # Add sources for experimental tumour samples distance from wt peptides
    final_exp_ligandome['source_peptide'] = final_exp_ligandome['source_peptide'].fillna(final_exp_ligandome['peptide'])

    # Aggregate to remove duplicate peptides
    final_exp_ligandome = (
        final_exp_ligandome.groupby(["peptide"])
        .agg(
            {
                'mhc_allele': lambda x: x.iloc[0],
                'source_peptide': lambda x: x.tolist(),
                'peptide_origin': lambda x: x.tolist(),
                'netmhcpan_rank_el': lambda x: x.iloc[0],
                'data_source':lambda x: x.tolist(),
                'uniprot_wt':lambda x: x.tolist(),
                'wt_alignment_start':lambda x: x.tolist(),
                'wt_alignment_end': lambda x: x.tolist()
            }
        )
        .reset_index()
    )

    # Create data sources dict
    data_sources = experimental_full.groupby('peptide').agg({'data_source':
                                                             lambda x: list(set(x.tolist()))}
                                                             ).reset_index()
    data_sources = {peptide: source for peptide, source in zip(data_sources['peptide'],
                                                               data_sources['data_source'])}

    # Cleanup columns
    final_exp_ligandome['mhc_allele'] = args.allele

    final_exp_ligandome['data_source'] = final_exp_ligandome['source_peptide'].apply(lambda x: [data_sources.get(y) for y in x])

    final_exp_ligandome['uniprot_wt'] = final_exp_ligandome['source_peptide'].apply(lambda x: [dict(zip(endogenous_experimental['peptide'],
                                                                                                     endogenous_experimental['uniprot_wt'])).get(y) for y in x])
    final_exp_ligandome['wt_alignment_start'] = final_exp_ligandome['source_peptide'].apply(lambda x: [dict(zip(endogenous_experimental['peptide'],
                                                                                                     endogenous_experimental['wt_alignment_start'])).get(y) for y in x])
    final_exp_ligandome['wt_alignment_end'] = final_exp_ligandome['source_peptide'].apply(lambda x: [dict(zip(endogenous_experimental['peptide'],
                                                                                                     endogenous_experimental['wt_alignment_end'])).get(y) for y in x])
    
    final_exp_ligandome.to_csv(args.output / f"{args.allele}ExperimentalLigandome.csv", index=False)
    
    # Calculate dominating set members
    if args.dominating_set:
        for threshold in args.edit_distance_thresholds:
            final_exp_ligandome[f'dominating_set_edit_distance_{threshold}'] = calculate_dominating_set_members(
                list(final_exp_ligandome['peptide']),
                edit_distance_threshold=threshold,
                tmpfile_dir=tmp_dir / 'for_sachica_inference'
            )
            final_exp_ligandome[f'dominating_set_edit_distance_{threshold}'] = final_exp_ligandome[f'dominating_set_edit_distance_{threshold}'].astype(int)

    final_exp_ligandome.to_csv(args.output / f"{args.allele}ExperimentalLigandome.csv", index=False)

if __name__ == "__main__":
    main()
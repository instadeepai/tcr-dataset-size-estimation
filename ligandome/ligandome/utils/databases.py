"""Module with database-specific functions for standardising peptide data."""
from __future__ import annotations

import pandas as pd
from pathlib import Path
from ligandome.utils.cleaning import clean_mhc_annotations

def query_all_databases(hla_ligand_atlas_data: Path,
                        iedb_mhc_mapping_data: Path,
                        iedb_healthy_data: Path,
                        iedb_viral_data: Path,
                        iedb_tumour_data: Path,
                        vdjdb_viral_data: Path,
                        netmhcpan_data: Path,
                        allele: str) -> pd.DataFrame:
    """Wrapper function for querying all sources for allele specific data.

    Args:
        hla_ligand_atlas_data (Path): Path to HLA Ligand Atlas export.
        iedb_mhc_mapping_data (Path): Path to IEDB allele mapping export.
        iedb_healthy_data (Path): Path to IEDB healthy data export.
        iedb_viral_data (Path): Path to IEDB viral data export.
        iedb_tumour_data (Path): Path to IEDB tumour data export.
        vdjdb_viral_data (Path): Path to VDJDB viral data export.
        netmhcpan_data (Path): Path to NetMHCPan training data export.
        allele (str): Allele name to query databases for.

    Returns:
        pd.DataFrame: Dataframe of allele specific peptides. 
    """    
    # clean and aggregate all data
    hla_data = parse_hla_ligand_data(hla_ligand_atlas_data, allele=allele)
    iedb_epitope_mapping = prepare_iedb_epitope_dict(iedb_mhc_mapping_data)
    iedb_healthy = parse_iedb_exported_data(iedb_healthy_data, 'Healthy Tissues', iedb_epitope_mapping, allele=allele)
    iedb_viral = parse_iedb_exported_data(iedb_viral_data, 'Viral Proteome', iedb_epitope_mapping, allele=allele)
    iedb_tumour = parse_iedb_exported_data(iedb_tumour_data, 'Tumour Sample', iedb_epitope_mapping, allele=allele)
    vdjdb_viral = parse_vdjdb_annotated_data(vdjdb_viral_data, 'Viral Proteome', allele=allele)
    net_mhcpan_data = parse_netmhcpan_training_data(netmhcpan_data, allele=allele)
    experimental_full = pd.concat([hla_data,
                                iedb_healthy,
                                iedb_viral,
                                iedb_tumour,
                                vdjdb_viral,
                                net_mhcpan_data])
    experimental_full['mhc_allele'] = allele

    return experimental_full

def parse_hla_ligand_data(datasheet: Path, allele: str) -> pd.DataFrame:
    """parse and standardise HLA Ligand Atlas raw data.

    Args:
        datasheet (Path): HLA Ligand atlas raw data.
        allele (str): Allele to select data for in simplified format
            (e.g. A0201, not A*02:01)

    Returns:
        pd.DataFrame: Dataframe of parsed HLA data.
    """
    df = pd.read_csv(datasheet)
    
    # take only MHC I
    df = df.loc[df['hla_class'].isin(['HLA-I','HLA-I+II'])]
    
    # subset to 9mers
    df['peptide length'] = df['peptide_sequence'].apply(lambda x: len(x))
    df = df.loc[df['peptide length'] == 9]
    
    # parse mhc alleles
    df['mhc_allele'] = df['allele']
    
    # standardise columns
    df['peptide'] = df['peptide_sequence']
    df['data_source'] = 'HLA Ligand Atlas'
    df['peptide_origin'] = 'Healthy Tissues'
    df = df[['peptide','mhc_allele','data_source','peptide_origin']]
    
    return df.loc[df['mhc_allele'].str.contains(allele)]

def parse_vdjdb_annotated_data(datasheet: Path, origin_label: str, allele: str) -> pd.DataFrame:
    """parse and standardise VDJDB raw data.

    Args:
        datasheet (Path): Path to CDJDB exported data.
        origin_label (str): Label for peptide_origin column.
        allele (str): Allele to select data for in simplified format
            (e.g. A0201, not A*02:01)

    Returns:
        pd.DataFrame: Dataframe of standardised parsed data.
    """
    df = pd.read_csv(datasheet)
    df.drop_duplicates(subset='Epitope', inplace=True)
    df['peptide length'] = df['Epitope'].apply(lambda x: len(x))
    df = df.loc[df['peptide length'] == 9]
    df['peptide'] = df['Epitope']
    df['data_source'] = 'VDJDB'
    df['peptide_origin'] = origin_label
    df['mhc_allele'] = df['MHC A']
    df['mhc_allele'] = df['mhc_allele'].apply(clean_mhc_annotations)
    df = df[['peptide','mhc_allele','data_source','peptide_origin']]
    return df.loc[df['mhc_allele'] == allele]

def parse_iedb_exported_data(datasheet: Path, origin_label: str, mhc_epitope_dict: dict[str, str], allele: str) -> pd.DataFrame:
    """parse and standardise IEDB raw data.

    Args:
        datasheet (Path): IEDB raw data.
        origin_label (str): Label to assign to data.
        mhc_epitope_dict (dict[str, str]): Mapping from epitope ID to MHC alleles.
        allele (str): Allele to select data for in simplified format
            (e.g. A0201, not A*02:01)

    Returns:
        pd.DataFrame: Dataframe of parsed IEDB data.
    """
    df = pd.read_csv(datasheet)
    df['EpitopeID'] = df['Epitope ID - IEDB IRI'].apply(lambda x: str(x).split('/')[-1]).astype(int)
    df['MHC_allele'] = df['EpitopeID'].map(mhc_epitope_dict)
    df = df.loc[df['Epitope - Name'].apply(lambda x: len(x)) == 9]
    df['peptide length'] = df['Epitope - Name'].apply(lambda x: len(x))
    df['mhc_allele'] = df['MHC_allele']
    df['peptide'] = df['Epitope - Name']
    df['data_source'] = 'IEDB'
    df['peptide_origin'] = origin_label
    df = df[['peptide','mhc_allele','data_source','peptide_origin']]

    return df.loc[df['mhc_allele'].str.contains(allele)]

def prepare_iedb_epitope_dict(datasheet: Path) -> dict[str, str]:
    """Prepare MHC allele annotations for IEDB epitopes.

    Args:
        datasheet (Path): Path to IEDB full MHC mapping sheet.

    Returns:
        dict[str, str]: Mapping from epitope ID to MHC allele list.
    """    
    # add mhc screening information
    mhc_healthy_key = pd.read_csv(datasheet)
    mhc_healthy_key = mhc_healthy_key[['Epitope IRI','mhc_allele']]
    mhc_healthy_key['mhc_allele'] = mhc_healthy_key['mhc_allele'].apply(clean_mhc_annotations)
    mhc_healthy_key = mhc_healthy_key.groupby(['Epitope IRI']).agg({'mhc_allele': lambda x: str(x.tolist())})
    mhc_healthy_key = mhc_healthy_key.reset_index()

    return dict(zip(mhc_healthy_key['Epitope IRI'], mhc_healthy_key['mhc_allele']))


def parse_netmhcpan_training_data(datasheet: Path, allele: str) -> pd.DataFrame:
    """parse and standardise data from NetMHCPan.

    Args:
        datasheet (Path): Path to NetMHCPan raw data.
        allele (str): Allele to select data for in simplified format
            (e.g. A0201, not A*02:01)

    Returns:
        pd.DataFrame: parsed and standardised dataframe.
    """
    df = pd.read_csv(datasheet)
    
    df = df.loc[df['allele'].str.contains('HLA')]
    df['peptideLength'] = df['peptide'].apply(len)
    df = df.loc[df['peptideLength'] == 9]
    df = df.loc[df['Target'] > 0.5]

    # keep all our relevant MHCs
    df['mhc_allele'] = df['allele']
    df['data_source'] = 'NetMHCPan'
    df['peptide_origin'] = 'Healthy Tissues'

    df = df[['peptide','mhc_allele','data_source','peptide_origin']]

    df['mhc_allele'] = df['mhc_allele'].apply(clean_mhc_annotations)

    df = df[['peptide','mhc_allele','data_source','peptide_origin']]
    return df.loc[df['mhc_allele'] == allele]
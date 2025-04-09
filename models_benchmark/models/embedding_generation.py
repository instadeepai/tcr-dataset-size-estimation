"""Functions to generate embeddings."""
from pathlib import Path

import hydra
import numpy as np
import pandas as pd
import torch
from omegaconf import DictConfig, OmegaConf

from models.utils.embedding import (
    apply_pca,
    atchley_embedding,
    embeddings_from_tokens,
    get_token,
    get_tokenizer_model,
)
from models.utils.generic import load_config


class EmbeddingGenerator:
    """Class to generate embeddings based on a csv with sequences and a configuration.
.
    Given sequences, the run command preprocess the sequence csv file, computes embeddings
    for CDR3b and peptide based on initial configuration.
    """

    def __init__(self, data_dir: Path) -> None:
        """Initialize an EmbeddingGenerator object.

        Args:
            data_dir (Path): Directory where embeddings and datasets are saved.
        """
        self.peptide_len = 15
        self.data_dir = data_dir
        self.device = (
            torch.device("cuda:0")
            if torch.cuda.is_available() 
            else torch.device("cpu")
        )

    @staticmethod
    def _preprocess_csv(file: Path) -> pd.DataFrame:
        """Process the dataset by removing samples with antigens longer than 15AA and HLAs without sequence mapping.

        Args:
            file (Path): Path to csv containing sequences.

        Returns:
            dataset (pd.DataFrame): DataFrame with cdr3b, peptide and labels in the processed dataset.
        """
        # Preprocess TCR files
        assert file.exists()

        dataset = pd.read_csv(file, header=0)

        dataset = dataset[["cdr3b", "peptide", "binder", "fold"]]
        dataset.columns = ["cdr3b", "peptide", "binder", "fold"]

        dataset["cdr3b"] = dataset["cdr3b"].apply(lambda s: s.replace("X", ""))
        dataset["peptide"] = dataset["peptide"].apply(lambda s: s.replace("X", ""))

        dataset["cdr3b"] = dataset["cdr3b"].apply(lambda s: s.replace("-", ""))
        dataset["peptide"] = dataset["peptide"].apply(lambda s: s.replace("-", ""))

        # Remove antigens longer than 15aa
        dataset = dataset[dataset.peptide.str.len() < 16]
        dataset = dataset[dataset.peptide.str.len() < 35]

        return dataset

    def _get_embeddings(self, df: pd.DataFrame) -> pd.DataFrame:
        """Compute sequence embeddings.

        CDR3b embeddings are obtained using TCR-Bert.
        Peptide embeddings are obtained with the Atchley Factors.

        Args:
            df (pd.DataFrame): Dataframe with cdr3b and peptide sequences.

        Returns:
            embeddings (np.ndarray): Concatenated cdr3b and peptide embeddings.
        """
        model, tokenizer, transformer_dim = get_tokenizer_model()

        cdr3b_token = get_token(tokenizer, df["cdr3b"])
        cdr3b_embeddings = embeddings_from_tokens(
            cdr3b_token, transformer_dim, model, batch_size=64, device=self.device
        )
        cdr3b_embeddings = apply_pca(cdr3b_embeddings)
        
        df['peptide_padded'] = df['peptide'].apply(lambda s: s + "X" * (self.peptide_len - len(s)))
        peptide_embeddings = df["peptide_padded"].apply(atchley_embedding)
        peptide_embeddings = np.stack(peptide_embeddings.to_numpy())
        peptide_embeddings = peptide_embeddings.reshape(-1, self.peptide_len*5)
        
        embeddings = np.hstack([cdr3b_embeddings, peptide_embeddings])
        return embeddings

    def run(
        self,
        csv_path: Path,
    ) -> pd.DataFrame:
        """Run embedding generation given a csv with sequences and labels.

        Embeddings and sequences are saved in designated output directory.

        Args:
            csv_path (Path): File path to csv containing sequences and labels.
        """
        df = EmbeddingGenerator._preprocess_csv(csv_path)

        embeddings = self._get_embeddings(df)

        df = df.reset_index(drop=True)
        idx_test = df.loc[df.fold==1].index.to_list()
        idx_train = df.loc[df.fold==2].index.to_list()

        df_test = df.loc[df.fold==1]
        df_train = df.loc[df.fold==2]

        embeddings_train = embeddings[idx_train]
        embeddings_test = embeddings[idx_test]

        train_dir = self.data_dir / "train"
        test_dir = self.data_dir / "test"

        train_dir.mkdir(exist_ok=True, parents=True)
        test_dir.mkdir(exist_ok=True, parents=True)

        df_train.to_csv(train_dir / "dataset.csv", index=False)
        df_test.to_csv(test_dir / "dataset.csv", index=False)

        np.save(train_dir / "embeddings.npy", embeddings_train)
        np.save(test_dir / "embeddings.npy", embeddings_test)


@hydra.main(config_path="config/", config_name="config")
def main(args: DictConfig) -> None:
    """Given csv files with sequences, compute embeddings of the dataset."""
    OmegaConf.set_struct(args, False)
    conf_data, _ = load_config(args)
    embedding_generator = EmbeddingGenerator(Path(conf_data.data_dir))
    embedding_generator.run(Path(conf_data.csv_path))


if __name__ == "__main__":
    main()

"""Dataset classes."""
from typing import Any, Dict, List, Tuple

import numpy as np
import torch
from torch.utils.data import Dataset


class ClassifierDataset:
    """Dataset class for classifier models."""

    def __init__(
        self,
        embeddings: np.ndarray,
        y_true: np.ndarray,
        cdr3b: List[str],
        peptide: List[str],
    ) -> None:
        """Initialize a ClassifierDataset object.

        Args:
            embeddings (np.ndarray): Combined CDR3b-peptide embeddings.
            y_true (np.ndarray): Binder label.
            cdr3b (List[str]): List of CDR3b sequences.
            peptide (List[str]): List of peptide sequences.
        """
        self.data = embeddings
        self.sequences = [(cdr3b_, peptide_) for (cdr3b_, peptide_) in zip(cdr3b, peptide)]
        self.cdr3b = cdr3b
        self.peptide = peptide
        self.n_tot = self.data.shape[0]
        self.y = y_true

    def __getitem__(self, item: Any) -> Tuple[np.ndarray, float]:
        """Load the data corresponding to the given index.

        The returned data contains the embeddings and the binding label of the corresponding item.

        Args:
            item (Any): Sample index, value will be 0 to self.len()-1.

        Returns:
            (Tuple[np.ndarray, float]): Loaded data (embeddings and binding label).
        """
        return self.data[item], self.y[item]

    def __len__(self) -> int:
        """Returns the length of the dataset."""
        return self.n_tot

    def get_sequence(self, item: Any) -> Tuple[str, str]:
        """Load sequences of a given index in the same way as getitem.

        Args:
            item: Sample index, value will be 0 to self.len()-1.

        Returns:
            (Tuple[str, str]): Sequences of the index (CDR3b, peptide).
        """
        return self.sequences[item]


class TokenDataset(Dataset):
    """Dataset class to compute transformer embeddings."""

    def __init__(self, data: torch.Tensor) -> None:
        """Init.

        Args:
            data (torch.Tensor): tokenized sequences.
        """
        self.data = data

    def __getitem__(self, item: Any) -> Dict:
        """Load the data corresponding to the given index.

        The returned data object contains input ids, attention and masks vectors.

        Args:
            item (Any): sample index, value will be 0 to self.len()-1.

        Returns:
            (Dict): Loaded data.
        """
        return {key: self.data[key][item] for key in self.data}

    def __len__(self) -> int:
        """Returns the length of the dataset."""
        return self.data["input_ids"].shape[0]

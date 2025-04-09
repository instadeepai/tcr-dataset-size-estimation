"""Utils to compute embeddings."""
from typing import List, Tuple

import numpy as np
import torch
from sklearn.decomposition import PCA
from torch.utils.data import DataLoader
from tqdm import tqdm
from transformers import AutoModelForMaskedLM, AutoTokenizer

from models.dataset import TokenDataset
from models.utils.constants import ATCHLEY_FACTORS, BERT_DIM

cpu_device = torch.device("cpu")


def atchley_embedding(seq: str) -> np.ndarray:
    """Convert sequences to Atchley factors embeddings.

    Args:
        seq (str): Sequence to embed.

    Returns:
        embeddings (np.ndarray): Embedded sequence.
    """
    embeddings = np.array([ATCHLEY_FACTORS[aa] for aa in seq])

    return np.array(embeddings)


def get_tokenizer_model() -> Tuple[AutoModelForMaskedLM, AutoTokenizer, int]:
    """Return pre-trained TCRBert language model, tokenizer and embedding dimension.

    Returns:
        model, tokenizer, transformer_dim: Pre-trained TCRBert model, tokenizer and transformer_dim
    """
    tokenizer = AutoTokenizer.from_pretrained("wukevin/tcr-bert")
    model = AutoModelForMaskedLM.from_pretrained("wukevin/tcr-bert")
    model.bert.encoder.layer = model.bert.encoder.layer[:8]
    model = model.bert
    transformer_dim = BERT_DIM

    return model, tokenizer, transformer_dim


def get_token(
    tokenizer: AutoTokenizer,
    seq_list: List[str],
    max_len: int = 35,
) -> torch.Tensor:
    """Convert sequences into tokens.

    Args:
        tokenizer (AutoTokenizer): Tokenizer object corresponding to the HF model.
        seq_list (List[str]): List of sequences to be converted.
        max_len (int): Maximum sequence length (for padding).

    Returns:
        tokens (torch.Tensor): Array containing tokenized sequences.
    """
    # For TCR bert, AAs are separated by spaces.
    seq_list = [" ".join(seq) for seq in seq_list]
    tokens = tokenizer(seq_list, padding="max_length", max_length=max_len, return_tensors="pt")
    return tokens


def embeddings_from_tokens(
    tokens: torch.Tensor,
    emb_dim: int,
    transformer: AutoModelForMaskedLM,
    batch_size: int = 64,
    device: torch.device = cpu_device,
) -> np.ndarray:
    """Compute transformer embeddings using tokenized sequences.

    Args:
        tokens (torch.Tensor): Array with tokenized sequences.
        emb_dim (int): Output embedding dimension of the model.
        transformer (AutoModelForMaskedLM): Pre-trained model.
        batch_size (int): Batch size for compute transformer embeddings.
        device (torch.device): Device on which to compute transformer embeddings.

    Returns:
        embeddings (np.ndarray): Output embeddings of the sequences.
    """
    dataloader = DataLoader(TokenDataset(tokens), batch_size=batch_size, shuffle=False)
    embeddings = torch.zeros(tokens["input_ids"].shape[0], emb_dim, device=device)

    transformer = transformer.to(device)
    transformer.eval()
    sample_index = 0
    for batch in tqdm(dataloader):
        batch = {k: batch[k].to(device) for k in batch}
        batch_emb = transformer(**batch).last_hidden_state.detach()

        for j, _ in enumerate(batch["input_ids"]):
            seq_len = torch.where(batch["attention_mask"][j] == 1)[0][-1].item() - 1
            seq_emb = batch_emb[j, 1 : seq_len + 1, :].mean(0)
            embeddings[sample_index] = seq_emb
            sample_index += 1

    return embeddings.cpu().numpy()


def apply_pca(cdr3b_embeddings: np.ndarray, n_components: int = 40) -> np.ndarray:
    """Reduce CDR3b embeddings dimension using PCA.

    Args:
        cdr3b_embeddings (np.ndarray): Array containing CDR3b embeddings.
        n_components (int): Number of components of the PCA. Default 40.

    Returns:
        cdr3b_embeddings_reduced (np.ndarray): CDR3b embeddings with reduced dimension.
    """
    cdr3b_embeddings = cdr3b_embeddings
    pca = PCA(n_components=n_components)
    cdr3b_embeddings_reduced = pca.fit_transform(cdr3b_embeddings)

    return cdr3b_embeddings_reduced

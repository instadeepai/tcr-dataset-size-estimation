"""Generic utils functions."""
import inspect
from pathlib import Path
from typing import List, Tuple, Union

import numpy as np
import pandas as pd

from omegaconf import DictConfig
from sklearn import ensemble
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier

from models.dataset import ClassifierDataset

MODEL_TYPE_TO_CLASS = {
    "RandomForestClassifier": RandomForestClassifier,
    "GradientBoostingClassifier": GradientBoostingClassifier,
}

def load_config(
    args: DictConfig,
) -> Tuple[DictConfig, DictConfig]:
    """Return config values.

    Returns:
        args.data (DictConfig): Dictionary with embedding functions' config.
        args.model (DictConfig): Dictionary with train/score functions' config.
    """
    return args.data, args.model


def create_classifier_dataset(
    data_dir: Path,
) -> ClassifierDataset:
    """Load pre-computed embeddings and create a dataset for the classification models.

    Args:
        data_dir (Path): Path to the folder containing the embeddings.

    Returns:
        dataset (ClassifierDataset): Dataset for the classification models.
    """
    embeddings = np.load(data_dir / "embeddings.npy")
    dataset = pd.read_csv(data_dir / "dataset.csv")
    y_true = dataset['binder'].to_numpy()
    cdr3b = dataset['cdr3b'].to_list()
    peptide = dataset['peptide'].to_list()

    return  ClassifierDataset(embeddings, y_true, cdr3b, peptide)


def get_model_params(model_type: str) -> List[str]:
    """Return list with arguments passed to a given model class.

    Args:
        model_type (str): Class of the model to create.

    Returns:
        params (List[str]): List of parameters to initialize the model.
    """
    try:
        model_class = MODEL_TYPE_TO_CLASS[model_type]
        params = list(inspect.signature(model_class.__init__).parameters)
        params.remove("self")

    except KeyError as model_class_inexistent:
        raise KeyError(
            "Model type must be in MODEL_TYPE_TO_CLASS."
        ) from model_class_inexistent

    return params


def create_model(config: DictConfig) -> Union[RandomForestClassifier, GradientBoostingClassifier]:
    """Create model object based on model type and associated properties.

    Args:
        config (DictConfig): Model configuration.

    Returns:
        model (Union[RandomForestClassifier, GradientBoostingClassifier]): Created model.
    """
    model_classes = dict(inspect.getmembers(ensemble, predicate=inspect.isclass))
    model_cls = model_classes[config.model_type]
    init_params_names = get_model_params(config.model_type)
    init_params = {
        init_param_name: getattr(config, init_param_name)
        for init_param_name in init_params_names
        if init_param_name != "self" and init_param_name in config
    }
    return model_cls(**init_params)


def save_labels_predictions(
    output_dir: Path,
    y_true: np.ndarray,
    preds: np.ndarray,
    sequences: List[Tuple[str, str]],
    model_type: str,
    fold_index: int = 1,
) -> None:
    """Save true labels, predictions and associated sequences.

    Args:
        output_dir (Path): Directory in which to save predictions and labels.
        y_true (np.ndarray): Array with true labels. Shape = (N,).
        preds (np.ndarray): Array with predictions. Shape = (N,).
        sequences (List[Tuple[str,str]]): List of sequences used (cdr3b, peptide).
        model_type (str): Type of model used.
        fold_index (int): Index of the fold.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(sequences, columns=["cdr3b", "peptide"])
    df['y_true'] = y_true
    df['y_score'] = preds

    subfolder = output_dir / model_type
    subfolder.mkdir(parents=True, exist_ok=True)
    df.to_csv(str(subfolder / f"fold_{fold_index}.csv"), index=False)

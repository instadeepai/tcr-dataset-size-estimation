"""Train a model on the train set and run inference on the test set"""
from pathlib import Path
from typing import Union

import hydra
import numpy as np
from omegaconf import DictConfig, OmegaConf
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier

from models.utils.constants import MODEL_ACRONYMS
from models.utils.generic import (
    create_classifier_dataset,
    create_model,
    load_config,
    save_labels_predictions
)


def training(
    model: Union[GradientBoostingClassifier, RandomForestClassifier],
    train_set: np.ndarray,
    test_set: np.ndarray,
    output_dir: Path,
    model_type: str,
) -> None:
    """Train a classifier model on the train set and score it on the test set.

    Args:
        model: Model to be trained.
        train_set (np.ndarray): Train set.
        test_set (np.ndarray): Test set.
        output_dir (Path): Directory where to save predictions.
        model_type (str): Type of model used.
    """
    model.fit(train_set.data, train_set.y)

    preds_test = model.predict_proba(test_set.data)[:,1]

    save_labels_predictions(
        output_dir, test_set.y, preds_test, test_set.sequences, MODEL_ACRONYMS[model_type]
    )


@hydra.main(config_path="config/", config_name="config")
def main(args: DictConfig) -> None:
    """Entry point to train a model on the train set and run it on the test set."""
    # Training variables
    OmegaConf.set_struct(args, False)
    conf_data, conf_model = load_config(args)
    train_dir, test_dir, output_dir = (
        Path(conf_data.data_dir) / "train",
        Path(conf_data.data_dir) / "test",
        Path(conf_data.output_dir)
    )
    output_dir.mkdir(exist_ok=True, parents=True)

    train_set = create_classifier_dataset(train_dir)
    test_set = create_classifier_dataset(test_dir)

    model = create_model(conf_model)

    training(
        model,
        train_set,
        test_set,
        output_dir,
        conf_model.model_type
    )

    return None


if __name__ == "__main__":
    main()

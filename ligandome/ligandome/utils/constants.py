"""Constants used for ligandome computations."""

from pathlib import Path

KMER_LENGTH: int = 9
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")
DATABASE_EXPORTS_PATH: Path = Path(__file__).parent.parent / 'database_exports'
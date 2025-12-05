import os
import threading
from typing import Optional

import joblib

# Path to the saved model from train_model.py
# (we expect train_model.py to save to models/interaction_rf.joblib)
MODEL_PATH = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    "models",
    "interaction_rf.joblib",
)

_model = None
_model_lock = threading.Lock()

FEATURE_ORDER = [
    "deltaG_predicted",
    "gc_content_srna",
    "gc_content_mrna",
    "seed_length",
    "has_seed",
    "hybrid_length",
]


def _load_model():
    """
    Lazy-load the trained model once per process.
    Returns None if the model file does not exist.
    """
    global _model
    if _model is not None:
        return _model

    with _model_lock:
        if _model is not None:
            return _model

        if not os.path.exists(MODEL_PATH):
            return None

        try:
            _model = joblib.load(MODEL_PATH)
        except Exception:
            _model = None

    return _model


def is_available() -> bool:
    """
    Returns True if an ML model is available and loaded.
    """
    return _load_model() is not None


def score_interaction(
    deltaG: Optional[float],
    gc_srna: Optional[float],
    gc_mrna: Optional[float],
    seed_length: Optional[float],
    has_seed: Optional[float],
    hybrid_length: Optional[float],
) -> Optional[float]:
    """
    Compute an ML score for a single interaction, using the feature vector:
        [deltaG, gc_srna, gc_mrna, seed_length, has_seed, hybrid_length]

    Returns:
        float score (e.g., 0–1) or None if model/features are unavailable.
    """
    model = _load_model()
    if model is None:
        return None

    # All features must be present
    feats = [deltaG, gc_srna, gc_mrna, seed_length, has_seed, hybrid_length]
    if any(v is None for v in feats):
        return None

    try:
        feat_vec = [[float(f) for f in feats]]
        score = model.predict(feat_vec)[0]
        return float(score)
    except Exception:
        return None

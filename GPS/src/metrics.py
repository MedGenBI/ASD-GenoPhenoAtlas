import numpy as np
import numpy as np
from sklearn.metrics import (
    auc,
    average_precision_score,
    roc_curve,
)


def bootstrap_auc(y_true, y_scores, n_bootstrap=1000, random_state=42):
    """Return scalar bootstrap summaries for AUROC and average precision."""
    y_true = np.asarray(y_true, dtype=int)
    y_scores = np.asarray(y_scores, dtype=float)
    if len(y_true) != len(y_scores):
        raise ValueError("Labels and scores must have the same length.")
    if n_bootstrap < 1:
        raise ValueError("n_bootstrap must be positive.")

    rng = np.random.default_rng(random_state)
    auroc_values = []
    auprc_values = []
    nauprc_values = []

    for _ in range(n_bootstrap):
        indices = rng.choice(len(y_true), len(y_true), replace=True)
        y_boot = y_true[indices]
        if np.unique(y_boot).size < 2:
            continue

        score_boot = y_scores[indices]
        fpr, tpr, _ = roc_curve(y_boot, score_boot)
        auroc_values.append(float(auc(fpr, tpr)))

        ap = float(average_precision_score(y_boot, score_boot))
        auprc_values.append(ap)
        prevalence = float(y_boot.mean())
        if prevalence > 0:
            nauprc_values.append(ap / prevalence)

    def summarize(values):
        if not values:
            return np.nan, np.nan, np.nan, np.nan
        array = np.asarray(values, dtype=float)
        return (
            float(array.mean()),
            float(np.median(array)),
            float(np.percentile(array, 2.5)),
            float(np.percentile(array, 97.5)),
        )

    auroc_mean, auroc_median, auroc_low, auroc_high = summarize(auroc_values)
    auprc_mean, auprc_median, auprc_low, auprc_high = summarize(auprc_values)
    nauprc_mean, nauprc_median, nauprc_low, nauprc_high = summarize(nauprc_values)

    return {
        "auroc_mean": auroc_mean,
        "auroc_median": auroc_median,
        "auroc_ci_lower": auroc_low,
        "auroc_ci_upper": auroc_high,
        "auprc_mean": auprc_mean,
        "auprc_median": auprc_median,
        "auprc_ci_lower": auprc_low,
        "auprc_ci_upper": auprc_high,
        "nauprc_mean": nauprc_mean,
        "nauprc_median": nauprc_median,
        "nauprc_ci_lower": nauprc_low,
        "nauprc_ci_upper": nauprc_high,
        "n_bootstrap_requested": int(n_bootstrap),
        "n_bootstrap_valid": int(len(auroc_values)),
        "n_bootstrap_skipped": int(n_bootstrap - len(auroc_values)),
        "n_nauprc_valid": int(len(nauprc_values)),
        "n_nauprc_skipped": int(n_bootstrap - len(nauprc_values)),
    }

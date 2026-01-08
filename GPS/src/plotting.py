# src/plotting.py

import matplotlib.pyplot as plt
from matplotlib import rcParams
import os


def setup_plot_style():

    rcParams['pdf.fonttype'] = 42
    rcParams['ps.fonttype'] = 42
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['Arial']
    rcParams['axes.linewidth'] = 0.5

def plot_global_roc(
    hc_data: dict, 
    ex_data: dict, 
    output_path: str,
    figsize: tuple = (3.4, 3.14)
):
    setup_plot_style()
    
    plt.figure(figsize=figsize)

    # 1. Plot High-confidence Set
    if hc_data:
        plt.plot(
            hc_data['base_fpr'], 
            hc_data['tpr_mean'], 
            color='#dd666e', 
            linewidth=1,
            label=f"High-confidence Set (AUC = {hc_data['auc_median']:.3f})"
        )

    # 2. Plot Extended Set
    if ex_data:
        plt.plot(
            ex_data['base_fpr'], 
            ex_data['tpr_mean'], 
            color='#64B5F6', 
            linewidth=1,
            label=f"Extended Set (AUC = {ex_data['auc_median']:.3f})"
        )

    # 3. Random classifier
    plt.plot([0, 1], [0, 1], 'k--', alpha=0.5, linewidth=1, label='Random Classifier')

    # Styling
    plt.xlabel('False Positive Rate', fontsize=7, fontweight='regular')
    plt.ylabel('True Positive Rate', fontsize=7, fontweight='regular')
    plt.xticks(fontsize=7)
    plt.yticks(fontsize=7)
    # plt.legend(fontsize=6, loc='lower right')
    plt.legend(fontsize=6)
    plt.grid(True, alpha=0.3)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.tight_layout()

    # Fine-tune spines
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
    ax.tick_params(axis='both', which='major', width=0.5)

    # Save
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    plt.savefig(output_path, format='pdf', dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    print(f"ROC Plot saved to {output_path}")
    plt.close()

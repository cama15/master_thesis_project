import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Load the dataframe
file_path = "/data/cephfs-1/work/groups/kühnen/users/cama15_c/pomc_cpg_data.csv"
df = pd.read_csv(file_path)

# Convert Start and End to integers if necessary
df["Start"] = df["Start"].astype(int)
df["End"] = df["End"].astype(int)

# Output directory
output_dir = "/data/cephfs-1/work/groups/kühnen/users/cama15_c/"

# Define color mapping
color_map = {
    "Naive": "#c51b8a",                    # dark pink
    "Formative": "#ae017e",                # slightly darker pink
    "POMC_neuron_formative": "#7a0177",    # deepest pink-purple
    "Primed": "#238b45",                   # dark green
    "POMC_neuron_primed": "#006d2c",       # deeper green
}

# Default colors
default_methylation_color = "blue"
default_coverage_color = "gray"

cell_states = df["Cell_state"].unique()

# Define the number of rows and columns for the subplot grid
num_plots = len(cell_states)
fig, axes = plt.subplots(num_plots, 1, figsize=(10, 5 * num_plots), sharex=True)

# If there's only one plot, axes might be a single object instead of an array, so handle that case
if num_plots == 1:
    axes = [axes]

axes[0].set_title("Methylation and Sequencing Coverage for Each Cell State at POMC Intron 2-Exon 3 CpG Region")

# Loop through each cell state and create a plot on the respective axes
for i, state in enumerate(cell_states):
    subset = df[df["Cell_state"] == state]
    ax1 = axes[i]  # Select the current axis for this plot

    # Midpoints for x-axis
    x = (subset["Start"] + subset["End"]) / 2
    methylation = subset["Avg_methylation"]
    coverage = subset["Coverage"]

    # Choose colors based on state
    methylation_color = color_map.get(state, default_methylation_color)

    # Plot methylation
    ax1.plot(x, methylation, color=methylation_color, label="Methylation Fraction", linewidth=2)
    ax1.set_ylabel(f"Avg Methylation Fraction\n({state})", color=methylation_color)
    ax1.tick_params(axis="y", labelcolor=methylation_color)

    # Set y-axis limits for methylation to always be between 0 and 1
    ax1.set_ylim(0, 1.0)

    # Secondary y-axis for coverage
    ax2 = ax1.twinx()
    ax2.bar(x, coverage, width=1.5, color=default_coverage_color, alpha=0.4, label="Coverage")
    ax2.set_ylabel("Coverage", color="gray")
    ax2.tick_params(axis="y", labelcolor="gray")

    # Set x-axis ticks and labels for the last plot only
    if i == num_plots - 1:
        ax1.set_xlabel("Genomic Position (chr2:25,161,68-chr2:25,161,767)")
        ax1.set_xticks([25161687, 25161690, 25161700, 25161717, 
          25161721, 25161729, 25161743, 25161765, 25161767])
        ax1.set_xticklabels(["-2", "-1", "+1", "+2", "+3", "+4", "+5", "+6", "+7"])

# Adjust layout to prevent overlap
fig.tight_layout()

# Save the individual plots as a single JPEG
individual_plots_filename = os.path.join(output_dir, "pomc_cpg_by_cell_state_methylation_coverage.jpeg")
plt.savefig(individual_plots_filename, format="jpeg", dpi=300)
plt.close(fig)

# Now create the combined plot with all lines from different cell states

fig_combined, ax_combined = plt.subplots(figsize=(12, 6))

xticks = [25161687, 25161690, 25161700, 25161717, 25161721, 
        25161729, 25161743, 25161765, 25161767
    ]

xticklabels = [
    "-2", "-1",
    "+1", "+2", "+3",
    "+4", "+5", "+6", 
    "+7"
]

for state in cell_states:
    subset = df[df["Cell_state"] == state].dropna(subset=["Avg_methylation"])

    x = (subset["Start"] + subset["End"]) / 2
    methylation = subset["Avg_methylation"]

    # Choose color
    color = color_map.get(state, default_methylation_color)

    # Plot the line
    ax_combined.plot(x, methylation, label=state, color=color, linewidth=2)

# Format combined plot
ax_combined.set_xlabel("Genomic Position (chr2:25,161,65-chr2:25,161,767)")
ax_combined.set_ylabel("Avg Methylation Fraction")
ax_combined.set_ylim(0, 1.0)
ax_combined.set_xticks(xticks)
ax_combined.set_xticklabels(xticklabels)
ax_combined.set_title("Methylation and Sequencing Coverage at the POMC Intron 2-Exon 3 CpG Region")
ax_combined.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))

# Save combined plot
combined_filename = os.path.join(output_dir, "pomc_cpg_methylation_lines_2.jpeg")
fig_combined.tight_layout()
fig_combined.savefig(combined_filename, format="jpeg", dpi=300)

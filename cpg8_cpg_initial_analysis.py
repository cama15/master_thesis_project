import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# Load the dataframe
file_path = "/data/cephfs-1/work/groups/kühnen/users/cama15_c/cpg8_filtered_data.csv"
df = pd.read_csv(file_path)

# Convert Start and End to integers if necessary
df["Start_pos"] = df["Start_pos"].astype(int)
df["End_pos"] = df["End_pos"].astype(int)

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

# Define genomic positions to highlight
highlighted_positions = [11102040, 11102051, 11102086, 11102094, 11102142, 11102150, 11102157, 11102180, 11102204]

highlighted_region_start = 11102040
highlighted_region_end = 11102204

# Filter the data to include only rows within the highlighted region
df_highlighted = df[(df["Start_pos"] >= highlighted_region_start) & (df["End_pos"] <= highlighted_region_end)]


# Define the number of rows and columns for the subplot grid
num_plots = len(cell_states)
fig, axes = plt.subplots(num_plots, 1, figsize=(10, 5 * num_plots), sharex=True)

# If there's only one plot, axes might be a single object instead of an array
if num_plots == 1:
    axes = [axes]

axes[0].set_title("Methylation and Sequencing Coverage at the Hypermethylated Control CpG 8 Region")

# Loop through each cell state and create a plot
for i, state in enumerate(cell_states):
    subset = df_highlighted[df_highlighted["Cell_state"] == state]
    ax1 = axes[i]  # Select the current axis for this plot

    # Midpoints for x-axis
    x = (subset["Start_pos"] + subset["End_pos"]) / 2
    methylation = subset["Meth_fraction"]
    coverage = subset["Coverage"]

    # Choose colors
    methylation_color = color_map.get(state, default_methylation_color)

    # Plot methylation line
    ax1.plot(x, methylation, color=methylation_color, label="Methylation Fraction", linewidth=2)
    ax1.set_ylabel(f"Methylation Fraction\n({state})", color=methylation_color)
    ax1.tick_params(axis="y", labelcolor=methylation_color)
    ax1.set_ylim(0, 1.0)

    # Secondary y-axis for coverage
    ax2 = ax1.twinx()
    ax2.bar(x, coverage, width=1.5, color=default_coverage_color, alpha=0.4, label="Coverage")
    ax2.set_ylabel("Coverage", color="gray")
    ax2.tick_params(axis="y", labelcolor="gray")

    # Set x-axis ticks and labels for the last plot only
    if i == num_plots - 1:
        ax1.set_xlabel("Genomic Position (chr19:11102040-11102206)")
        ax1.set_xticks(highlighted_positions)
        ax1.set_xticklabels([str(pos) for pos in highlighted_positions], rotation=45, ha="right")

# Adjust layout to prevent overlap
fig.tight_layout()

# Save the individual plots as a single JPEG
individual_plots_filename = os.path.join(output_dir, "cpg8_by_cell_state_methylation_coverage.jpeg")
plt.savefig(individual_plots_filename, format="jpeg", dpi=300)
plt.close(fig)

# Now create the combined plot with all lines from different cell states
fig_combined, ax_combined = plt.subplots(figsize=(12, 6))

for state in cell_states:
    subset = df_highlighted[df_highlighted["Cell_state"] == state].dropna(subset=["Meth_fraction"])

    x = (subset["Start_pos"] + subset["End_pos"]) / 2
    methylation = subset["Meth_fraction"]

    # Choose color
    color = color_map.get(state, default_methylation_color)

    # Plot the line
    ax_combined.plot(x, methylation, label=state, color=color, linewidth=2)

# Format combined plot
ax_combined.set_xlabel("Genomic Position (chr19:11102040-11102206)")
ax_combined.set_ylabel("Methylation Fraction")
ax_combined.set_ylim(0, 1.0)
ax_combined.set_xticks(highlighted_positions)
ax_combined.set_xticklabels([str(pos) for pos in highlighted_positions], rotation=45, ha="right")
ax_combined.set_title("Methylation and Sequencing Coverage at the Hypermethylated Control CpG 10 Region")
ax_combined.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))

# Save combined plot
combined_filename = os.path.join(output_dir, "cpg8_methylation_lines.jpeg")
fig_combined.tight_layout()
fig_combined.savefig(combined_filename, format="jpeg", dpi=300)

"""
==============================================================================
Date: February 11, 2025
Version: 1.0

Description:
This script generates visualizations and 3D models for the manuscript by 
Putti et al. (in preparation). The script includes:
- Data processing and statistical analysis of synaptic distributions.
- Visualization of neuron functional morphotypes and neurotransmitter classifications.
- Identification and display of significantly active neurons based on functional 
  imaging data.
- Import and visualization of neuron meshes in Blender, with appropriate color 
  coding based on neuronal classification.

Dependencies:
- Python 3.8+
- Conda environment (install using: `conda env create --file make_videos.yaml`)
- Required Python packages: h5py, pandas, numpy, matplotlib, scipy, seaborn, 
  pathlib, Blender (for 3D visualization)

Usage:
- Modify file paths in the `ROOT_PATH` variable to match your local setup.
- Run the script in Python for statistical analysis and figure generation.
- Run Blender-specific sections in the Blender Python terminal for 3D visualization.

File Structure:
- The script is divided into sections:
  1. Data Processing and Loading
  2. Axo-Type Distributions Analysis
  3. Functional Morphotypes Visualization
  4. Neurotransmitter Classification
  5. Active Neuron Detection and Visualization
  6. Blender Import and 3D Visualization

Output:
- Figures are saved as high-resolution PDFs in the `ROOT_PATH` directory.
- Active neuron data and significance values are returned for further analysis.
- Blender-imported neuron meshes are colored and labeled based on their functional 
  classification.

License:
- This code is intended for research purposes in conjunction with the 
  manuscript by Putti et al. (in preparation). Redistribution or modification 
  should be discussed with the author.

Contact:
Jonathan Boulanger-Weill  
Harvard University & Institut de la Vision, Sorbonne Université Paris 
Email: jonathanboulangerweill@harvard.edu  

==============================================================================
"""

import os
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy import stats

# Define file paths
ROOT_PATH = Path("/Users/jonathanboulanger-weill/Dropbox (Harvard University)/hb_connectome/hindbrain_structure_function/clem_zfish1/collab_rgcs/")
HDF5_PATH = "/Users/jonathanboulanger-weill/Dropbox (Harvard University)/hb_connectome/hindbrain_structure_function/clem_zfish1/function/all_cells_bars.h5"
EXCEL_FILE_PATH = ROOT_PATH / "rgc_axons_output_020525.csv"
SYP_FILE_PATH = ROOT_PATH / "rgc_axons_syp_020525.csv"
OUTPUT_FILE = ROOT_PATH / "rgcs_outputs.pdf"

# Set global font to Arial
plt.rcParams['font.family'] = 'Arial'

# Define RGC IDs
rgcs_ids = ['576460752701268702', '576460752654570799', '576460752665937704', '576460752643481529']

# Load datasets
df_traced = pd.read_csv(EXCEL_FILE_PATH, dtype=str)  # Traced synapses
df_total = pd.read_csv(SYP_FILE_PATH, dtype=str)  # Total synapses (traced + untraced)

# --------------------------- 1st PLOT: Axo-Type Distributions ---------------------------

# Define axo-type categories
axo_categories = ['axo-axonic', 'axo-somatic']
default_category = 'axo-dendritic'

# Extract relevant data
df_traced['connectivity'] = df_traced['connectivity'].astype(str).str.split().str[0]
filtered_traced_df = df_traced[df_traced['connectivity'].isin(rgcs_ids)]

# Count occurrences per RGC for traced synapses
axo_counts = {rgc: {axo: 0 for axo in axo_categories + [default_category]} for rgc in rgcs_ids}

for rgc in rgcs_ids:
    rgc_subset = filtered_traced_df[filtered_traced_df['connectivity'] == rgc]
    
    for axo in axo_categories:
        axo_counts[rgc][axo] = rgc_subset['comment'].str.contains(axo, na=False).sum()
    
    # Assign remaining cases to axo-dendritic
    total_rows = len(rgc_subset)
    specific_axo_count = sum(axo_counts[rgc][axo] for axo in axo_categories)
    axo_counts[rgc][default_category] = total_rows - specific_axo_count

axo_type_df = pd.DataFrame.from_dict(axo_counts, orient='index')

# Compute total synapses (traced + untraced)
df_total['RGC'] = df_total['RGC'].astype(str).str.split().str[0]
rgc_counts = {rgc: df_total[df_total['RGC'] == rgc].shape[0] for rgc in rgcs_ids}
rgc_count_df = pd.DataFrame(list(rgc_counts.items()), columns=['RGC ID', 'Total Synapses']).set_index('RGC ID')

# Compute untraced synapses
rgc_count_df['Untraced Synapses'] = rgc_count_df['Total Synapses'] - axo_type_df.sum(axis=1)
rgc_count_df.drop(columns=['Total Synapses'], inplace=True)

# Merge traced and untraced synapse counts
combined_df = axo_type_df.join(rgc_count_df, how='left')

# Define colors
axo_type_colors = {
    'axo-axonic': '#E64B35',  # Red
    'axo-somatic': '#F39C12',  # Orange
    'axo-dendritic': '#1F77B4',  # Blue
    'Untraced Synapses': '#000000'  # Black
}

# --------------------------- 2nd PLOT: Functional Morphotypes ---------------------------

# Define output types
output_types = ['axo-axonic', 'nsPVIN', 'nsPVPN', 'bsPVIN', 'bsPVPN', 'msPVIN', 'msPVPN', 'tsPVIN', 'tsPVPN', 'unknown']

# Define colors
functional_colors_dict = {
    'axo-axonic': '#000000', 'nsPVIN': '#F39C12', 'nsPVPN': '#F39C12', 'bsPVIN': '#1F77B4',
    'bsPVPN': '#1F77B4', 'msPVIN': '#D62728', 'msPVPN': '#D62728', 'tsPVIN': '#2CA02C', 'tsPVPN': '#2CA02C',
    'unknown': '#7F8C8D'
}

df_traced['connectivity_id'] = df_traced['connectivity'].astype(str).str.split().str[0]
filtered_df = df_traced[df_traced['connectivity_id'].isin(rgcs_ids)].drop_duplicates(subset='dendrite_id')

# Count functional morphotypes
functional_morphotype_counts = {rgc: {output: 0 for output in output_types} for rgc in rgcs_ids}

for rgc in rgcs_ids:
    rgc_subset = filtered_df[filtered_df['connectivity_id'] == rgc]
    for output in output_types:
        functional_morphotype_counts[rgc][output] = rgc_subset['comment'].str.contains(output, na=False).sum()

functional_morphotype_df = pd.DataFrame.from_dict(functional_morphotype_counts, orient='index')

# For clarify, set all values = 1 to NaN
# functional_morphotype_df[functional_morphotype_df <= 1] = 0

# --------------------------- 3rd PLOT: Neurotransmitter Classifications ---------------------------

neurotransmitter_labels = ['Excitatory', 'Inhibitory', 'Unknown']
neurotransmitter_colors = ['#27AE60', '#E91E63', '#7F8C8D']  # Green, Pink, Gray

# Count neurotransmitters
neurotransmitter_counts = {rgc: {'excitatory': 0, 'inhibitory': 0, 'unknown': 0} for rgc in rgcs_ids}

for rgc in rgcs_ids:
    rgc_subset = filtered_df[filtered_df['connectivity_id'] == rgc]
    neurotransmitter_distribution = rgc_subset['neurotransmitter classifier'].value_counts().to_dict()
    neurotransmitter_counts[rgc]['excitatory'] = neurotransmitter_distribution.get('excitatory', 0)
    neurotransmitter_counts[rgc]['inhibitory'] = neurotransmitter_distribution.get('inhibitory', 0)
    neurotransmitter_counts[rgc]['unknown'] = neurotransmitter_distribution.get('unknown', 0)

neurotransmitter_df = pd.DataFrame.from_dict(neurotransmitter_counts, orient='index')

# --------------------------- PLOTTING ---------------------------

# Create figure with specified height (41.2 mm = 1.622 inches)
fig, axes = plt.subplots(1, 4, figsize=(8, 1.622))  # Width: 9 inches, Height: 1.622 inches

# X positions
x_positions = np.arange(len(rgcs_ids))
bar_width = 0.7

# Plot Axo-Type Distributions
bottom_values = np.zeros(len(rgcs_ids))
for category in axo_type_colors:
    axes[0].bar(x_positions, combined_df[category], width=bar_width, label=category,
                color=axo_type_colors[category], alpha=0.6, bottom=bottom_values)
    bottom_values += combined_df[category]

axes[0].set_title("Axo-Type Distributions", fontsize=7)
axes[0].set_ylabel("Number of synapses", fontsize=7)
axes[0].set_xlabel("Input RGCs", fontsize=7)
axes[0].legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.05, 1), frameon=False)

# Plot Functional Morphotypes
bottom_values = np.zeros(len(rgcs_ids))
legend_handles_dict = {}  # Dictionary to store unique legend handles

for output in output_types:
    values = functional_morphotype_df[output].values  # Extract values for current category
    mask = values > 0  # Mask to select only RGCs with nonzero values

    if np.any(mask):  # Only proceed if at least one value is nonzero
        if "PVPN" in output:
            bar = axes[1].bar(x_positions[mask], values[mask], width=bar_width,
                               edgecolor=functional_colors_dict[output], fill=False, 
                               bottom=bottom_values[mask], linewidth=1)
        else:
            bar = axes[1].bar(x_positions[mask], values[mask], width=bar_width, 
                              edgecolor=functional_colors_dict[output], 
                              color=functional_colors_dict[output], alpha=0.7, bottom=bottom_values[mask])

        bottom_values[mask] += values[mask]  # Update stacking only for non-zero positions
        
        # Add to legend dictionary only once
        if output not in legend_handles_dict:
            legend_handles_dict[output] = bar[0]

# Convert dictionary values (bars) to a list for the legend
legend_handles = list(legend_handles_dict.values())

axes[1].set_title("Functional Morphotypes", fontsize=7)
axes[1].set_ylabel("Number of neurons", fontsize=7)
axes[1].set_xlabel("Input RGCs", fontsize=7)
# Add legend to the plot
axes[1].legend(legend_handles, legend_handles_dict.keys(), fontsize=6, loc="upper left", bbox_to_anchor=(1.05, 1))

# Plot Neurotransmitter Classifications
bottom_values = np.zeros(len(rgcs_ids))
for i, neurotransmitter in enumerate(neurotransmitter_labels):
    axes[2].bar(x_positions, neurotransmitter_df[neurotransmitter.lower()], width=bar_width,
                color=neurotransmitter_colors[i], alpha=0.8, bottom=bottom_values)
    bottom_values += neurotransmitter_df[neurotransmitter.lower()]

axes[2].set_title("Neurotransmitter Classification", fontsize=7)
axes[2].set_ylabel("Number of neurons", fontsize=7)
axes[2].set_xlabel("Input RGCs", fontsize=7)
axes[2].legend(neurotransmitter_labels, fontsize=6, loc="upper left", bbox_to_anchor=(1.05, 1), frameon=False)

# New Histogram for Layer RGCs
layer_rgcs = [4, 6, 4, 12]  # Given values

axes[3].bar(x_positions, layer_rgcs, width=bar_width, color='gray', alpha=0.8)
axes[3].set_title("Layer RGCs", fontsize=7)
axes[3].set_ylabel("Number of RGCs", fontsize=7)
axes[3].set_xlabel("Layer", fontsize=7)

# Match y-axis for plots 2 and 3
max_y = max(functional_morphotype_df.max().sum(), neurotransmitter_df.max().sum())
axes[1].set_ylim(0, max_y)
axes[2].set_ylim(0, max_y)

for ax in axes:
    ax.set_xticks(x_positions)
    ax.set_xticklabels([1, 2, 3, 4], fontsize=7)
    ax.tick_params(axis='both', labelsize=7)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set width of left and bottom spines to 0.8
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)

    # Set tick width to 0.8 for both major and minor ticks
    ax.tick_params(axis='both', labelsize=7, width=0.8, length=3)  # length=3 keeps them readable

plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches='tight', format='pdf')
plt.show()

# --------------------------- 4th PLOT: Activity ---------------------------

# Color dictionary
COLOR_CELL_TYPE_DICT_LINES = {
    "bsPVPN": (81/255, 54/255, 187/255, 1),   # Deep Blue (#5136BB)
    "nsPVIN": (243/255, 156/255, 18/255, 1),  # Bright Orange (#F39C12)
    "nsPVPN": (248/255, 125/255, 210/255, 1),  # Vivid Orange (#F87D15)
    "undetermined": (127/255, 140/255, 141/255, 1)  # Muted Gray (#7F8C8D)
}

def plot_neuron_activity(file_path, neuron_idx):
    """
    Plots the average and individual trial activity dynamics for a specified neuron.

    Parameters:
    - file_path: str, path to the HDF5 file containing neuron data
    - neuron_idx: int, index of the neuron to plot
    """
    # Open the HDF5 file
    with h5py.File(file_path, "r") as hdf_file:
        # Access the group for the specified neuron
        neuron_group = hdf_file[f"neuron_{neuron_idx}"]

        # Get average activity for left and right stimuli
        avg_activity_left = neuron_group["average_activity_left"][()]
        avg_activity_right = neuron_group["average_activity_right"][()]
        
        # Get individual trial data for left and right stimuli
        trials_left = neuron_group["neuronal_activity_trials_left"][()]
        trials_right = neuron_group["neuronal_activity_trials_right"][()]

        # Smooth using a Savitzky-Golay filter
        smooth_avg_activity_left = savgol_filter(avg_activity_left, 20, 3)
        smooth_avg_activity_right = savgol_filter(avg_activity_right, 20, 3)
        smooth_trials_left = savgol_filter(trials_left, 20, 3, axis=1)
        smooth_trials_right = savgol_filter(trials_right, 20, 3, axis=1)

        # Define time axis in seconds
        dt = 0.5  # Time step is 0.5 seconds
        time_axis = np.arange(len(avg_activity_left)) * dt

        # Plot smoothed average activity with thick lines
        fig, ax = plt.subplots()
        plt.plot(time_axis, smooth_avg_activity_left, color='blue', alpha=0.7, linewidth=3, label='Smoothed Average Left')
        plt.plot(time_axis, smooth_avg_activity_right, color='red', alpha=0.7, linestyle='--', linewidth=3, label='Smoothed Average Right')
        
        # Plot individual trial data with thin black lines for left and dashed black lines for right
        for trial_left, trial_right in zip(smooth_trials_left, smooth_trials_right):
            plt.plot(time_axis, trial_left, color='black', alpha=0.3, linewidth=1)
            plt.plot(time_axis, trial_right, color='black', alpha=0.3, linestyle='--', linewidth=1)
        
        # Overlay shaded rectangle for stimulus epoch
        plt.axvspan(20, 60, color='gray', alpha=0.1, label='Stimulus Epoch')

        cluster_ID=neuron_group["cluster"][()]
        cluster_vals=neuron_group["cluster_values"][()]
        cluster_trials=neuron_group["cluster_values_trials"][()]
        cluster_max_trials=neuron_group["cluster_max_trials"][()]

        plt.title(f'Average and Individual Trial for Neuron {neuron_idx} \n (C: {cluster_ID}) \n (V: {cluster_vals}) \n (VT: {cluster_trials}) \n (VT: {cluster_max_trials})' )
        plt.xlabel('Time (seconds)')
        plt.ylabel('Activity')
        
        # Set font of legend text to Arial
        legend = plt.legend()
        for text in legend.get_texts():
            text.set_fontfamily('Arial')
            
        # Set aspect ratio to 1
        ratio = 1.0
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()
        ax.set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)
        plt.show()

def plot_active_neurons(cell_ids, cell_types, hdf5_file_path, cell_type, out_dir, dt=0.5, alpha_threshold=0.05):
    with h5py.File(hdf5_file_path, "r") as hdf_file:
        # Set figure size (height fixed to 1.564 inches, adjust width proportionally)
        fig, ax = plt.subplots(figsize=(4, 1.564))  # Adjust width (4 in) as needed
        
        # Set font to Arial globally
        plt.rcParams["font.family"] = "Arial"
        plt.rcParams["axes.linewidth"] = 0.8  # Set axis thickness

        # Initialize arrays
        all_max_arrays = []
        active_neuron_indices = []  # Store indices of significant neurons
        active_neuron_pvalues = []  # Store p-values of significant neurons

        # Map cell IDs to their types
        cell_id_to_type = {cell_ids[i]: cell_types[i] for i in range(len(cell_ids))}

        for i, idx in enumerate(cell_ids):
            neuron_group = hdf_file[f"neuron_{idx}"]
            avg_activity_left = neuron_group["average_activity_left"][()]
            avg_activity_right = neuron_group["average_activity_right"][()]
            smooth_avg_activity_left = savgol_filter(avg_activity_left, 20, 3)
            smooth_avg_activity_right = savgol_filter(avg_activity_right, 20, 3)

            # Find the maximum values in each array
            max_left = np.max(smooth_avg_activity_left[40:120])
            max_right = np.max(smooth_avg_activity_right[40:120])

            # Determine which array has the maximum value
            max_array = smooth_avg_activity_left if max_left > max_right else smooth_avg_activity_right

            # Pool activity 
            all_max_arrays.append(max_array)

        # Convert list to NumPy array
        all_max_arrays = np.array(all_max_arrays)  # Shape: (num_neurons, num_timepoints)
        
        # Define Time Windows (indices)
        baseline_window = all_max_arrays[:, 0:40]  # Pre-stimulus (20s before stimulus)
        early_stim_window = all_max_arrays[:, 40:80]  # 20s early during stimulus
        late_stim_window = all_max_arrays[:, 80:120]  # 20s late during stimulus

        # --- Identify Significant Neurons ---
        for i in range(len(cell_ids)):
            baseline = baseline_window[i]
            early_stim = early_stim_window[i]
            late_stim = late_stim_window[i]

            # Check normality for early_stim vs baseline
            _, p_early_normality = stats.shapiro(early_stim - baseline)
            _, p_late_normality = stats.shapiro(late_stim - baseline)

            # Choose test based on normality
            if p_early_normality > 0.05:
                _, p_early = stats.ttest_rel(early_stim, baseline, alternative='greater')
            else:
                _, p_early = stats.wilcoxon(early_stim, baseline, alternative='greater')

            if p_late_normality > 0.05:
                _, p_late = stats.ttest_rel(late_stim, baseline, alternative='greater')
            else:
                _, p_late = stats.wilcoxon(late_stim, baseline, alternative='greater')

            # If neuron is significantly active in either phase, include it
            if p_early < alpha_threshold or p_late < alpha_threshold:
                active_neuron_indices.append(i)
                active_neuron_pvalues.append((cell_ids[i], p_early, p_late))

        # --- Plot Only Significant Neurons ---
        num_active_neurons = len(active_neuron_indices)

        if num_active_neurons > 0:
            legend_handles = []
            used_types = set()  # To track which types are already in legend

            for i in active_neuron_indices:
                neuron_type = cell_id_to_type[cell_ids[i]]
                color = COLOR_CELL_TYPE_DICT_LINES[neuron_type]  # Assign color by neuron type
                time_axis = np.arange(len(all_max_arrays[i])) * dt
                
                plt.plot(time_axis, all_max_arrays[i], color=color, alpha=0.7, linestyle='-', linewidth=1)

                # Add to legend only once per type
                if neuron_type not in used_types:
                    legend_handles.append(plt.Line2D([0], [0], color=color, lw=2, label=neuron_type))
                    used_types.add(neuron_type)

            # Compute and plot the mean trace of active neurons
            mean_active_array = np.nanmean(all_max_arrays[active_neuron_indices], axis=0)
            plt.plot(time_axis, mean_active_array, color='black', alpha=0.7, linestyle='-', linewidth=2, label="Mean Active Response")

            # Add legend
            plt.legend(handles=legend_handles, loc='upper right', fontsize=8, frameon=False)

        # Overlay shaded rectangle for stimulus epoch
        plt.axvspan(20, 60, color='gray', alpha=0.1)

        # Add a dashed line at y=0 during stimulation
        plt.axhline(y=0, xmin=20/len(time_axis), xmax=60/len(time_axis), color='black', linestyle='dashed', linewidth=0.8)
        
        # Customize axis
        ax.spines['top'].set_linewidth(0.8)
        ax.spines['right'].set_linewidth(0.8)
        ax.spines['left'].set_linewidth(0.8)
        ax.spines['bottom'].set_linewidth(0.8)

        # Calculate correct scale bar height (20% of normalized scale)
        scale_bar_height = 20  # Since normalization is in %

        # Remove the axes and add the scale bars
        ax.plot([0, 10], [-5, -5], color='k', lw=0.8)  # Time scale bar (10 sec)
        ax.text(5, -7, '10 sec', ha='center', fontfamily='Arial', fontsize=8)

        # Adapted scale bar for normalized activity (now using 20%)
        ax.plot([-2, -2], [0, scale_bar_height], color='k', lw=0.8)  
        ax.text(-2.5, scale_bar_height / 2, '20%', va='center', fontfamily='Arial', rotation=90, fontsize=8)

        # Set aspect ratio to 1 and remove the axis lines
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()
        ax.set_aspect(abs((x_right - x_left) / (y_low - y_high)))
        ax.set_axis_off()  # Remove the axis

        # Set plot title with active neuron count
        plt.title(f"Active Neurons ({num_active_neurons}) - {cell_type}", fontsize=10, fontweight='bold')
        plt.show()

        # Save the figure
        filename = f"{cell_type}_active_neurons.pdf"
        output_path = os.path.join(out_dir, filename)
        fig.savefig(output_path, dpi=1200, bbox_inches='tight')
        print(f"Figure saved successfully at: {output_path}")

        return num_active_neurons, active_neuron_pvalues

# Single neuron activity     
neuron_idx = 12432
plot_neuron_activity(HDF5_PATH, neuron_idx) 

# Significantly active neurons    
SPV_outputs_idx = [5759, 4949, 15015, 10099, 11672, 10764, 12432, 8958, 13673, 7258, 5506, 2551, 13544, 7237] 
SPV_outputs_types = ['nsPVIN', 'nsPVIN', 'nsPVIN', 'nsPVIN', 'nsPVIN', 'nsPVPN', 'undetermined', 'nsPVIN', 'bsPVPN', 'undetermined', 'undetermined', 'nsPVIN', 'bsPVPN', 'bsPVPN'] 

plot_active_neurons(SPV_outputs_idx, SPV_outputs_types, HDF5_PATH, "nspvin_bars", ROOT_PATH, dt=0.5, alpha_threshold=0.01)

# --------------------------- BLENDER IMPORT ---------------------------
# %%

# Note: this code works only in the Blender Terminal 

import os
import bpy
import pandas as pd

# Define your color dictionary for cell types
color_cell_type_dict = {
    "SAC": (231/255, 154/255, 221/255, 0.7),    # Yellow-orange
    "SGC": (233/255, 170/255, 160/255, 0.7),    # Peach
    "SFGS": (153/255, 204/255, 236/255, 0.7),   # Light blue
    "SO": (55/255, 104/255, 98/255, 0.7),       # Dark green
    "nsPVIN": (220/255, 130/255, 10/255, 0.7),  # Darker Orange-Yellow
    "nsPVPN": (255/255, 180/255, 50/255, 0.7),  # Lighter Golden Yellow
    "bsPVIN": (31/255, 119/255, 180/255, 0.7)   # Gray for unknown types
}

def pull_segments_ids(df, search_strings, column_to_search=3):
    """
    Function to search for substrings in a specific column of a DataFrame
    and return the appropriate column for the filtered rows.

    Parameters:
    - df: pandas DataFrame to search.
    - search_strings: list of substrings to search for in the DataFrame.
    - column_to_search: index of the column to search for substrings (default is 3).

    Returns:
    - results: a dictionary where keys are the search strings and values are the found segment_ids.
    """
    results = {}

    for search_string in search_strings:
        # Filter rows where the search string is present
        filtered_df = df[df.iloc[:, column_to_search].str.contains(search_string, na=False)]

        # Choose the display column dynamically
        segment_ids = []
        for _, row in filtered_df.iterrows():
            if row.iloc[0] == "cell":  # If first column is 'cell', use column 5
                segment_ids.append(row.iloc[5])
            elif row.iloc[0] == "axon":  # If first column is 'axon', use column 7
                segment_ids.append(row.iloc[7])

        results[search_string] = segment_ids  # Store as list for easy access

    return results

def create_and_assign_material(obj, color):
    if obj is None:
        print(f"Object not found")
        return
    
    # Create a new material
    mat_name = f"{obj.name}_Material"
    mat = bpy.data.materials.new(name=mat_name)
    
    # Disable nodes to assign color directly
    mat.use_nodes = False
    
    # Set the diffuse color directly (RGBA format)
    mat.diffuse_color = color
    
    # Assign the material to the object
    if len(obj.data.materials) == 0:
        obj.data.materials.append(mat)  # If no materials, add the new one
    else:
        obj.data.materials[0] = mat  # Replace the existing material with the new one

# Function to load and color cell meshes based on segment IDs found

def find_and_load_cell_meshes_with_colors(root_dir, segment_ids_found, color_cell_type_dict, default_color=(0.5, 0.5, 0.5, 1.0)):
    """
    Load axon, dendrite, and soma meshes based on segment_ids and associate colors based on the search_string.
    This function will load the meshes into Blender and assign colors directly.
    
    Parameters:
        root_dir (str): Root directory containing neuron mesh files.
        segment_ids_found (dict): Dictionary where keys are search strings, and values are lists of segment IDs.
        color_cell_type_dict (dict): Dictionary mapping search strings to RGBA colors.
        default_color (tuple): RGBA default color in case search_string has no assigned color.
    """

    # Ensure root directory is valid
    if not os.path.exists(root_dir):
        print(f"Error: The root directory {root_dir} does not exist.")
        return

    # Walk through the directory tree
    for root, dirs, files in os.walk(root_dir):
        for folder_name in dirs:
            # Extract the segment_id from the folder name (assuming it's at the end)
            folder_segment_id = folder_name.split('_')[-1]
            
            try:
                folder_segment_id = int(folder_segment_id)
            except ValueError:
                continue  # Skip if folder name doesn't contain a valid segment_id

            # Find matching segment ID in segment_ids_found
            assigned_color = None
            for search_string, segment_ids in segment_ids_found.items():
                if folder_segment_id in segment_ids:
                    assigned_color = color_cell_type_dict.get(search_string, default_color)
                    folder_path = os.path.join(root, folder_name)

                    # Load meshes (axon, dendrite, soma) if they exist
                    for part in ["axon", "dendrite", "soma"]:
                        obj_file = f"{folder_name}_{part}.obj"
                        obj_path = os.path.join(folder_path, obj_file)

                        if os.path.exists(obj_path):
                            bpy.ops.wm.obj_import(filepath=obj_path)
                            print(f"Loaded {part} mesh: {obj_file}")

                            # Assign material color
                            imported_obj = bpy.context.selected_objects[0]
                            create_and_assign_material(imported_obj, assigned_color)

                    # Break after processing the matching folder
                    break

# Example usage
root_cells = '/Users/jonathanboulanger-weill/Dropbox (Harvard University)/hb_connectome/hindbrain_structure_function/clem_zfish1/collab_rgcs/traced_axons_neurons'

search_strings = ['SAC', 'SFGS', 'SO', 'SGC', 'nsPVIN', 'nsPVPN', 'bsPVIN']

path_all_cells = '/Users/jonathanboulanger-weill/Dropbox (Harvard University)/hb_connectome/hindbrain_structure_function/clem_zfish1/xls_spreadsheets/rgcs_elena.xlsx'
df_all_cells = pd.read_excel(path_all_cells)

segment_ids_found = pull_segments_ids(df_all_cells, search_strings)

# Load the meshes and apply colors directly in Blender
find_and_load_cell_meshes_with_colors(root_cells, segment_ids_found, color_cell_type_dict)

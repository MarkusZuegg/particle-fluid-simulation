# Fully created by Google Ai studio, for grpahing only
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import sys

# --- Script Configuration ---
CSV_FILE_PATH = 'data.csv'

# Choose which methods to plot (leave empty to plot all)
methods_to_plot = ['SoA', 'SoA grid', 'AoS grid', 'O-1', 'O-2']

# Choose which x-values to exclude
x_values_to_exclude = []
# --------------------------

# --- Data Loading and Processing ---
try:
    all_data = pd.read_csv(CSV_FILE_PATH)
except FileNotFoundError:
    print(f"Error: The file '{CSV_FILE_PATH}' was not found.")
    sys.exit(1)

# Stage 1: Filter by method name
if methods_to_plot:
    filtered_data = all_data[all_data['method'].isin(methods_to_plot)].copy()
else:
    filtered_data = all_data.copy()

# Stage 2: Filter by excluding specified x-values
if x_values_to_exclude:
    filtered_data = filtered_data[~filtered_data['#particles'].isin(x_values_to_exclude)]

if filtered_data.empty:
    print("Error: No data remains after filtering.")
    sys.exit(1)

# Calculate stats
run_columns = ['run1', 'run2', 'run3']
filtered_data['mean'] = filtered_data[run_columns].mean(axis=1)
filtered_data['std_dev'] = filtered_data[run_columns].std(axis=1)

# --- Plotting ---
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(12, 8))
markers = {'SoA': 'o', 'SoA grid': 's', 'AoS grid': '^', 'O1': 'D', 'O2': 'X'}

# Loop through each method to plot it
for method_name, group in filtered_data.groupby('method'):

    # === NEW: Robustness Check ===
    # Check if there is any valid data to plot for this group before proceeding.
    # .notna() finds non-empty values, .any() checks if there's at least one.
    if group['mean'].notna().any():
        ax.errorbar(
            group['#particles'],
            group['mean'],
            yerr=group['std_dev'],
            fmt='-' + markers.get(method_name, 'd'),
            capsize=5,
            label=method_name
        )
    else:
        # If no valid data, print a warning and skip this method.
        print(f"Warning: Skipping method '{method_name}' as it contains no valid data to plot.")
    # ==============================


# --- Final Plot Customization ---
ax.set_xscale('log')
ax.set_yscale('log')

unique_x_values = sorted(filtered_data['#particles'].unique())
ax.set_xticks(unique_x_values)
ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')


ax.grid(True, which='major', linestyle='-', linewidth='0.5', color='gray')
ax.grid(True, which='minor', linestyle=':', linewidth='0.5', color='lightgray')
ax.minorticks_on()

ax.set_title('Performance Comparison of Particle Implementations', fontsize=16)
ax.set_xlabel('# of Particles (log scale)', fontsize=12)
ax.set_ylabel('Performance Metric / Time (log scale)', fontsize=12)
ax.legend(fontsize=10)

plt.tight_layout()
plt.show()
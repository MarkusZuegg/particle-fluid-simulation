import matplotlib.pyplot as plt
import matplotlib.ticker as mticker # Import the ticker module for formatting
import numpy as np
import pandas as pd
import sys

# --- Script Configuration ---
CSV_FILE_PATH = 'data.csv'

# V V V EDIT THIS LIST TO CHOOSE WHICH METHODS TO PLOT V V V
methods_to_plot = ['SoA grid', 'AoS grid']
# --------------------------

# --- Data Loading and Processing ---
try:
    all_data = pd.read_csv(CSV_FILE_PATH)
except FileNotFoundError:
    print(f"Error: The file '{CSV_FILE_PATH}' was not found.")
    sys.exit(1)

# Filter the data based on the methods_to_plot list
if methods_to_plot:
    filtered_data = all_data[all_data['method'].isin(methods_to_plot)].copy()
else:
    filtered_data = all_data.copy()

if filtered_data.empty:
    print("Error: None of the specified methods were found in the data file.")
    sys.exit(1)

run_columns = ['run1', 'run2', 'run3']
filtered_data['mean'] = filtered_data[run_columns].mean(axis=1)
filtered_data['std_dev'] = filtered_data[run_columns].std(axis=1)

# --- Plotting ---
plt.style.use('seaborn-v0_8-whitegrid')
fig, ax = plt.subplots(figsize=(12, 8)) # Made figure slightly larger for clarity
markers = {'SoA': 'o', 'SoA grid': 's', 'AoS grid': '^'}

for method_name, group in filtered_data.groupby('method'):
    ax.errorbar(
        group['#particles'],
        group['mean'],
        yerr=group['std_dev'],
        fmt='-' + markers.get(method_name, 'd'),
        capsize=5,
        label=method_name
    )

# --- Final Plot Customization ---

# Set scales to logarithmic
ax.set_xscale('log')
ax.set_yscale('log')

# === NEW: Detailed Axis Configuration ===

# 1. Force all unique x-values to be shown as major ticks
unique_x_values = sorted(filtered_data['#particles'].unique())
ax.set_xticks(unique_x_values)

# 2. Format the x-axis tick labels as plain integers, not scientific notation
ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter())
# To prevent labels from overlapping, you can rotate them slightly
# plt.setp(ax.get_xticklabels(), rotation=30, horizontalalignment='right')


# 3. Add a more detailed grid by enabling it for both major and minor ticks
# The `which='both'` parameter is key here.
ax.grid(True, which='major', linestyle='-', linewidth='0.5', color='gray')
ax.grid(True, which='minor', linestyle=':', linewidth='0.5', color='lightgray')

# Ensure minor ticks are visible (they are needed for the minor grid)
ax.minorticks_on()

# =======================================

# Add titles, labels, and legend
ax.set_title('Performance Comparison of Particle Implementations', fontsize=16)
ax.set_xlabel('# of Particles (log scale)', fontsize=12)
ax.set_ylabel('Time microseconds (log scale)', fontsize=12)
ax.legend(fontsize=10)

# Display the plot
plt.tight_layout()
plt.show()
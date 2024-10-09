#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import glob
import os

# Set the directory where the .xlsx files are located
data_dir = "/fh/fast/_IRC/FHIL/user/dgratz/BM02/data"

# Get a list of all .xlsx files in the directory
file_list = glob.glob(os.path.join(data_dir, "*.xlsx"))

# Define colors for each prefix
color_dict = {
    "B2_FB": "red",
    "B2_FL": "blue",
    "B2_GEMX": "green",
    "B2_NextGEM": "orange",
    "B2_Parse": "purple",
    "B2_Scale": "yellow"
}

# Define line styles for each condition (F1A, F1B, etc.)
line_style_dict = {
    "F1A": "-",
    "F1B": "--",
    "F5A": "-.",
    "F5B": ":"
}

# Create a plot
plt.figure(figsize=(12, 8))

# Iterate over each file
for file in file_list:
  # Extract the base name (without extension) and split by underscores
  base_name = os.path.basename(file).replace('.xlsx', '')
  parts = base_name.split('_')
    
    # Extract prefix (e.g., "B2_FB") and condition (e.g., "F1A")
  prefix = "_".join(parts[:2])  # First two parts for prefix
  condition = parts[2]          # Third part for condition
    
    # Read the Excel file (specify engine='openpyxl' to handle .xlsx files)
  data = pd.read_excel(file, engine='openpyxl')
    
    # Plot the data, using color for prefix and line style for condition
  plt.plot(data["reads"], data['genes'], label=f"{prefix}_{condition}", color=color_dict.get(prefix, 'black'), linestyle=line_style_dict.get(condition, '-'))

# Set titles and labels
plt.title("Saturation Curves for All Samples")
plt.xlabel("Mean Reads per Cell")
plt.ylabel("Median Genes per Cell")
plt.xticks(ticks=[5000, 10000, 15000, 20000, 25000, 30000], labels=["5k", "10k", "15k", "20k", "25k", "30k"])
plt.ylim(0, 3000)  # Adjust as needed



# Add a legend outside the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Save the figure with DPI=200 before showing
plt.savefig("/fh/fast/_IRC/FHIL/user/dgratz/BM02/figures/saturation_curves_all_samples.png", dpi=200, bbox_inches='tight')

# Display the plot
plt.tight_layout()
plt.show()

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import os
from sklearn.preprocessing import StandardScaler

data_path = # Specify the path to paraquet files

gene_names = # INPUT GENE NAMES

selected_genes = # SELECT YOUR GENE OF INTEREST

tissues_of_interest = #select your tissue of interest

# Age group mapping for each tissue
group_labels_dict = {
    "Brain-Cortex": ["20-49"] * 9 + ["50-59"] * 18 + ["60-79"] * 55,
    # Add more tissues and their age groups if needed
}

# Initialize lists to store heatmap data
heatmap_data_list = []
age_groups_list = []

# Loop through each tissue file
for file_name in os.listdir(data_path):
    if file_name.endswith(".parquet"): 
        tissue_name = file_name.replace(".parquet", "")
        
        print(f"Found tissue file: {tissue_name}")
        
        if tissue_name not in tissues_of_interest:
            print(f"Skipping tissue: {tissue_name}")
            continue
        
        tissue_data = pd.read_parquet(os.path.join(data_path, file_name))
        
        if tissue_data.shape[0] != len(gene_names):
            raise ValueError(f"Error: The number of gene names does not match the number of rows in {tissue_name}")
        
        tissue_data.index = gene_names
        
        # Subset the data for the selected genes
        tissue_data = tissue_data.loc[selected_genes, :]
        
        # Log-transform
        tissue_data_log10 = np.log10(tissue_data + 1e-9)  # Add a small constant to avoid log(0)
        
        # Z-score normalization
        scaler = StandardScaler()
        tissue_data_zscore = pd.DataFrame(scaler.fit_transform(tissue_data_log10.T).T,
                                          index=tissue_data.index, columns=tissue_data.columns)
        
        # Check for any NA values and remove rows/columns that are fully NA
        tissue_data_zscore = tissue_data_zscore.dropna()
        
        # Add tissue labels and store heatmap data
        heatmap_data_list.append(tissue_data_zscore)
        age_groups_list.extend(group_labels_dict[tissue_name])

# Combine all tissue data into one DataFrame
combined_data = pd.concat(heatmap_data_list, axis=1)

# Transpose the combined data for vertical heatmap
combined_data_T = combined_data.T

# Create a DataFrame for annotations
annotation_col = pd.DataFrame({'AgeGroup': pd.Categorical(age_groups_list, categories=['20-49', '50-59', '60-79'])},
                               index=combined_data_T.index)

# Reorder combined_data_T DataFrame by age group
ordered_age_groups = ['20-49', '50-59', '60-79']
ordered_indices = annotation_col['AgeGroup'].map(lambda x: ordered_age_groups.index(x))
combined_data_T = combined_data_T.iloc[np.argsort(ordered_indices)]

annotation_colors = {'20-49': '#ff9999', '50-59': '#66b3ff', '60-79': '#99ff99'}
num_colors = [list(annotation_colors.keys()).index(age) for age in annotation_col['AgeGroup']]
annotation_image = np.array(num_colors).reshape(-1, 1)

# Set up the figure and axes using gridspec
fig = plt.figure(figsize=(8, 4))  # Adjust figure size to fit both plots
gs = fig.add_gridspec(1, 2, width_ratios=[0.1, 4], wspace=0.01)  # Further reduce width ratio for annotation bar

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

# Plot the annotation data
cmap = mcolors.ListedColormap(list(annotation_colors.values()))
bounds = list(range(len(annotation_colors) + 1))
norm = mcolors.BoundaryNorm(bounds, cmap.N)

# Set the aspect ratio to ensure the annotation bar is thin
ax1.imshow(annotation_image, cmap=cmap, norm=norm, aspect='auto')
ax1.set_xticks([])  # Remove x-ticks
ax1.set_yticks([])  # Remove y-ticks
ax1.set_title('Age Groups', pad=20)

# Create the heatmap with square cells
num_rows, num_cols = combined_data_T.shape
aspect_ratio = num_cols / num_rows  # Aspect ratio to make cells square

# Adjust the figure size dynamically to maintain square cells
fig.set_size_inches(num_cols * 0.3, num_rows * 0.3)  # Adjust 0.3 as necessary

sns.heatmap(
    combined_data_T,
    cmap="coolwarm",
    annot=True,  
    fmt=".2f",  
    linewidths=0.2,
    cbar_kws={"label": "Expression"},  
    xticklabels=True,  
    yticklabels=False,  
    ax=ax2
)

# Add a custom color bar for age groups as a legend
patches_list = [patches.Patch(color=color, label=age) for age, color in annotation_colors.items()]
ax2.legend(handles=patches_list, title="Age Groups", bbox_to_anchor=(1.2, 1), loc='upper left')

ax2.set_xlabel('Genes', labelpad=20)
ax2.xaxis.set_label_position('top')  

ax2.xaxis.set_ticks_position('none')  # Hide ticks
ax2.xaxis.set_tick_params(labeltop=True, labelbottom=False)  

# Save the heatmap as PDF and PNG with larger dimensions
result_directory = "results_brain_heatmap"
os.makedirs(result_directory, exist_ok=True)

plt.savefig(os.path.join(result_directory, 'heatmap_apoe_sod1_age_groups.pdf'), format='pdf', bbox_inches='tight')
plt.savefig(os.path.join(result_directory, 'heatmap_apoe_sod1_age_groups.png'), format='png', bbox_inches='tight', dpi=300)

plt.show()

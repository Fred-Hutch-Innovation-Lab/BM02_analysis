#!/usr/bin/env python3

import scanpy as sc
import sc3s
import pandas as pd
import numpy as np
from pathlib import Path
from pyprojroot import here

def process_h5ad_file(h5ad_path):
    """Process a single h5ad file with SC3s clustering."""
    print(f"Processing {h5ad_path}")
    
    adata = sc.read_h5ad(h5ad_path)
    adata.obsm["X_pca"] = adata.obsm["X_PCA_allgenes"]
    
    sc3s.tl.consensus(adata, n_clusters=list(range(5, 30, 1)))
    
    # Create DataFrame with cell names and cluster assignments
    results_df = pd.DataFrame(index=adata.obs.index)
    
    # Add cluster assignments for each k
    for k in range(5, 30):
        col_name = f'sc3s_{k}'
        if col_name in adata.obs.columns:
            # Convert categorical to numeric
            results_df[f'k{k}'] = adata.obs[col_name].cat.codes
    
    return results_df

def main():
    # Get the directory containing h5ad files
    h5ad_dir = here('rds/3p/sc3_anndata')
    output_dir = here('rds/3p/sc3s_results')
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process all h5ad files
    for h5ad_file in Path(h5ad_dir).glob('*.h5ad'):
        kit_name = h5ad_file.stem
        print(f"Processing kit: {kit_name}")
        
        try:
            results_df = process_h5ad_file(h5ad_file)
            
            # Save results for this kit
            output_file = output_dir / f"{kit_name}_sc3s_results.csv"
            results_df.to_csv(output_file)
            
            print(f"Saved results for {kit_name} to {output_file}")
            
        except Exception as e:
            print(f"Error processing {kit_name}: {str(e)}")
            continue

if __name__ == "__main__":
    main() 
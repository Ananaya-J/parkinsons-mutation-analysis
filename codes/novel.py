import pandas as pd
import numpy as np

def clean_af_column(series):
    """Helper function to clean AF columns"""
    return pd.to_numeric(series.replace('.', '1.0'), errors='coerce').fillna(1.0)

def identify_novel_mutations(vcf_file, output_file="novel_mutations1.csv", population_freq_threshold=0.01):
    """
    Identify potentially novel mutations from an annotated VCF file
    
    Parameters:
    vcf_file (str): Path to annotated VCF file (e.g., from ANNOVAR)
    output_file (str): Name of output CSV file
    population_freq_threshold (float): Maximum allele frequency to consider novel
    
    Returns:
    DataFrame: Filtered novel mutations
    """
    # Read the annotated VCF file
    data = pd.read_csv(vcf_file, sep="\t")
    
    print("Available columns:", list(data.columns))
    
    # Clean AF columns
    af_columns = [
        'AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw',
        'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe',
        'AF_fin', 'AF_asj', 'AF_oth', 'non_topmed_AF_popmax',
        'non_neuro_AF_popmax', 'non_cancer_AF_popmax', 'controls_AF_popmax'
    ]
    
    # Convert all AF columns to numeric, replacing '.' with 1.0 (not novel)
    for col in af_columns:
        data[col] = clean_af_column(data[col])
    
    # Create mask for novel variants
    novel_mask = (
        # Check main AF columns
        (data['AF'] < population_freq_threshold) &
        (data['AF_popmax'] < population_freq_threshold) &
        
        # Filter for functional impact
        data['ExonicFunc.refGene'].isin([
            'nonsynonymous_SNV',
            'frameshift_deletion',
            'frameshift_insertion',
            'stopgain',
            'stoploss'
        ])
    )
    
    # Add population-specific AF filters
    population_af_cols = ['AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 'AF_asj', 'AF_oth']
    for col in population_af_cols:
        novel_mask &= (data[col] < population_freq_threshold)
    
    # Extract novel mutations
    novel_mutations = data[novel_mask].copy()
    
    # Sort by AF (lower is more novel)
    novel_mutations = novel_mutations.sort_values(
        by=['AF_popmax', 'AF'],
        ascending=[True, True]
    )
    
    # Save results
    try:
        novel_mutations.to_csv(output_file, index=False)
        print(f"\nResults saved to: {output_file}")
    except Exception as e:
        print(f"Error saving file: {str(e)}")
        
    return novel_mutations

# Usage example
if __name__ == "__main__":
    try:
        # Read and process the file
        novel_mutations = identify_novel_mutations("annotated.hg19_multianno.txt")
        
        # Print summary
        print(f"\nTotal novel mutations identified: {len(novel_mutations)}")
        
        # Print top mutations
        important_cols = [
            'Chr', 'Start', 'End', 'Ref', 'Alt', 
            'Gene.refGene', 'ExonicFunc.refGene', 
            'AF', 'AF_popmax'
        ]
        
        print("\nTop 10 most novel mutations:")
        print(novel_mutations[important_cols].head(10))
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        print("Please check if your input file exists and has the correct format.")

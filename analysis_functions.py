import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import numpy as np

def analyze_gene_exon_distances(df):
    """Analyze the distances between genes and their first exons"""
    results = []
    
    # Get all genes
    genes = df[df['type'] == 'gene']
    
    for _, gene in genes.iterrows():
        gene_id = gene['ID']
        
        # Find mRNAs for this gene
        mrnas = df[(df['type'] == 'mRNA') & (df['Parent'] == gene_id)]
        
        for _, mrna in mrnas.iterrows():
            mrna_id = mrna['ID']
            
            # Find exons for this mRNA
            exons = df[(df['type'] == 'exon') & (df['Parent'] == mrna_id)]
            
            if not exons.empty:
                # Sort exons by position
                if gene['strand'] == '+':
                    first_exon = exons.sort_values('start').iloc[0]
                else:
                    first_exon = exons.sort_values('end', ascending=False).iloc[0]
                
                # Calculate distance
                if gene['strand'] == '+':
                    distance = first_exon['start'] - gene['start']
                else:
                    distance = gene['end'] - first_exon['end']
                
                results.append({
                    'gene_id': gene_id,
                    'mrna_id': mrna_id,
                    'gene_start': gene['start'],
                    'gene_end': gene['end'],
                    'first_exon_start': first_exon['start'],
                    'first_exon_end': first_exon['end'],
                    'distance': distance
                })
    
    return pd.DataFrame(results)

def analyze_intron_sizes(df):
    """Calculate intron sizes for all genes"""
    results = []
    
    # Get all mRNAs
    mrnas = df[df['type'] == 'mRNA']
    
    for _, mrna in mrnas.iterrows():
        mrna_id = mrna['ID']
        
        # Find exons for this mRNA
        exons = df[(df['type'] == 'exon') & (df['Parent'] == mrna_id)]
        
        if len(exons) > 1:
            # Sort exons by position
            exons = exons.sort_values('start')
            
            # Calculate introns
            for i in range(len(exons) - 1):
                current_exon = exons.iloc[i]
                next_exon = exons.iloc[i + 1]
                
                intron_start = current_exon['end'] + 1
                intron_end = next_exon['start'] - 1
                intron_size = intron_end - intron_start + 1
                
                if intron_size > 0:  # Ensure it's a valid intron
                    results.append({
                        'mrna_id': mrna_id,
                        'intron_number': i + 1,
                        'intron_start': intron_start,
                        'intron_end': intron_end,
                        'intron_size': intron_size
                    })
    
    return pd.DataFrame(results)

def add_predefined_analyses(df):
    """Add predefined analyses to the Streamlit app"""
    st.header("Pre-defined Analyses")
    
    # Create a dropdown for different analyses
    analysis_type = st.selectbox(
        "Select analysis type (select, let load, press run):",
        ["Gene-Exon Distances", "Intron Sizes", "CDS vs. Gene Length Ratio"]
    )
    
    if analysis_type == "Gene-Exon Distances":
        if st.button("Run Gene-Exon Distance Analysis"):
            with st.spinner("Analyzing gene-exon distances..."):
                results = analyze_gene_exon_distances(df)
                
                if not results.empty:
                    st.subheader("Gene to First Exon Distances")
                    st.dataframe(results)
                    
                    # Create a histogram of distances
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ax.hist(results['distance'], bins=30)
                    ax.set_xlabel("Distance (bp)")
                    ax.set_ylabel("Count")
                    ax.set_title("Distribution of Gene to First Exon Distances")
                    st.pyplot(fig)
                    
                    # Add download button
                    csv = results.to_csv(index=False)
                    st.download_button(
                        "Download results as CSV",
                        csv,
                        "gene_exon_distances.csv",
                        "text/csv"
                    )
                else:
                    st.info("No gene-exon relationships found in the data")
    
    elif analysis_type == "Intron Sizes":
        if st.button("Run Intron Size Analysis"):
            with st.spinner("Analyzing intron sizes..."):
                results = analyze_intron_sizes(df)
                
                if not results.empty:
                    st.subheader("Intron Sizes")
                    st.dataframe(results)
                    
                    # Summary statistics
                    st.write(f"Total introns: {len(results)}")
                    st.write(f"Average intron size: {results['intron_size'].mean():.2f} bp")
                    st.write(f"Median intron size: {results['intron_size'].median():.2f} bp")
                    st.write(f"Largest intron: {results['intron_size'].max()} bp")
                    
                    # Create a histogram
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ax.hist(results['intron_size'], bins=30)
                    ax.set_xlabel("Intron Size (bp)")
                    ax.set_ylabel("Count")
                    ax.set_title("Distribution of Intron Sizes")
                    st.pyplot(fig)
                    
                    # Add download button
                    csv = results.to_csv(index=False)
                    st.download_button(
                        "Download results as CSV",
                        csv,
                        "intron_sizes.csv",
                        "text/csv"
                    )
                else:
                    st.info("No introns found in the data")
    
    elif analysis_type == "CDS vs. Gene Length Ratio":
        if st.button("Run CDS/Gene Length Analysis"):
            with st.spinner("Analyzing CDS vs gene lengths..."):
                # Get all genes
                genes = df[df['type'] == 'gene']
                
                results = []
                
                for _, gene in genes.iterrows():
                    gene_id = gene['ID']
                    gene_length = gene['end'] - gene['start'] + 1
                    
                    # Find all CDS regions for this gene
                    cds_regions = df[(df['type'] == 'CDS') & 
                                    (df['attributes'].str.contains(gene_id))]
                    
                    if not cds_regions.empty:
                        total_cds_length = sum(cds_regions['end'] - cds_regions['start'] + 1)
                        
                        # Calculate ratio
                        ratio = (total_cds_length / gene_length) * 100
                        
                        results.append({
                            'gene_id': gene_id,
                            'gene_length': gene_length,
                            'total_cds_length': total_cds_length,
                            'cds_percentage': ratio
                        })
                
                if results:
                    result_df = pd.DataFrame(results)
                    
                    st.subheader("CDS vs. Gene Length Analysis")
                    st.dataframe(result_df)
                    
                    # Summary statistics
                    st.write(f"Average CDS percentage: {result_df['cds_percentage'].mean():.2f}%")
                    
                    # Create a histogram
                    fig, ax = plt.subplots(figsize=(10, 6))
                    ax.hist(result_df['cds_percentage'], bins=20)
                    ax.set_xlabel("CDS Percentage of Gene Length")
                    ax.set_ylabel("Count")
                    ax.set_title("Distribution of CDS Percentage")
                    st.pyplot(fig)
                    
                    # Add download button
                    csv = result_df.to_csv(index=False)
                    st.download_button(
                        "Download results as CSV",
                        csv,
                        "cds_gene_ratio.csv",
                        "text/csv"
                    )
                else:
                    st.info("No gene-CDS relationships found in the data")


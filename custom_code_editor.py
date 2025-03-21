# custom_code_editor.py
# Add this as a separate file to keep your main app clean

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import io
from contextlib import redirect_stdout

def add_custom_code_editor(active_df, comparison_df=None, comparison_mode=False):
    """
    Adds a custom code editor to the Streamlit app
    
    Parameters:
    - active_df: The currently active dataframe
    - comparison_df: The comparison dataframe (if in comparison mode)
    - comparison_mode: Whether comparison mode is enabled
    """
    st.header("Custom Code Analysis")
    st.write("Write and run custom Python code to analyze the genome data")
    
    # Show available dataframe names
    available_dfs = {
        "df": "Current active dataframe",
    }
    
    if comparison_mode and comparison_df is not None:
        available_dfs["comparison_df"] = "Comparison dataframe"
    
    st.info("Available dataframes: " + ", ".join([f"`{k}` ({v})" for k, v in available_dfs.items()]))
    
    # Add documentation for common operations
    with st.expander("Common GFF3 Analysis Examples"):
        st.markdown("""
        ### Examples of useful GFF3 analyses:
        
        **1. Count features by type:**
        ```python
        feature_counts = df['type'].value_counts()
        print(feature_counts)
        ```
        
        **2. Calculate gene lengths:**
        ```python
        gene_df = df[df['type'] == 'gene']
        gene_df['length'] = gene_df['end'] - gene_df['start'] + 1
        print(f"Average gene length: {gene_df['length'].mean():.2f} bp")
        print(f"Median gene length: {gene_df['length'].median():.2f} bp")
        print(f"Longest gene: {gene_df['length'].max()} bp")
        ```
        
        **3. Find distance between genes and first exons:**
        ```python
        results = []
        
        # Get sample gene IDs (first 5)
        gene_ids = df[df['type'] == 'gene']['ID'].unique()[:5]
        
        for gene_id in gene_ids:
            # Get gene info
            gene = df[df['ID'] == gene_id].iloc[0]
            
            # Find mRNAs for this gene
            mrnas = df[df['Parent'] == gene_id]
            
            for _, mrna in mrnas.iterrows():
                mrna_id = mrna['ID']
                
                # Find exons for this mRNA
                exons = df[(df['type'] == 'exon') & (df['Parent'] == mrna_id)]
                
                if not exons.empty:
                    # Sort exons by position
                    if gene['strand'] == '+':
                        exons = exons.sort_values('start')
                    else:
                        exons = exons.sort_values('end', ascending=False)
                    
                    first_exon = exons.iloc[0]
                    
                    # Calculate distance
                    if gene['strand'] == '+':
                        distance = first_exon['start'] - gene['start']
                    else:
                        distance = gene['end'] - first_exon['end']
                    
                    results.append({
                        'gene_id': gene_id,
                        'mrna_id': mrna_id,
                        'distance': distance
                    })
        
        result_df = pd.DataFrame(results)
        print(result_df)
        ```
        
        **4. Compare gene counts between files (in comparison mode):**
        ```python
        if 'comparison_df' in locals():
            # Count genes in both files
            gene_count = len(df[df['type'] == 'gene'])
            comp_gene_count = len(comparison_df[comparison_df['type'] == 'gene'])
            
            print(f"Current file: {gene_count} genes")
            print(f"Comparison file: {comp_gene_count} genes")
            print(f"Difference: {gene_count - comp_gene_count} genes")
        ```
        
        **5. Plot exon length distribution:**
        ```python
        exons = df[df['type'] == 'exon']
        exons['length'] = exons['end'] - exons['start'] + 1
        
        # Create a histogram
        fig, ax = plt.subplots()
        ax.hist(exons['length'], bins=50)
        ax.set_xlabel('Exon Length (bp)')
        ax.set_ylabel('Count')
        ax.set_title('Distribution of Exon Lengths')
        
        # Display the plot
        st.pyplot(fig)
        ```
        """)
    
    # Code editor
    custom_code = st.text_area(
        "Enter your Python code below:",
        height=300,
        help="Use 'df' to refer to the active dataframe. If in comparison mode, use 'comparison_df' for the comparison dataframe."
    )
    
    # Create a two-column layout for output options
    col1, col2 = st.columns(2)
    
    with col1:
        output_option = st.radio(
            "Output type:",
            ["Text output", "DataFrame", "Plot"],
            help="Select how you want to display the results"
        )
    
    with col2:
        if output_option == "DataFrame":
            max_rows = st.number_input("Max rows to display:", min_value=1, value=10)
        elif output_option == "Plot":
            st.info("Use matplotlib to create plots and they'll be displayed automatically")
    
    # Run button
    if st.button("Run Code"):
        if not custom_code.strip():
            st.warning("Please enter some code to run")
            return
        
        try:
            # Create a namespace with the active dataframe(s)
            local_namespace = {"df": active_df, "pd": pd, "np": np, "plt": plt, "st": st}
            
            if comparison_mode and comparison_df is not None:
                local_namespace["comparison_df"] = comparison_df
            
            # Create a buffer for capturing print outputs
            stdout_buffer = io.StringIO()
            
            # Execute the code and capture output
            with redirect_stdout(stdout_buffer):
                exec(custom_code, globals(), local_namespace)
            
            stdout_output = stdout_buffer.getvalue()
            
            # Display results based on selected output type
            if output_option == "Text output":
                if stdout_output:
                    st.text_area("Output:", value=stdout_output, height=250, disabled=True)
                else:
                    st.info("No text output produced")
                    
            elif output_option == "DataFrame":
                # Check if result_df was created
                if "result_df" in local_namespace:
                    st.subheader("Result DataFrame:")
                    st.dataframe(local_namespace["result_df"].head(max_rows))
                    
                    # Add download button
                    csv = local_namespace["result_df"].to_csv(index=False)
                    st.download_button(
                        "Download results as CSV",
                        csv,
                        "analysis_results.csv",
                        "text/csv",
                        key='download-csv'
                    )
                else:
                    st.warning("No 'result_df' variable found in the code. Create a DataFrame called 'result_df' to display it here.")
                    
                    # Still show any text output
                    if stdout_output:
                        st.text_area("Text output:", value=stdout_output, height=150, disabled=True)
            
            elif output_option == "Plot":
                # Check if plt was used (as a very simple heuristic)
                if "plt" in custom_code:
                    # The plot should have been displayed using st.pyplot(fig)
                    if "st.pyplot" not in custom_code:
                        st.warning("No plot displayed. Make sure to use 'st.pyplot(fig)' in your code to display plots.")
                else:
                    st.warning("No matplotlib plotting code detected. Use plt.subplots() to create figures.")
                
                # Still show any text output
                if stdout_output:
                    st.text_area("Text output:", value=stdout_output, height=150, disabled=True)
            
            st.success("Code executed successfully")
            
        except Exception as e:
            st.error(f"Error executing code: {str(e)}")
            # Show traceback for debugging
            import traceback
            st.text_area("Error traceback:", value=traceback.format_exc(), height=250, disabled=True)


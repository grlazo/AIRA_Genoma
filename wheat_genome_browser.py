import streamlit as st
import pandas as pd
import re
import tempfile
import os
import requests
from langchain_ollama import OllamaLLM
from langchain_experimental.agents import create_pandas_dataframe_agent

from custom_code_editor import add_custom_code_editor
from analysis_functions import add_predefined_analyses

# Set page config must be the first Streamlit command
st.set_page_config(page_title="Genome Feature Finder", layout="wide")

# GO Term Handler class for Gene Ontology functionality
class GOTermHandler:
    """
    Handler for Gene Ontology (GO) term lookups and searches
    """
    def __init__(self):
        self.go_term_cache = {}  # Cache for GO term lookups
        self.go_semantic_cache = {}  # Cache for semantic term to GO ID mappings
        # Basic mapping for common aluminum/metal-related GO terms
        self.metal_go_terms = {
            "GO:0046872": "metal ion binding",
            "GO:0046914": "transition metal ion binding",
            "GO:0030145": "manganese ion binding",
            "GO:0005507": "copper ion binding",
            "GO:0008270": "zinc ion binding",
            "GO:0005506": "iron ion binding",
            "GO:0043167": "ion binding",
            "GO:0046873": "metal ion transmembrane transporter activity",
            "GO:0015082": "di-, tri-valent inorganic cation transmembrane transporter activity",
            "GO:0015093": "ferrous iron transmembrane transporter activity",
            "GO:0015680": "intracellular aluminum ion transport",
            "GO:0015682": "ferric iron transport",
            "GO:0015683": "ferrous iron transport",
            "GO:0006829": "zinc II ion transport",
            "GO:0006826": "iron ion transport",
            "GO:0006824": "cobalt ion transport",
            "GO:0006825": "copper ion transport",
            "GO:0006823": "aluminum ion transport",
            "GO:0046688": "response to aluminum ion",
            "GO:0010044": "response to aluminum stress",
        }
    
    def search_go_terms_in_ontology_db(self, search_term):
        """
        Search for GO terms in public databases
        
        This function queries the QuickGO API to search for GO terms
        related to the given search term
        """
        try:
            # Use QuickGO API to search for GO terms
            # Documentation: https://www.ebi.ac.uk/QuickGO/api/index.html
            url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query={search_term}"
            response = requests.get(url, headers={"Accept": "application/json"})
            if response.status_code == 200:
                data = response.json()
                results = []
                if 'results' in data and len(data['results']) > 0:
                    for result in data['results']:
                        go_id = result.get('id')
                        go_name = result.get('name')
                        go_aspect = result.get('aspect', '')  # P (process), F (function), C (component)
                        if go_id and go_name:
                            results.append({
                                'go_id': go_id,
                                'go_name': go_name,
                                'go_aspect': go_aspect
                            })
                            # Update cache
                            self.go_term_cache[go_id] = go_name
                            self.go_semantic_cache[go_name.lower()] = go_id
                return results
            else:
                return []
        except Exception as e:
            st.warning(f"Error searching GO terms: {str(e)}")
            return []
    
    def get_go_term_info(self, go_id):
        """Get information about a specific GO term"""
        if go_id in self.go_term_cache:
            return self.go_term_cache[go_id]
        
        try:
            # Use QuickGO API to get GO term information
            url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}"
            response = requests.get(url, headers={"Accept": "application/json"})
            if response.status_code == 200:
                data = response.json()
                if 'results' in data and len(data['results']) > 0:
                    go_name = data['results'][0].get('name', '')
                    self.go_term_cache[go_id] = go_name
                    return go_name
            return go_id  # Return the ID if name can't be found
        except Exception as e:
            st.warning(f"Error fetching GO term info: {str(e)}")
            return go_id
    
    def search_semantic_matches(self, search_term, df):
        """Search for genes with GO terms semantically related to the search term"""
        matching_go_ids = []
        
        # Check if the search term is a GO ID
        if search_term.upper().startswith("GO:"):
            matching_go_ids.append(search_term.upper())
        else:
            # Try to find matching GO terms
            results = self.search_go_terms_in_ontology_db(search_term)
            matching_go_ids = [result['go_id'] for result in results]
            
            # Also check our predefined list for partial matches
            for go_id, go_name in self.metal_go_terms.items():
                if search_term.lower() in go_name.lower():
                    matching_go_ids.append(go_id)
        
        # Look for matches in the dataframe
        found_rows = []
        
        # Check Ontology_term column if it exists
        if 'Ontology_term' in df.columns:
            for go_id in matching_go_ids:
                matches = df[df['Ontology_term'].str.contains(go_id, na=False)]
                if not matches.empty:
                    found_rows.append(matches)
        
        # Also check attributes column for GO terms
        if 'attributes' in df.columns:
            for go_id in matching_go_ids:
                matches = df[df['attributes'].str.contains(go_id, na=False)]
                if not matches.empty:
                    found_rows.append(matches)
        
        # Also check Note column for semantic matches
        if 'Note' in df.columns:
            matches = df[df['Note'].str.contains(search_term, case=False, na=False)]
            if not matches.empty:
                found_rows.append(matches)
        
        # Combine all matches
        if found_rows:
            result_df = pd.concat(found_rows).drop_duplicates()
            return result_df
        else:
            return pd.DataFrame()

def parse_gff3_file(file_path):
    """
    Parse a GFF3 file into a pandas DataFrame with proper column names
    and extract attributes into separate columns.
    
    :param file_path: Path to GFF3 file
    :return: DataFrame with parsed GFF3 data
    """
    # Standard GFF3 columns
    columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    
    # Read the file, skip comment lines
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) == 9:
                data.append(fields)
    
    # Create DataFrame with standard columns
    df = pd.DataFrame(data, columns=columns)
    
    # Convert numeric columns
    df['start'] = pd.to_numeric(df['start'])
    df['end'] = pd.to_numeric(df['end'])
    df['score'] = df['score'].replace('.', None)
    df['score'] = pd.to_numeric(df['score'], errors='coerce')
    
    # Parse attributes into a dictionary
    def parse_attributes(attr_string):
        attrs = {}
        for pair in attr_string.split(';'):
            if '=' in pair:
                key, value = pair.split('=', 1)
                attrs[key.strip()] = value.strip()
        return attrs
    
    # Extract attributes into a new column as dictionaries
    df['attrs_dict'] = df['attributes'].apply(parse_attributes)
    
    # Extract GO terms from attributes
    df = extract_go_terms_from_attributes(df)
    
    # Find common attributes (those that appear in at least 10% of rows)
    total_rows = len(df)
    attr_keys = {}
    for attrs in df['attrs_dict']:
        for key in attrs:
            attr_keys[key] = attr_keys.get(key, 0) + 1
    
    common_attrs = [key for key, count in attr_keys.items() 
                   if count >= total_rows * 0.1]
    
    # Add common attributes as separate columns
    for attr in common_attrs:
        df[attr] = df['attrs_dict'].apply(lambda x: x.get(attr, ''))
    
    return df

# Function to extract GO terms from GFF3 attributes
def extract_go_terms_from_attributes(df):
    """Extract GO terms from the attributes column and create a separate column"""
    if 'attributes' not in df.columns:
        return df
    
    # Extract GO terms from attributes using regex
    go_terms = []
    for attr in df['attributes']:
        terms = []
        # Look for Ontology_term=GO:XXXX
        if 'Ontology_term=' in attr:
            # Extract all GO terms
            terms_part = attr.split('Ontology_term=')[1].split(';')[0]
            # Split multiple GO terms if they exist
            terms.extend(terms_part.split(','))
        go_terms.append(','.join(terms) if terms else '')
    
    # Add as a new column
    df['Ontology_term'] = go_terms
    return df

def get_feature_statistics(df):
    """Get basic statistics about the GFF3 features"""
    stats = {}
    
    # Count of feature types
    stats['feature_counts'] = df['type'].value_counts()
    
    # Count of genes
    stats['gene_count'] = len(df[df['type'] == 'gene'])
    
    # Average gene length
    gene_df = df[df['type'] == 'gene']
    stats['avg_gene_length'] = (gene_df['end'] - gene_df['start']).mean()
    
    # Feature density (features per Mb)
    max_position = df['end'].max()
    stats['feature_density'] = len(df) / (max_position / 1_000_000)
    
    # Count features by chromosome/sequence
    stats['seq_counts'] = df['seqid'].value_counts()
    
    return stats

def run_genome_app():
    """Main Streamlit app function"""
    st.title("Genome Feature Finder")
    st.markdown("### Query and explore GFF3 genome annotations")
    
    # Sidebar for file upload
    st.sidebar.title("Data Import")
    
    # Multi-file uploader
    uploaded_files = st.sidebar.file_uploader(
        "Upload GFF3 files", 
        type=["gff", "gff3"],
        accept_multiple_files=True,
        help="Upload one or more GFF3 files containing genome annotations."
    )
    
    # Store dataframes in session state
    if 'dataframes' not in st.session_state:
        st.session_state.dataframes = {}
        st.session_state.active_df = None
        st.session_state.tmp_files = []
    
    # Process uploaded files
    if uploaded_files:
        for uploaded_file in uploaded_files:
            file_name = uploaded_file.name
            
            # Skip if already processed
            if file_name in st.session_state.dataframes:
                continue
                
            # Create a temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix=".gff3") as tmp:
                tmp.write(uploaded_file.getvalue())
                tmp_path = tmp.name
                st.session_state.tmp_files.append(tmp_path)
            
            # Process the file
            with st.spinner(f"Processing {file_name}..."):
                try:
                    df = parse_gff3_file(tmp_path)
                    st.session_state.dataframes[file_name] = df
                    st.success(f"Successfully processed {file_name} - {len(df)} features loaded")
                except Exception as e:
                    st.error(f"Error processing {file_name}: {str(e)}")
    
    # File selection dropdown (only show if files are available)
    if st.session_state.dataframes:
        st.sidebar.title("File Selection")
        file_selector = st.sidebar.selectbox(
            "Select file to analyze",
            list(st.session_state.dataframes.keys())
        )
        
        # Comparison mode (if multiple files uploaded)
        if len(st.session_state.dataframes) > 1:
            comparison_mode = st.sidebar.checkbox("Enable comparison mode")
            if comparison_mode:
                comparison_file = st.sidebar.selectbox(
                    "Select file to compare with",
                    [f for f in st.session_state.dataframes.keys() if f != file_selector]
                )
        else:
            comparison_mode = False
        
        # Set active dataframe
        active_df = st.session_state.dataframes[file_selector]
        st.session_state.active_df = active_df
        
        # If comparison mode is enabled
        if comparison_mode:
            comparison_df = st.session_state.dataframes[comparison_file]
            st.sidebar.write(f"Comparing: {file_selector} with {comparison_file}")
        
        # Main UI with tabs for different functions
        main_tabs = st.tabs(["Genome Overview", "Data Explorer", "Gene Ontology", "Query", "Analyses", "Code Editor"])
        
        # Genome Overview Tab
        with main_tabs[0]:
            st.header(f"Genome Overview: {file_selector}")
            
            # Create columns for statistics
            col1, col2 = st.columns(2)
            
            # Get statistics
            stats = get_feature_statistics(active_df)
            
            # Display statistics in columns
            with col1:
                st.metric("Total Features", len(active_df))
                st.metric("Total Genes", stats['gene_count'])
                if 'mRNA' in active_df['type'].values:
                    st.metric("Total mRNAs", len(active_df[active_df['type'] == 'mRNA']))
                
            with col2:
                if stats['avg_gene_length']:
                    st.metric("Average Gene Length", f"{stats['avg_gene_length']:.0f} bp")
                st.metric("Feature Density", f"{stats['feature_density']:.2f} features/Mb")
            
            # Feature type distribution
            st.subheader("Feature Distribution")
            st.bar_chart(stats['feature_counts'])
            
            # Comparison stats if enabled
            if comparison_mode:
                st.header(f"Comparison with: {comparison_file}")
                comp_stats = get_feature_statistics(comparison_df)
                
                # Create comparison metrics
                col1, col2 = st.columns(2)
                with col1:
                    st.metric("Features Difference", 
                             len(active_df) - len(comparison_df),
                             delta_color="normal")
                    st.metric("Genes Difference", 
                             stats['gene_count'] - comp_stats['gene_count'], 
                             delta_color="normal")
                
                with col2:
                    if stats['avg_gene_length'] and comp_stats['avg_gene_length']:
                        st.metric("Avg Gene Length Difference", 
                                 f"{stats['avg_gene_length'] - comp_stats['avg_gene_length']:.0f} bp", 
                                 delta_color="normal")
                    st.metric("Feature Density Difference", 
                             f"{stats['feature_density'] - comp_stats['feature_density']:.2f}", 
                             delta_color="normal")
                
                # Feature comparison
                st.subheader("Feature Type Comparison")
                # Combine feature counts
                combined_counts = pd.DataFrame({
                    file_selector: stats['feature_counts'],
                    comparison_file: comp_stats['feature_counts']
                }).fillna(0)
                
                st.bar_chart(combined_counts)
        
        # Data Explorer Tab
        with main_tabs[1]:
            # Tabs for different exploration views
            explorer_tabs = st.tabs(["Data Table", "Search & Filter", "Gene Details"])
            
            with explorer_tabs[0]:
                # Sample of the data
                st.subheader("Sample Data")
                st.dataframe(active_df.head(10))
                
                # Column information
                st.subheader("Column Information")
                col_info = {
                    'seqid': 'Chromosome or scaffold identifier',
                    'source': 'Source of the annotation',
                    'type': 'Feature type (gene, mRNA, exon, etc...)',
                    'start': 'Start position of the feature',
                    'end': 'End position of the feature',
                    'score': 'Score (often not used, shown as ".")',
                    'strand': 'Strand (+ or -)',
                    'phase': 'Phase (0, 1, 2, or ".")',
                    'attributes': 'Key-value pairs containing feature attributes'
                }
                for col, desc in col_info.items():
                    st.write(f"**{col}**: {desc}")
            
            with explorer_tabs[1]:
                # Search and filter options
                st.subheader("Search & Filter")
                
                # Feature type filter
                feature_types = ['All Types'] + sorted(active_df['type'].unique().tolist())
                selected_type = st.selectbox("Feature Type", feature_types)
                
                # Gene ID search
                if 'ID' in active_df.columns:
                    gene_search = st.text_input("Search by ID", 
                                               placeholder="Enter gene ID (e.g., Tm.TA299.r1.1AG0001000)")
                
                # Region filter
                col1, col2 = st.columns(2)
                with col1:
                    region_start = st.number_input("Region Start", min_value=1, value=1)
                with col2:
                    region_end = st.number_input("Region End", min_value=1, value=int(active_df['end'].max()))
                
                # Apply filters
                filtered_df = active_df.copy()
                
                if selected_type != 'All Types':
                    filtered_df = filtered_df[filtered_df['type'] == selected_type]
                
                if 'ID' in active_df.columns and gene_search:
                    filtered_df = filtered_df[filtered_df['ID'].str.contains(gene_search, case=False, na=False)]
                
                filtered_df = filtered_df[(filtered_df['start'] >= region_start) & 
                                         (filtered_df['end'] <= region_end)]
                
                # Show filtered results
                st.write(f"Showing {len(filtered_df)} features")
                st.dataframe(filtered_df)
            
            with explorer_tabs[2]:
                # Gene details view
                st.subheader("Gene Details")
                
                if 'ID' in active_df.columns:
                    # Get all gene IDs
                    gene_ids = sorted(active_df[active_df['type'] == 'gene']['ID'].dropna().unique().tolist())
                    
                    if gene_ids:
                        selected_gene = st.selectbox("Select Gene", gene_ids)
                        
                        # Get gene details
                        gene_row = active_df[active_df['ID'] == selected_gene].iloc[0]
                        
                        # Display gene details
                        st.write(f"**Gene ID:** {selected_gene}")
                        st.write(f"**Location:** {gene_row['seqid']}:{gene_row['start']}-{gene_row['end']} ({gene_row['strand']})")
                        st.write(f"**Length:** {gene_row['end'] - gene_row['start'] + 1} bp")
                        
                        if 'Note' in active_df.columns:
                            note = gene_row.get('Note', '')
                            if note:
                                st.write(f"**Description:** {note}")
                        
                        # Get all features related to this gene
                        if 'Parent' in active_df.columns:
                            # Get direct children (usually mRNAs)
                            children = active_df[active_df['Parent'] == selected_gene]
                            
                            # Get grand-children (exons, CDS, etc.)
                            child_ids = children['ID'].dropna().tolist()
                            grandchildren = pd.DataFrame()
                            
                            for child_id in child_ids:
                                grandchildren = pd.concat([grandchildren, active_df[active_df['Parent'] == child_id]])
                            
                            # Display structure
                            st.subheader("Gene Structure")
                            
                            if not children.empty:
                                st.write(f"**Transcripts:** {len(children)} mRNA(s)")
                                
                                # For each transcript, show exons, CDS, etc.
                                for i, transcript in children.iterrows():
                                    transcript_id = transcript['ID']
                                    transcript_features = active_df[active_df['Parent'] == transcript_id]
                                    
                                    exons = transcript_features[transcript_features['type'] == 'exon']
                                    cds = transcript_features[transcript_features['type'] == 'CDS']
                                    utrs = transcript_features[transcript_features['type'].str.contains('UTR', na=False)]
                                    
                                    st.write(f"**Transcript {transcript_id}:**")
                                    st.write(f"- Exons: {len(exons)}")
                                    st.write(f"- CDS regions: {len(cds)}")
                                    st.write(f"- UTRs: {len(utrs)}")
                                    
                                    # Display the features in a table
                                    feature_table = transcript_features[['type', 'start', 'end', 'strand']]
                                    feature_table['length'] = feature_table['end'] - feature_table['start'] + 1
                                    feature_table = feature_table.sort_values('start')
                                    
                                    st.dataframe(feature_table)
                                    
                            # Comparison if enabled
                            if comparison_mode and 'ID' in comparison_df.columns:
                                # Check if this gene exists in comparison file
                                comp_gene = comparison_df[comparison_df['ID'] == selected_gene]
                                if not comp_gene.empty:
                                    st.subheader(f"Gene Comparison: {comparison_file}")
                                    comp_gene_row = comp_gene.iloc[0]
                                    
                                    # Compare basic stats
                                    gene_length = gene_row['end'] - gene_row['start'] + 1
                                    comp_gene_length = comp_gene_row['end'] - comp_gene_row['start'] + 1
                                    
                                    st.write(f"**Location in {comparison_file}:** {comp_gene_row['seqid']}:{comp_gene_row['start']}-{comp_gene_row['end']} ({comp_gene_row['strand']})")
                                    st.write(f"**Length Difference:** {gene_length - comp_gene_length} bp")
                                    
                                    # Compare transcripts if available
                                    if 'Parent' in comparison_df.columns:
                                        comp_children = comparison_df[comparison_df['Parent'] == selected_gene]
                                        st.write(f"**Transcript count in {comparison_file}:** {len(comp_children)}")
                                        
                                        # Compare exon counts, etc.
                                        if not comp_children.empty and not children.empty:
                                            exon_count = sum(len(active_df[active_df['Parent'] == child_id]) for child_id in children['ID'])
                                            comp_exon_count = sum(len(comparison_df[comparison_df['Parent'] == child_id]) for child_id in comp_children['ID'])
                                            st.write(f"**Exon count difference:** {exon_count - comp_exon_count}")
                                else:
                                    st.info(f"This gene does not exist in {comparison_file}")
        
        # Gene Ontology Tab
        with main_tabs[2]:
            st.header("Gene Ontology Search")
            
            # Initialize GO handler
            go_handler = GOTermHandler()
            
            # Search options
            search_method = st.radio("Search method", ["GO ID", "Biological term"])
            
            if search_method == "GO ID":
                go_id = st.text_input("Enter GO ID (e.g., GO:0046872)")
                if go_id:
                    if not go_id.upper().startswith("GO:"):
                        st.warning("GO IDs should start with 'GO:'")
                    else:
                        go_term = go_handler.get_go_term_info(go_id.upper())
                        st.write(f"Term: {go_term}")
                        
                        if st.button("Search"):
                            results = go_handler.search_semantic_matches(go_id, active_df)
                            if not results.empty:
                                st.write(f"Found {len(results)} features with GO term {go_id}:")
                                relevant_cols = ['seqid', 'type', 'start', 'end', 'ID']
                                if 'Note' in results.columns:
                                    relevant_cols.append('Note')
                                if 'Ontology_term' in results.columns:
                                    relevant_cols.append('Ontology_term')
                                st.dataframe(results[relevant_cols])
                                
                                # Compare with other file if comparison mode enabled
                                if comparison_mode:
                                    comp_results = go_handler.search_semantic_matches(go_id, comparison_df)
                                    if not comp_results.empty:
                                        st.write(f"Found {len(comp_results)} features with GO term {go_id} in {comparison_file}:")
                                        st.dataframe(comp_results[relevant_cols])
                                        
                                        # Compare the counts
                                        st.write(f"**Difference in feature count:** {len(results) - len(comp_results)}")
                                    else:
                                        st.info(f"No features found with GO term {go_id} in {comparison_file}")
                            else:
                                st.info(f"No features found with GO term {go_id}")
            else:
                bio_term = st.text_input("Enter biological term (e.g., aluminum binding)")
                if bio_term:
                    if st.button("Search for related genes"):
                        # Get related GO terms
                        with st.spinner("Searching for GO terms..."):
                            go_terms = go_handler.search_go_terms_in_ontology_db(bio_term)
                        
                        if go_terms:
                            st.write(f"Found {len(go_terms)} related GO terms:")
                            terms_df = pd.DataFrame(go_terms)
                            st.dataframe(terms_df)
                            
                            # Search for genes with these terms
                            with st.spinner("Searching for genes..."):
                                results = go_handler.search_semantic_matches(bio_term, active_df)
                            
                            if not results.empty:
                                st.write(f"Found {len(results)} features related to '{bio_term}':")
                                relevant_cols = ['seqid', 'type', 'start', 'end', 'ID']
                                if 'Note' in results.columns:
                                    relevant_cols.append('Note')
                                if 'Ontology_term' in results.columns:
                                    relevant_cols.append('Ontology_term')
                                st.dataframe(results[relevant_cols])
                                
                                # Compare with other file if comparison mode enabled
                                if comparison_mode:
                                    comp_results = go_handler.search_semantic_matches(bio_term, comparison_df)
                                    if not comp_results.empty:
                                        st.write(f"Found {len(comp_results)} features related to '{bio_term}' in {comparison_file}:")
                                        st.dataframe(comp_results[relevant_cols])
                                    else:
                                        st.info(f"No features found related to '{bio_term}' in {comparison_file}")
                            else:
                                st.info(f"No genes found related to '{bio_term}'")
                        else:
                            st.info(f"No GO terms found related to '{bio_term}'")
        
        # Query Tab
        with main_tabs[3]:
            st.header("Ask Genoma!")
            st.write("Create a prompt about the genome data")
            
            # Model selection
            model_name = st.sidebar.selectbox(
                "Select Ollama Model",
                ["llama3.1:latest", "deepseek-r1:14b", "qwen2.5:7b"],
                index=1
            )
            temperature = st.sidebar.slider("Temperature", 0.0, 1.0, 0.1)
            
            # Query input
            query = st.text_input("Enter your prompt about the genome data:", 
                                placeholder="e.g., How many genes are on chromosome Tm_TA299_1A?")
            
            if query:
                with st.spinner("Processing prompt..."):
                    try:
                        # Create LLM agent with proper configuration
                        llm = OllamaLLM(model=model_name, temperature=temperature)
                        agent = create_pandas_dataframe_agent(
                            llm, 
                            active_df, 
                            verbose=True,
                            allow_dangerous_code=True,  # Allow code execution
                            handle_parsing_errors=True  # Handle parsing errors
                        )
                        
                        # Additional processing to handle common search cases
                        if "aluminum" in query.lower() or "aluminium" in query.lower():
                            st.write("Performing specialized search for aluminum-related genes...")
                            # Direct search in Note and other relevant columns
                            search_cols = ['Note', 'interPro', 'pFAM', 'Ontology_term']
                            search_terms = ['aluminum', 'aluminium', 'Al', 'metal']
                            
                            results = pd.DataFrame()
                            for col in search_cols:
                                if col in active_df.columns:
                                    for term in search_terms:
                                        matched = active_df[active_df[col].str.contains(term, case=False, na=False)]
                                        results = pd.concat([results, matched])
                            
                            if len(results) > 0:
                                st.write(f"Found {len(results)} features related to aluminum/metal interaction:")
                                st.dataframe(results[['seqid', 'type', 'start', 'end', 'ID', 'Note']])
                                
                                # Compare with other file if comparison mode enabled
                                if comparison_mode:
                                    comp_results = pd.DataFrame()
                                    for col in search_cols:
                                        if col in comparison_df.columns:
                                            for term in search_terms:
                                                matched = comparison_df[comparison_df[col].str.contains(term, case=False, na=False)]
                                                comp_results = pd.concat([comp_results, matched])
                                    
                                    if len(comp_results) > 0:
                                        st.write(f"Found {len(comp_results)} features related to aluminum/metal interaction in {comparison_file}:")
                                        st.dataframe(comp_results[['seqid', 'type', 'start', 'end', 'ID', 'Note']])
                        
                        # Run the query using the agent
                        response = agent.run(query)
                        st.write(response)
                    except Exception as e:
                        st.error(f"Error: {str(e)}")
                        st.info("Make sure Ollama is running with the selected model available.")

        # Analyses
        with main_tabs[4]:
            add_predefined_analyses(active_df)

        # Code Editor
        with main_tabs[5]:
            # If in comparison mode, pass both dataframes
            if comparison_mode:
                add_custom_code_editor(active_df, comparison_df, comparison_mode=True)
            else:
                add_custom_code_editor(active_df)


    # Clean up temporary files when app is closed
    for tmp_file in st.session_state.tmp_files:
        if os.path.exists(tmp_file):
            try:
                os.unlink(tmp_file)
            except:
                pass  # Ignore cleanup errors

if __name__ == "__main__":
    run_genome_app()

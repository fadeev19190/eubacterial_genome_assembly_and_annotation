import streamlit as st
import pandas as pd
import plotly.express as px
import glob
import plotly.graph_objects as go

# Set page configuration
st.set_page_config(
    page_title="Eubacterial Genome Assembly and Annotation",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS for styling
st.markdown("""
    <style>
        .main {
            background-color: #f8f9fa;
            color: #212529;
        }
        .sidebar .sidebar-content {
            background-color: #f8f9fa;
        }
        .reportview-container .main .block-container {
            padding: 20px;
            max-width: 1200px;
            margin: 0 auto;
        }
        .stButton > button {
            background-color: #4CAF50;
            color: white;
            border-radius: 5px;
            padding: 10px;
            border: none;
            cursor: pointer;
        }
        .stButton > button:hover {
            background-color: #45A049;
        }
        h1, h2, h3, h4, h5 {
            font-family: "Arial", sans-serif;
        }
    </style>
""", unsafe_allow_html=True)

# Title
st.title("Eubacterial Genome Assembly and Annotation")

st.markdown("""Student: Artem Fadeev

            Professor: Jan Pačes""")

# Introduction
st.markdown("""
### Project Overview
This project demonstrates the process of assembling and annotating a eubacterial genome. The workflow consists of:

1. Quality checking input sequences.
2. Cleaning and trimming sequences.
3. Genome assembly using advanced computational methods.
4. Comparative analysis of assemblies.
5. Gene prediction and annotation.
6. Visual representation of findings.

The methods and tools highlighted in this application reflect standard practices in bioinformatics, with additional customization for educational and research purposes.
""")

# Step 1: Quality Check
st.header("1. Quality Check of Input Sequences")

st.markdown("To ensure high-quality input data, we performed a quality check using FastQC.")

# Display FastQC Reports
st.markdown("**FastQC Reports (Before Trimming):**")
report_files = sorted(glob.glob("fastqc_before/*.html"))
report_file_names = [file.split('/')[-1] for file in report_files]
selected_report = st.selectbox("Select a FastQC report to view:", report_file_names)
selected_report_path = next((file for file in report_files if selected_report in file), "")

if selected_report_path:
    with open(selected_report_path, 'r') as file:
        report_content = file.read()
    st.components.v1.html(report_content, height=600, scrolling=True)

# Step 2: Data Cleaning
st.header("2. Data Cleaning and Trimming")

st.markdown("""
Cleaning and trimming are essential to remove adapters and low-quality sequences, ensuring reliable downstream analysis. Tools used include:

- **Trim Galore** for initial cleaning.
- **Cutadapt** for trimming adapters and adjusting sequence lengths.

**Trimming Rules:**
- For single-end sequences: Removed sequences shorter than 50 bp.
- For paired-end sequences: Trimmed sequences shorter than 50 bp.
- For mate-pair reads (3kb): Trimmed sequences between 50–400 bp and removed adapters.
""")

st.markdown("**FastQC Reports (After Trimming):**")
trimmed_report_files = sorted(glob.glob("fastqc_after/*.html"))
report_file_names_trimmed = [file.split('/')[-1] for file in trimmed_report_files]
selected_trimmed_report = st.selectbox("Select a trimmed FastQC report to view:", report_file_names_trimmed)
selected_report_path_trimmed = next((file for file in trimmed_report_files if selected_trimmed_report in file), "")

if selected_report_path_trimmed:
    with open(selected_report_path_trimmed, 'r') as file:
        report_content = file.read()
    st.components.v1.html(report_content, height=600, scrolling=True)

# Step 3: Genome Assembly
st.header("3. Genome Assembly")

st.markdown("""
Genome assembly combines sequencing reads into longer contiguous sequences. This project used two methods:

1. **Overlap-Layout-Consensus (OLC)** via Newbler.
2. **De Bruijn Graphs** via Velvet and SOAPdenovo.

**Assembly Variants:**
- **Newbler:** Single-end (overlap) and paired-end (no overlap).
- **SOAPdenovo:** Without mate-pair reads.
- **Newbler with Mate Pairs:** Incorporates 3kb mate-pair reads.
""")

# Step 4: Comparison of Assemblies
st.header("4. Comparison of Assemblies")

comparison_data = {
    "Assembly": ["Newbler", "Newbler with Long Reads", "SOAP-denovo"],
    "Total Contigs": [258, 247, 7327],
    "Total Length (bp)": [7800775, 7800100, 8877951],
    "Average Contig Length (bp)": [30235.56, 31579.35, 1211.68],
    "N50 (bp)": [100802, 107247, 6117],
    "L50": [25, 23, 436]
}
comparison_df = pd.DataFrame(comparison_data)
st.table(comparison_df)

st.markdown("""
### Key Metrics:
- **N50:** Indicates assembly continuity. Higher values signify longer contigs.
- **L50:** Indicates fragmentation. Lower values mean fewer contigs cover 50% of the genome.

**Conclusion:** The assembly using Newbler with long reads performed best, showing higher continuity and lower fragmentation.
""")

# Quast Report
st.markdown("**Assembly Quality Report:**")
report_path = "report.html"
try:
    with open(report_path, 'r') as file:
        report_content = file.read()
    # Add inline styling for a white background
    styled_report_content = f"""
    <html>
    <head>
        <style>
            body {{
                background-color: white !important;
                color: black !important;
            }}
        </style>
    </head>
    <body>
        {report_content}
    </body>
    </html>
    """
    st.components.v1.html(styled_report_content, height=600, scrolling=True)
except FileNotFoundError:
    st.error("The file 'report.html' was not found.")


# Step 5: Gene Prediction and Annotation
st.header("5. Gene Prediction and Annotation")

st.markdown("""
**Gene Prediction:** Using **Glimmer**, we identified genes and coding sequences (CDS). Predicted genes were further annotated using BLAST searches against reference databases.

**Visualization:** Gene length distribution for scaffolds is shown below.
""")

# Gene Length Visualization
file_path = "run_tres.predict"

def extract_gene_lengths(file_path):
    gene_lengths = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith(">"):
                fields = line.strip().split()
                start = int(fields[1])
                end = int(fields[2])
                gene_length = abs(end - start) + 1
                gene_lengths.append(gene_length)
    return gene_lengths

predicted_gene_lengths = extract_gene_lengths(file_path)
df = pd.DataFrame(predicted_gene_lengths, columns=["Predicted Gene Length"])
df = df[df["Predicted Gene Length"] <= 100000]

fig = px.box(df, y="Predicted Gene Length", title="Gene Length Distribution",
             labels={"Predicted Gene Length": "Gene Length (bp)"},
             height=600, width=800)
fig.update_layout(
    title_font_size=20,
    yaxis_title_font_size=16,
    xaxis_tickfont_size=14,
    template='plotly_white'
)

st.plotly_chart(fig, use_container_width=True)

st.markdown("### BLAST Results:")

st.markdown("""Based on the latest blasta results I can conclude, that the sequences belongs to a Burkholderia cenocepacia.
I limited the query to this organism only to save time and improve performance""")

# Define column names for the BLAST results
column_names = [
    "Query ID", "Subject ID", "Percent Identity", "Alignment Length",
    "E-value", "Bit Score", "Subject Accession", "Subject Title"
]

file_path = 'sorted_file.tsv'
try:
    # Load the data and assign column names
    df = pd.read_csv(file_path, sep='\t', names=column_names, header=None)
    # Display the dataframe with wider width
    st.dataframe(df, use_container_width=True)
except FileNotFoundError:
    st.error("The file 'sorted_file.tsv' was not found.")

# Challenges
st.markdown("### Challenges and Difficulties:")

st.markdown("""
1. At the beginning of the project, I struggled with understanding certain commands and likely made some mistakes due to inexperience.
2. The results could have been improved with more time spent on cleaning sequences before assembly.
3. The virtual machine (VM) was slow, particularly during library installations, which impacted workflow efficiency.
4. Managing files across multiple environments was confusing and led to occasional errors.
5. For future projects, I plan to set up a dedicated environment on my local machine to avoid such issues.
6. I need to familiarize myself better with the key commands and programs. For example, I am still struggling to figure out how to use the IGV program effectively.
""")

# Conclusion
st.markdown("### Conclusion")
st.markdown("""
This project provided a hands-on understanding of genome assembly and annotation workflows. 
Despite challenges, the results demonstrate the value of computational tools in bioinformatics and highlight areas for improvement in future projects.
""")

# Display three images in a row
col1, col2, col3 = st.columns(3)

with col1:
    st.image('1.jpg', use_container_width=True, caption="Image 1")
with col2:
    st.image('2.jpg', use_container_width=True, caption="Image 2")
with col3:
    st.image('3.jpg', use_container_width=True, caption="Image 3")

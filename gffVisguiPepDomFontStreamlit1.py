import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import streamlit as st
from io import StringIO
import random

# Function to parse the GFF file
def parse_gff(uploaded_file):
    domain_data = []
    protein_names = set()
    
    # Read the contents of the uploaded file
    content = uploaded_file.getvalue().decode("utf-8")
    
    for line in content.splitlines():
        if line.startswith("##"):  # Skip comment lines
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9:
            continue  # Skip malformed lines
        
        protein_name = fields[0]
        feature_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        attributes = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
        
        if feature_type == "polypeptide":
            protein_names.add(protein_name)
        
        if feature_type == "protein_match":
            domain_name = attributes.get("Name", "Unknown")
            domain_data.append((protein_name, start, end, domain_name))
    
    return domain_data, protein_names

# Function to generate random color in hex format
def random_color():
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

# Function to plot the protein domains with selected colors
def plot_domains(domain_data, selected_proteins, selected_domains, shape_choice, domain_colors):
    filtered_data = [
        entry for entry in domain_data 
        if entry[0] in selected_proteins and entry[3] in selected_domains
    ]
    
    df = pd.DataFrame(filtered_data, columns=["Protein", "Start", "End", "Domain"])

    if df.empty:
        st.info("No domain data found for the selected proteins and domains!")
        return

    fig, ax = plt.subplots(figsize=(10, len(selected_proteins) * 1.5))

    protein_height = 0.4  
    protein_gap = 1.0     

    for i, protein in enumerate(selected_proteins):
        y_position = i * (protein_height + protein_gap)
        protein_rows = df[df['Protein'] == protein]

        for _, row in protein_rows.iterrows():
            width = row["End"] - row["Start"]
            alpha_value = 0.6  # Set transparency for overlapping domains

            domain_color = domain_colors.get(row["Domain"], "#1f77b4")  # Default to blue if not specified

            if shape_choice == "Rectangle":
                ax.add_patch(mpatches.Rectangle((row["Start"], y_position),
                                                 width, protein_height,
                                                 color=domain_color,
                                                 alpha=alpha_value,  # Set transparency
                                                 label=row["Domain"] if row["Domain"] not in ax.get_legend_handles_labels()[1] else ""))
            elif shape_choice == "Rounded Rectangle":
                ax.add_patch(FancyBboxPatch((row["Start"], y_position),
                                            width, protein_height,
                                            boxstyle="round,pad=0.02,rounding_size=0.15",
                                            color=domain_color,
                                            alpha=alpha_value,  # Set transparency
                                            label=row["Domain"] if row["Domain"] not in ax.get_legend_handles_labels()[1] else ""))
            elif shape_choice == "Oval":
                ax.add_patch(mpatches.Ellipse(((row["Start"] + row["End"]) / 2, y_position + protein_height / 2),
                                              width, protein_height,
                                              color=domain_color,
                                              alpha=alpha_value,  # Set transparency
                                              label=row["Domain"] if row["Domain"] not in ax.get_legend_handles_labels()[1] else ""))

    ax.set_xlim(0, df["End"].max() + 10)
    ax.set_ylim(-protein_gap, len(selected_proteins) * (protein_height + protein_gap))

    ax.set_yticks([(i * (protein_height + protein_gap)) + (protein_height / 2) for i in range(len(selected_proteins))])
    ax.set_yticklabels(selected_proteins)
    ax.set_xlabel("Position on Protein Sequence")
    ax.set_title("Protein Domains Visualization")

    handles = [mpatches.Patch(color=domain_colors[domain], label=domain, alpha=alpha_value) for domain in selected_domains]
    ax.legend(handles=handles, title="Domains", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

    plt.grid(axis='x')
    plt.tight_layout()
    st.pyplot(fig)  # Use Streamlit's method to show the figure

# Streamlit application
st.title("Protein Domain Visualizer")

uploaded_file = st.file_uploader("Select GFF File", type="gff")

if uploaded_file is not None:
    domain_data, protein_names = parse_gff(uploaded_file)

    selected_proteins = st.multiselect("Select Protein(s):", options=list(protein_names))

    if selected_proteins:
        selected_domains = st.multiselect("Select Domain(s):", options=list(set(entry[3] for entry in domain_data)))

        if selected_domains:
            shape_choice = st.selectbox("Select Domain Shape:", ["Rectangle", "Rounded Rectangle", "Oval"])

            # Initialize color storage in session state if not already done
            if 'domain_colors' not in st.session_state:
                st.session_state.domain_colors = {domain: random_color() for domain in selected_domains}

            # Ensure all selected domains have a color (avoid resetting on re-selection)
            for domain in selected_domains:
                if domain not in st.session_state.domain_colors:
                    st.session_state.domain_colors[domain] = random_color()

            st.write("You can customize the domain colors below:")
            for domain in selected_domains:
                st.session_state.domain_colors[domain] = st.color_picker(f"Pick a color for {domain}", st.session_state.domain_colors[domain])

            if st.button("Visualize Selected Proteins"):
                plot_domains(domain_data, selected_proteins, selected_domains, shape_choice, st.session_state.domain_colors)

import streamlit as st
import pandas as pd
import requests
import gseapy as gp
from Bio import Entrez
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import urllib.parse
import plotly.express as px


Entrez.email = "blanca.puechegranados@usp.ceu.es"
GRAPHQL_URL = "https://dgidb.org/api/graphql"
ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/id/"

st.markdown("""
    <style>
    .main {
        background: linear-gradient(to bottom right, #e3f2fd, #fce4ec);
        font-family: 'Arial', sans-serif;
        padding: 20px;
    }
    .stButton > button {
        background-color: #6a1b9a;
        color: white;
        font-weight: bold;
        border-radius: 10px;
    }
    .stDownloadButton > button {
        background-color: #2e7d32;
        color: white;
        font-weight: bold;
        border-radius: 10px;
    }
    </style>
""", unsafe_allow_html=True)

# Methods
def get_disease_name(mesh_id):
    try:
        search_handle = Entrez.esearch(db="mesh", term=mesh_id)
        search_record = Entrez.read(search_handle)
        search_handle.close()
        if not search_record['IdList']:
            return None
        uid = search_record['IdList'][0]
        summary_handle = Entrez.esummary(db="mesh", id=uid)
        summary_record = Entrez.read(summary_handle)
        summary_handle.close()
        return summary_record[0]['DS_MeshTerms'][0]
    except:
        return None

def get_gene_name_from_ensembl(ensembl_id):
    response = requests.get(f"{ENSEMBL_LOOKUP_URL}{ensembl_id}?content-type=application/json")
    if response.status_code == 200:
        data = response.json()
        return data.get("display_name", "Not Found")
    return "Not Found"

def fetch_gene_names(df):
    df["Gene Name"] = df["Gene"].apply(get_gene_name_from_ensembl)
    df_filtered = df[df["Gene Name"] != "Not Found"]

    # Sum log2 fold changes for repeated genes
    grouped = df_filtered.groupby(["Gene", "Gene Name"], as_index=False).agg({
        "log_2 fold change": "sum"
    })

    return grouped


def analyze_pathways(df):
    gene_list = df["Gene Name"].dropna().unique().tolist()
    enr = gp.enrichr(gene_list=gene_list, gene_sets="Reactome_2016", organism="Human", outdir=None)
    if enr.results.empty:
        return None, None
    #top_pathways = enr.results.sort_values("Adjusted P-value").head(5)
    top_pathways = enr.results.sort_values("Adjusted P-value").head(5).copy()

    top_pathways["Reactome Link"] = top_pathways["Term"].apply(
        lambda term: f'<a href="https://reactome.org/content/query?q={urllib.parse.quote(term)}" target="_blank">{term}</a>'
    )

    top_pathway_genes = top_pathways.iloc[0]["Genes"].split(";")
    important_genes = df[df["Gene Name"].isin(top_pathway_genes)].copy()
    important_genes["abs_fc"] = important_genes["log_2 fold change"].abs()
    important_genes = important_genes.sort_values("abs_fc", ascending=False)
    return top_pathways, important_genes


def get_drug_targets_dgidb_graphql(gene_names):
    all_results = []
    for gene_name in gene_names:
        graphql_query = f"""
        {{
          genes(names: ["{gene_name}"]) {{
            nodes {{
              interactions {{
                drug {{
                  name
                  conceptId
                }}
                interactionScore
                interactionTypes {{
                  type
                  directionality
                }}
                interactionAttributes {{
                  name
                  value
                }}
                publications {{
                  pmid
                }}
                sources {{
                  sourceDbName
                }}
              }}
            }}
          }}
        }}
        """
        response = requests.post(GRAPHQL_URL, json={"query": graphql_query})
        if response.status_code == 200:
            try:
                data = response.json()
                if 'data' in data and 'genes' in data['data'] and len(data['data']['genes']['nodes']) > 0:
                    interactions = data['data']['genes']['nodes'][0].get('interactions', [])
                    for interaction in interactions:
                        drug_name = interaction['drug'].get('name', 'Unknown')
                        score = interaction.get('interactionScore', 'N/A')
                        types = interaction.get('interactionTypes', [])
                        interaction_type = types[0].get('type', 'N/A') if types else 'N/A'
                        directionality = types[0].get('directionality', 'N/A') if types else 'N/A'
                        sources = interaction.get('sources', [])
                        source = sources[0]['sourceDbName'] if sources else 'N/A'
                        pmids = interaction.get('publications', [])
                        pmid = pmids[0]['pmid'] if pmids else 'N/A'

                        all_results.append({
                            'Gene': gene_name,
                            'Drug': drug_name,
                            'Interaction Type': interaction_type,
                            'Directionality': directionality,
                            'Source': source,
                            'PMID': pmid,
                            'Interaction Score': score
                        })
            except:
                continue
    
    return pd.DataFrame(all_results)

def drug_with_links(df):
    df_with_links = df.copy()
    df_with_links["Drug"] = df_with_links["Drug"].apply(
        lambda drug: f'<a href="https://dgidb.org/results?searchType=drug&searchTerms={urllib.parse.quote(drug)}" target="_blank">{drug}</a>'
    )
    return df_with_links


# Interface

st.markdown(
    """
    <style>
    body {
        color: #000000;
    }
    .stApp {
        background-image: linear-gradient(to bottom, #94c3dc, #f0f0f0, #ffffff);
        background-size: cover;
    }    
    </style>
    """,
    unsafe_allow_html=True
)

st.title("🧬 Disease-Gene-Drug Analysis Dashboard")

# Step 1: Enter MeSH ID
mesh_id = st.text_input("🔍 Enter MeSH ID (e.g., D003920 for Diabetes Mellitus):")

if mesh_id:
    with st.spinner("Fetching disease information..."):
        disease = get_disease_name(mesh_id)
    if disease:
        st.success(f"🎯 Disease identified: **{disease}**")
        
        # Step 2: File upload only after disease is found
        uploaded_file = st.file_uploader("📁 Upload Differential Expression File (TSV from Expression Atlas)", type=["tsv"])
        
        if uploaded_file:
            df_raw = pd.read_csv(uploaded_file, sep="\t")
            st.subheader("📊 Raw Uploaded Data")
            st.dataframe(df_raw.head())
            
            st.subheader("⚙️ Threshold Settings")
            log2fc_threshold = st.slider("Minimum absolute log2 fold change", min_value=0.0, max_value=5.0, value=1.0, step=0.1)
            pval_threshold = st.slider("Adjusted p-value threshold", min_value=0.0001, max_value=0.1, value=0.05, step=0.001)

            
            st.subheader("📊 Volcano Plot (log2FC vs. -log10 p-value)")

            if "Adjusted p-value" in df_raw.columns:
                try:
                    df_raw["neg_log10_p"] = -np.log10(df_raw["Adjusted p-value"])
                    df_raw["Significant"] = (df_raw["Adjusted p-value"] < pval_threshold) & (df_raw["log_2 fold change"].abs() > log2fc_threshold)

                    fig = px.scatter(
                        df_raw,
                        x="log_2 fold change",
                        y="neg_log10_p",
                        color="Significant",
                        hover_data=["Gene"],
                        labels={
                            "log_2 fold change": "log₂ Fold Change",
                            "neg_log10_p": "-log₁₀ Adjusted p-value"
                        },
                        title="Volcano Plot (Interactive)"
                    )

                    fig.add_hline(y=-np.log10(pval_threshold), line_dash="dash", line_color="blue")
                    fig.add_vline(x=log2fc_threshold, line_dash="dash", line_color="green")
                    fig.add_vline(x=-log2fc_threshold, line_dash="dash", line_color="green")

                    st.plotly_chart(fig, use_container_width=True)

                except Exception as e:
                    st.warning(f"Could not generate interactive volcano plot: {e}")
            else:
                st.info("No 'p_value' column found for volcano plot.")
    

            with st.spinner("🔄 Mapping Ensembl IDs to gene names..."):
                df_selected = fetch_gene_names(df_raw)
            # Convert the Ensembl gene IDs into HTML links
            df_selected_with_links = df_selected.copy()
            df_selected_with_links["Gene"] = df_selected_with_links["Gene"].apply(
                lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
            )
            df_selected_with_links = df_selected_with_links.sort_values(by="log_2 fold change", ascending=False)

            st.markdown("### 🧬 Gene Table with Links to Ensembl")
            st.write(df_selected_with_links.to_html(escape=False, index=False), unsafe_allow_html=True)
            st.download_button("📥 Download Genes CSV", df_selected_with_links.to_csv(index=False), "genes.csv", "text/csv")


            with st.spinner("🧠 Performing pathway analysis..."):
                top_pathways, important_genes = analyze_pathways(df_selected)
                top_pathways["-log10(Adj P)"] = -np.log10(top_pathways["Adjusted P-value"])
            
            
            if top_pathways is not None:
                st.subheader("📈 Top Enriched Reactome Pathways")
                #st.dataframe(top_pathways)
                st.write(top_pathways[["Reactome Link", "Adjusted P-value", "-log10(Adj P)"]].to_html(escape=False, index=False), unsafe_allow_html=True)


                selected_pathway = st.selectbox(
                    "🔍 Select a pathway to see the top genes",
                top_pathways["Term"].tolist(),
                key="pathway"
                )

                pathway_genes = important_genes[important_genes["Gene Name"].isin(
                    top_pathways[top_pathways["Term"] == selected_pathway]["Genes"].iloc[0].split(";")
                )]
    
                st.subheader(f"🔥 Important genes in pathway: {selected_pathway}")
                pathway_genes["Gene"] = pathway_genes["Gene"].apply(
                    lambda gene_id: f'<a href="https://www.ensembl.org/Multi/Search/Results?q={gene_id}" target="_blank">{gene_id}</a>'
                )
                st.write(pathway_genes.to_html(escape=False, index=False), unsafe_allow_html=True)
                st.download_button("📥 Download Important Genes CSV", pathway_genes.to_csv(index=False), "genes_pathway.csv", "text/csv")

            else:
                st.warning("⚠️ No enriched pathways found.")

            with st.spinner("💊 Searching drug-gene interactions..."):
                selected_pathway = st.selectbox(
                    "🔍 Select a pathway to see drugs",
                top_pathways["Term"].tolist(),
                key="drugs"
                )
                
                pathway_genes = important_genes[important_genes["Gene Name"].isin(
                    top_pathways[top_pathways["Term"] == selected_pathway]["Genes"].iloc[0].split(";")
                )]
                                
                drug_df = get_drug_targets_dgidb_graphql(pathway_genes["Gene Name"].tolist())

            if not drug_df.empty:
                st.subheader("💊 Drug-Gene Interactions")
                
                drug_df_with_links = drug_with_links(drug_df)
    
                st.write(drug_df_with_links.to_html(escape=False, index=False), unsafe_allow_html=True)
    
                st.download_button("📥 Download Drug Interactions CSV", drug_df.to_csv(index=False), "drug_interactions.csv", "text/csv")
            else:
                st.warning("🚫 No drug-gene interactions found.")
    else:
        st.error("❌ Could not find a disease with the provided MeSH ID.")

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np

st.set_page_config(page_title="VCF Mutation Explorer", layout="wide")

st.title("🧬 VCF Mutation Explorer Dashboard")

files = st.file_uploader("Upload VCF files", accept_multiple_files=True)

mutations = []

if files:

    for file in files:
        sample = file.name

        for line in file:
            line = line.decode("utf-8")

            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")

            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            info = cols[7]

            gene = "Unknown"
            effect = "Unknown"

            if "ANN=" in info:
                ann = info.split("ANN=")[1].split("|")
                if len(ann) > 3:
                    effect = ann[1]
                    gene = ann[3]

            mutations.append([sample, chrom, pos, ref, alt, gene, effect])

    df = pd.DataFrame(
        mutations,
        columns=["Sample","Chrom","Position","REF","ALT","Gene","Effect"]
    )

    samples = df["Sample"].unique()

    # ======================
    # DATASET STATISTICS
    # ======================

    st.subheader("Dataset Statistics")

    c1, c2, c3, c4 = st.columns(4)

    c1.metric("Samples", df["Sample"].nunique())
    c2.metric("Total Mutations", len(df))
    c3.metric("Unique Mutations", df.drop_duplicates(["Chrom","Position","REF","ALT"]).shape[0])
    c4.metric("Genes Detected", df["Gene"].nunique())

    # ======================
    # MUTATION TABLE
    # ======================

    st.subheader("Mutation Table")
    st.dataframe(df)

    # ======================================================
    # COLISTIN RESISTANCE GENE PANEL
    # ======================================================

    st.subheader("Colistin Resistance Gene Panel")

    resistance_genes = {
        "mgrB": "ACOXO5_RS08420",
        "phoP": "ACOXO5_RS15300",
        "phoQ": "ACOXO5_RS15305",
        "pmrA": "ACOXO5_RS17135",
        "pmrB": "ACOXO5_RS17140",
        "pmrD": "ACOXO5_RS07870",
        "arnA": "ACOXO5_RS27340",
        "arnB": "ACOXO5_RS27330",
        "arnC": "ACOXO5_RS27335",
        "arnT": "ACOXO5_RS27350"
    }

    resistance_df = df[df["Gene"].isin(resistance_genes.values())].copy()

    gene_map = {v:k for k,v in resistance_genes.items()}
    resistance_df["Gene_Name"] = resistance_df["Gene"].map(gene_map)

    if not resistance_df.empty:

        st.write("Detected Mutations in Resistance Genes")
        st.dataframe(resistance_df)

        gene_counts = (
            resistance_df.groupby("Gene_Name")["Sample"]
            .nunique()
            .reset_index(name="Sample_Count")
        )

        fig_res = px.bar(
            gene_counts,
            x="Gene_Name",
            y="Sample_Count",
            color="Sample_Count",
            color_continuous_scale="Reds",
            title="Samples with Mutations in Resistance Genes"
        )

        st.plotly_chart(fig_res)

        st.subheader("Mutation Effects in Resistance Genes")

        effect_counts = (
            resistance_df.groupby(["Gene_Name","Effect"])
            .size()
            .reset_index(name="Count")
        )

        fig_effect = px.bar(
            effect_counts,
            x="Gene_Name",
            y="Count",
            color="Effect",
            barmode="stack"
        )

        st.plotly_chart(fig_effect)



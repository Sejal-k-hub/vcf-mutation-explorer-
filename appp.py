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

    c1,c2,c3,c4 = st.columns(4)

    c1.metric("Samples", df["Sample"].nunique())
    c2.metric("Total Mutations", len(df))
    c3.metric("Unique Mutations", df.drop_duplicates(["Chrom","Position","REF","ALT"]).shape[0])
    c4.metric("Genes Detected", df["Gene"].nunique())

    # ======================
    # DOWNLOAD TABLE
    # ======================

    st.download_button(
        "Download Full Mutation Table",
        df.to_csv(index=False),
        file_name="mutation_table.csv"
    )

    # ======================
    # MUTATION TABLE
    # ======================

    st.subheader("Mutation Table")
    st.dataframe(df)

    # ======================
    # MUTATION IMPACT
    # ======================

    st.subheader("Mutation Impact Summary")

    impact_counts = df["Effect"].value_counts().reset_index()
    impact_counts.columns = ["Effect","Count"]

    fig_imp = px.bar(
        impact_counts,
        x="Effect",
        y="Count",
        color="Count",
        color_continuous_scale=px.colors.sequential.Turbo,
        template="plotly_dark"
    )

    st.plotly_chart(fig_imp, use_container_width=True)

    # ======================
    # PIE CHART
    # ======================

    st.subheader("Mutation Type Distribution")

    fig_pie = px.pie(
        impact_counts,
        names="Effect",
        values="Count",
        color_discrete_sequence=px.colors.qualitative.Bold
    )

    st.plotly_chart(fig_pie, use_container_width=True)

    # ======================
    # TOP MUTATED GENES
    # ======================

    st.subheader("Top Mutated Genes")

    gene_counts = (
        df.groupby("Gene")
        .size()
        .reset_index(name="Mutation_Count")
        .sort_values("Mutation_Count",ascending=False)
        .head(20)
    )

    fig_gene = px.bar(
        gene_counts,
        x="Gene",
        y="Mutation_Count",
        color="Mutation_Count",
        color_continuous_scale=px.colors.sequential.Viridis
    )

    st.plotly_chart(fig_gene, use_container_width=True)

    # ======================
    # MUTATION LOAD
    # ======================

    st.subheader("Mutation Load per Sample")

    sample_counts = df.groupby("Sample").size().reset_index(name="Mutation_Count")

    fig_sample = px.bar(
        sample_counts,
        x="Sample",
        y="Mutation_Count",
        color="Mutation_Count",
        color_continuous_scale=px.colors.sequential.Plasma
    )

    st.plotly_chart(fig_sample, use_container_width=True)

    # ======================
    # GENE HEATMAP
    # ======================

    st.subheader("Gene Mutation Heatmap")

    top_genes = df["Gene"].value_counts().head(25).index

    heatmap_data = pd.crosstab(
        df[df["Gene"].isin(top_genes)]["Gene"],
        df["Sample"]
    )

    fig_heat = px.imshow(
        heatmap_data,
        color_continuous_scale=px.colors.sequential.Inferno,
        aspect="auto",
        text_auto=True
    )

    st.plotly_chart(fig_heat, use_container_width=True)

    # ======================
    # EFFECT HEATMAP
    # ======================

    st.subheader("Mutation Effect Heatmap")

    effect_heat = pd.crosstab(df["Effect"], df["Sample"])

    fig_effect_heat = px.imshow(
        effect_heat,
        color_continuous_scale=px.colors.sequential.Magma,
        aspect="auto",
        text_auto=True
    )

    st.plotly_chart(fig_effect_heat, use_container_width=True)

    # ======================================================
    # COLISTIN RESISTANCE GENES
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
            color_continuous_scale=px.colors.sequential.Reds
        )

        st.plotly_chart(fig_res, use_container_width=True)

        heat = pd.crosstab(resistance_df["Gene_Name"], resistance_df["Sample"])

        fig_res_heat = px.imshow(
            heat,
            color_continuous_scale=px.colors.sequential.Reds,
            aspect="auto",
            text_auto=True
        )

        st.plotly_chart(fig_res_heat, use_container_width=True)

    else:
        st.info("No resistance mutations detected")

    # ======================
    # RESISTANCE PREDICTION
    # ======================

    st.subheader("Colistin Resistance Prediction")

    predictions = []

    for sample in samples:

        sample_mut = resistance_df[resistance_df["Sample"] == sample]

        if not sample_mut.empty:
            prediction = "Potential Resistant"
        else:
            prediction = "No Mutation Detected"

        predictions.append([sample, prediction])

    pred_df = pd.DataFrame(predictions, columns=["Sample","Prediction"])

    st.dataframe(pred_df)

    pred_counts = pred_df["Prediction"].value_counts().reset_index()
    pred_counts.columns = ["Prediction","Count"]

    fig_pred = px.bar(
        pred_counts,
        x="Prediction",
        y="Count",
        color="Prediction",
        color_discrete_map={
            "Potential Resistant":"red",
            "No Mutation Detected":"green"
        }
    )

    st.plotly_chart(fig_pred, use_container_width=True)

    # ======================
    # SAMPLE COMPARISON
    # ======================

    st.subheader("Sample Mutation Comparison")

    s1 = st.selectbox("Select Sample 1", samples)
    s2 = st.selectbox("Select Sample 2", samples)

    if s1 and s2:

        mut1 = set(df[df["Sample"]==s1]["Position"])
        mut2 = set(df[df["Sample"]==s2]["Position"])

        shared = len(mut1.intersection(mut2))
        unique1 = len(mut1-mut2)
        unique2 = len(mut2-mut1)

        comp_df = pd.DataFrame({
            "Category":["Shared","Unique Sample1","Unique Sample2"],
            "Count":[shared,unique1,unique2]
        })

        fig_comp = px.bar(
            comp_df,
            x="Category",
            y="Count",
            color="Category",
            color_discrete_sequence=px.colors.qualitative.Bold
        )

        st.plotly_chart(fig_comp, use_container_width=True)

    # ======================
    # GENE EXPLORER
    # ======================

    st.subheader("Gene Mutation Explorer")

    gene_search = st.text_input("Enter Gene Name")

    if gene_search:

        gene_df = df[df["Gene"].str.contains(gene_search, case=False, na=False)]

        st.write("Mutations found:", len(gene_df))

        st.dataframe(gene_df)

        freq = gene_df.groupby("Position").size().reset_index(name="Count")

        fig5 = px.bar(
            freq,
            x="Position",
            y="Count",
            color="Count",
            color_continuous_scale=px.colors.sequential.Viridis
        )

        st.plotly_chart(fig5, use_container_width=True)

        st.subheader("Gene Mutation Lollipop Plot")

        fig_lollipop = go.Figure()

        fig_lollipop.add_trace(
            go.Scatter(
                x=freq["Position"],
                y=[1]*len(freq),
                mode="markers",
                marker=dict(
                    size=freq["Count"]*6,
                    color=freq["Count"],
                    colorscale="Turbo",
                    showscale=True
                )
            )
        )

        fig_lollipop.update_layout(
            xaxis_title="Genome Position",
            yaxis=dict(showticklabels=False)
        )

        st.plotly_chart(fig_lollipop, use_container_width=True)

        gene_mut_search = st.text_input("Search mutation inside this gene")

        if gene_mut_search:

            gene_results = gene_df[
                gene_df.apply(
                    lambda row: gene_mut_search.lower() in str(row).lower(),
                    axis=1
                )
            ]

            st.dataframe(gene_results)

            st.download_button(
                "Download Gene Results",
                gene_results.to_csv(index=False),
                file_name="gene_mutations.csv"
            )

from pathlib import Path
import streamlit as st
import pandas as pd
import requests
import numpy as np
import plotly.graph_objects as go

# Try Plotly
try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

API_BASE = "https://srna-platform-backend.onrender.com"  # changed from http://backend:8080


# -----------------------------------------------------------------------------
# Streamlit page config
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="sRNA–mRNA Interaction Dashboard",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <h1 style='text-align:center; color:#1f4e79;'>
    🧬 sRNA–mRNA Interaction Intelligence Platform
    </h1>
    """,
    unsafe_allow_html=True,
)


# -----------------------------------------------------------------------------
# Helpers: API + plotting + formatting
# -----------------------------------------------------------------------------
def api_post(path: str, payload: dict):
    """Safely POST to backend and show user-friendly errors."""
    url = f"{API_BASE}{path}"
    try:
        resp = requests.post(url, json=payload, timeout=60)
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        st.error(
            f"❌ Could not reach backend at {url}. "
            f"Is uvicorn running on port 8080?\n\nError: {e}"
        )
        return None


def extract_hybrid_from_raw_row(raw_row: str):
    """
    IntaRNA CSV raw row:
      id1;start1;end1;id2;start2;end2;subseqDP;hybridDP;E
    We want subseqDP & hybridDP to build duplex viewers.
    """
    if not raw_row:
        return None, None
    parts = raw_row.split(";")
    if len(parts) < 8:
        return None, None
    subseq_dp = parts[-3]  # subseqDP
    hybrid_dp = parts[-2]  # hybridDP
    return subseq_dp, hybrid_dp


def build_ascii_duplex(subseq_dp: str, hybrid_dp: str) -> str:
    """ASCII duplex representation."""
    if not subseq_dp or not hybrid_dp:
        return "No duplex details available."

    try:
        target_seq, srna_seq = subseq_dp.split("&")
    except ValueError:
        target_seq = subseq_dp
        srna_seq = ""

    lines = []
    lines.append(f"target: 5' {target_seq} 3'")
    lines.append(f"hybrid:     {hybrid_dp}")
    if srna_seq:
        lines.append(f"sRNA  : 3' {srna_seq} 5'")
    return "\n".join(lines)


def make_dg_bar_chart(df: pd.DataFrame):
    if not HAS_PLOTLY:
        st.info("Install `plotly` to see ΔG bar chart.")
        return

    if "deltaG" not in df.columns:
        st.info("No ΔG values present in DataFrame.")
        return

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=df["rank"],
            y=df["deltaG"],
            text=[
                f"ML={round(x,3)}" if pd.notna(x) else "ML=NA"
                for x in df.get("ml_score", [])
            ],
            hovertemplate="Rank %{x}<br>ΔG=%{y}<extra></extra>",
        )
    )
    fig.update_layout(
        title="ΔG binding energies per predicted interaction",
        xaxis_title="Interaction rank (IntaRNA)",
        yaxis_title="ΔG (kcal/mol)",
    )
    st.plotly_chart(fig, use_container_width=True)


def make_violin_plot(df: pd.DataFrame):
    """Violin/box plot for ΔG and ML scores."""
    if not HAS_PLOTLY:
        st.info("Install `plotly` to see violin plots.")
        return

    fig = go.Figure()
    if "deltaG" in df.columns:
        fig.add_trace(
            go.Violin(
                y=df["deltaG"],
                name="ΔG",
                box_visible=True,
                meanline_visible=True,
            )
        )
    if "ml_score" in df.columns and df["ml_score"].notna().any():
        fig.add_trace(
            go.Violin(
                y=df["ml_score"],
                name="ML score",
                box_visible=True,
                meanline_visible=True,
            )
        )

    if not fig.data:
        st.info("No numeric columns available for violin plot.")
        return

    fig.update_layout(
        title="Distribution of ΔG and ML scores (violin plot)",
        yaxis_title="Value",
    )
    st.plotly_chart(fig, use_container_width=True)


def make_arc_figure(df: pd.DataFrame):
    """Gene-browser style curved arcs over mRNA coordinates."""
    if not HAS_PLOTLY:
        st.info("Install `plotly` to see arc plot.")
        return

    fig = go.Figure()

    for _, row in df.iterrows():
        start = row.get("target_start")
        end = row.get("target_end")
        if start is None or end is None:
            continue
        if start > end:
            start, end = end, start

        xs = np.linspace(start, end, 40)
        mid = 0.5 * (start + end)
        radius = (end - start) / 2.0
        ys = np.sqrt(np.maximum(radius**2 - (xs - mid) ** 2, 0.0))

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                name=f"rank {row.get('rank')}, ΔG={row.get('deltaG')}",
            )
        )

    fig.update_layout(
        title="Gene-browser style arcs over mRNA coordinates",
        xaxis_title="mRNA coordinate (target_start → target_end)",
        yaxis_title="Arc height (visual only)",
        showlegend=True,
    )

    st.plotly_chart(fig, use_container_width=True)

def make_graphical_duplex(subseq_dp: str, hybrid_dp: str):
    """
    Simple graphical duplex:
    - target strand at y=0
    - sRNA strand at y=1
    - lines connecting 'paired' bases from hybrid_dp.
    """
    if not HAS_PLOTLY:
        st.info("Install `plotly` for graphical duplex plot.")
        return

    if not subseq_dp or not hybrid_dp:
        st.info("No duplex details available for plot.")
        return

    try:
        target_seq, srna_seq = subseq_dp.split("&")
    except ValueError:
        target_seq = subseq_dp
        srna_seq = ""

    n = max(len(target_seq), len(hybrid_dp))

    xs = list(range(1, n + 1))
    fig = go.Figure()

    # Target strand
    fig.add_trace(
        go.Scatter(
            x=xs[: len(target_seq)],
            y=[0] * len(target_seq),
            mode="lines+markers",
            name="target",
        )
    )

    # sRNA strand
    if srna_seq:
        fig.add_trace(
            go.Scatter(
                x=xs[: len(srna_seq)],
                y=[1] * len(srna_seq),
                mode="lines+markers",
                name="sRNA",
            )
        )

    # Lines for paired positions (very simple heuristic: any '(' or ')' is a "paired" site)
    for i, ch in enumerate(hybrid_dp):
        if ch in ("(", ")"):
            x = i + 1
            fig.add_shape(
                type="line",
                x0=x,
                y0=0,
                x1=x,
                y1=1,
                line=dict(width=1),
            )

    fig.update_layout(
        title="Graphical duplex (schematic)",
        xaxis_title="Position in subseq",
        yaxis=dict(showticklabels=False),
        showlegend=True,
    )

    st.plotly_chart(fig, use_container_width=True)


def highlight_sequence_html(seq: str,
                            seed_info: dict | None,
                            bind_start: int | None,
                            bind_end: int | None) -> str:
    """
    HTML sequence highlight:
    - yellow background for seed region
    - cyan background for binding region
    (1-based coordinates)
    """
    if not seq:
        return "<i>No sequence</i>"

    seed_info = seed_info or {}
    seed_has = bool(seed_info.get("has_seed"))
    seed_len = seed_info.get("length") or 0
    seed_srna_start = seed_info.get("srna_start")
    seed_mrna_start = seed_info.get("mrna_start")

    # We don't know if this is sRNA or mRNA; caller passes correct start positions.
    # We'll treat seed_start as whichever is relevant (caller chooses).
    seed_start = seed_srna_start  # or mrna start separately by calling function accordingly

    html_chars = []
    for i, base in enumerate(seq, start=1):
        styles = []
        if seed_has and seed_start and seed_len and seed_start <= i < seed_start + seed_len:
            styles.append("background-color: #fff3b0")  # pale yellow
        if bind_start and bind_end and bind_start <= i <= bind_end:
            styles.append("background-color: #b3e6ff")  # pale cyan

        if styles:
            style_str = ";".join(styles)
            html_chars.append(f"<span style='{style_str}'>{base}</span>")
        else:
            html_chars.append(base)

    return "".join(html_chars)


# --- Dot-bracket parsing (RNAcofold scaffold) --------------------------------
def parse_dot_bracket_pairs(structure: str):
    """
    Parse a dot-bracket string into base-pair indices (1-based).
    Standard stack-based algorithm.
    """
    stack = []
    pairs = []
    for i, ch in enumerate(structure, start=1):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack:
                j = stack.pop()
                pairs.append((j, i))
    return pairs


def make_dotbracket_arc(sequence: str, structure: str):
    """
    Render a dot-bracket structure as base-pair arcs (RNAcofold scaffold).
    This only activates if the backend returns 'cofold' info in /explain.
    """
    if not HAS_PLOTLY:
        st.info("Install `plotly` to see dot-bracket arc plot.")
        return

    if not sequence or not structure:
        st.info("No sequence / structure for RNAcofold arc.")
        return

    pairs = parse_dot_bracket_pairs(structure)
    fig = go.Figure()

    # Plot sequence as baseline
    xs = list(range(1, len(sequence) + 1))
    fig.add_trace(
        go.Scatter(
            x=xs,
            y=[0] * len(xs),
            mode="lines",
            name="sequence",
        )
    )

    # Add arcs for base pairs
    for i, j in pairs:
        xs_arc = np.linspace(i, j, 40)
        mid = 0.5 * (i + j)
        radius = (j - i) / 2.0
        ys_arc = np.sqrt(np.maximum(radius**2 - (xs_arc - mid) ** 2, 0.0))
        fig.add_trace(
            go.Scatter(
                x=xs_arc,
                y=ys_arc,
                mode="lines",
                line=dict(width=1),
                showlegend=False,
            )
        )

    fig.update_layout(
        title="RNAcofold / dot-bracket arc plot (scaffold)",
        xaxis_title="Position in concatenated sequence",
        yaxis_title="Arc height (visual only)",
    )
    st.plotly_chart(fig, use_container_width=True)

# -----------------------------------------------------------------------------
# MAIN INPUT PANEL (center, not sidebar)
# -----------------------------------------------------------------------------
st.markdown("###  Sequence Input Panel")

c1, c2 = st.columns(2)
with c1:
    srna_name = st.text_input("sRNA Name", "MicF")
    srna_seq = st.text_area(
        "sRNA Sequence (5'→3')",
        "AUGCAGGUUCAUCAACACCUUAGGAAACCCUGCUUUGCACUCAUUGGUGAGUUUGGGUUUAACCCAAACUCACCAAUGAGCAA",
        height=120,
    )
with c2:
    mrna_name = st.text_input("mRNA Name", "ompF")
    mrna_seq = st.text_area(
        "mRNA Sequence (5'→3')",
        "AGACACATAAAGACACCAAACUCUCAUCA",
        height=120,
    )

c3, c4, c5 = st.columns([1.2, 1.2, 2])
with c3:
    demo_mode = st.checkbox("Relax constraints (demo mode)", True)
with c4:
    max_hits = st.slider("Max predictions", 1, 10, 5)
with c5:
    species = st.selectbox(
        "Organism / context (metadata only)",
        [
            "Escherichia coli",
            "Salmonella enterica",
            "Pseudomonas aeruginosa",
            "Staphylococcus aureus",
            "Mycobacterium tuberculosis",
            "Other / Unknown",
        ],
    )

# Extra context metadata (not yet sent to backend, but can be used later)
context = st.selectbox(
    "Biological context (metadata)",
    [
        "General",
        "Virulence",
        "Stress response",
        "Metabolism",
        "Biofilm / persistence",
        "Host-pathogen interaction",
        "Other",
    ],
)

export_csv = st.checkbox("Enable CSV export", True)

st.write("---")

# -----------------------------------------------------------------------------
# TABS: Home | Predict | Explain | Batch
# -----------------------------------------------------------------------------
tab_home, tab_predict, tab_explain, tab_batch, tab_validation = st.tabs(
    [" Home", " Predict", " Explain", "Batch", "Validation"]
)

# -----------------------------------------------------------------------------
# HOME TAB
# -----------------------------------------------------------------------------
with tab_home:
    st.subheader("Platform Overview")
    st.markdown(
        """
        This dashboard integrates:
        - IntaRNA thermodynamic predictions  
        - ML-based re-ranking of interactions  
        - GC/seed feature extraction  
        - Duplex ASCII + simple graphical visualization  
        - Gene-browser style arc plots  
        - CSV / Excel batch processing (CSV + XLSX)  
        - Scaffold for RNAcofold dot-bracket arc visualization  
        """
    )
    st.info(f"Backend expected at: **{API_BASE}** (/predict_intarna, /explain, /batch_predict)")


# -----------------------------------------------------------------------------
# PREDICT TAB
# -----------------------------------------------------------------------------
with tab_predict:
    st.subheader(" Run IntaRNA + ML prediction")

    if st.button(" Run Prediction"):
        if not srna_seq.strip() or not mrna_seq.strip():
            st.error("Please provide both sRNA and mRNA sequences.")
        else:
            # Metadata currently only used on UI side; not sent to backend yet
            payload = {
                "srna_name": srna_name,
                "mrna_name": mrna_name,
                "srna": srna_seq,
                "mrna": mrna_seq,
                "max_hits": max_hits,
                "demo": demo_mode,
            }
            result = api_post("/predict_intarna", payload)
            if not result:
                st.stop()

            # Global features
            st.markdown("#### Global Features")
            colA, colB, colC = st.columns(3)
            colA.metric("GC sRNA", f"{round(result.get('gc_content_srna', 0), 3)}")
            colB.metric("GC mRNA", f"{round(result.get('gc_content_mrna', 0), 3)}")
            colC.write("Seed features")
            seed_features = result.get("seed_features", {}) or {}
            colC.json(seed_features)

            interactions = result.get("interactions", [])
            if not interactions:
                st.warning("No interactions returned by IntaRNA.")
                st.stop()

            df = pd.DataFrame(interactions)

            st.markdown("####  Interactions Table")
            st.dataframe(df, use_container_width=True)

            # ΔG bar chart
            st.markdown("####  ΔG Binding Energies")
            make_dg_bar_chart(df)

            # Violin plot
            st.markdown("####  Violin Plot of ΔG and ML scores")
            make_violin_plot(df)

            # Gene-browser style arcs
            st.markdown("####  Gene-browser Arc View (mRNA)")
            make_arc_figure(df)

            # Sequence highlighting (seed + binding region) based on TOP interaction
            st.markdown("####  Sequence Highlighting (Top Interaction)")
            top = df.iloc[0]
            srna_bind_start = top.get("srna_start")
            srna_bind_end = top.get("srna_end")
            mrna_bind_start = top.get("target_start")
            mrna_bind_end = top.get("target_end")

            # Highlight seed differently for sRNA vs mRNA
            srna_seed_html = highlight_sequence_html(
                srna_seq,
                {
                    "has_seed": seed_features.get("has_seed"),
                    "length": seed_features.get("length"),
                    "srna_start": seed_features.get("srna_start"),
                },
                srna_bind_start,
                srna_bind_end,
            )
            mrna_seed_html = highlight_sequence_html(
                mrna_seq,
                {
                    "has_seed": seed_features.get("has_seed"),
                    "length": seed_features.get("length"),
                    "srna_start": seed_features.get("mrna_start"),  # reuse field
                },
                mrna_bind_start,
                mrna_bind_end,
            )

            cS, cM = st.columns(2)
            with cS:
                st.markdown("**sRNA (seed/binding highlighted)**", unsafe_allow_html=True)
                st.markdown(f"<code>{srna_seed_html}</code>", unsafe_allow_html=True)
            with cM:
                st.markdown("**mRNA (seed/binding highlighted)**", unsafe_allow_html=True)
                st.markdown(f"<code>{mrna_seed_html}</code>", unsafe_allow_html=True)

            # Per-interaction cards
            st.markdown("####  Per-Interaction Cards")

            for i, row in df.iterrows():
                with st.expander(
                    f"Rank {row['rank']}: ΔG={row['deltaG']} | ML={row.get('ml_score', 'NA')}",
                    expanded=(i == 0),
                ):
                    c1, c2, c3 = st.columns(3)
                    c1.metric("Target start", row.get("target_start"))
                    c2.metric("Target end", row.get("target_end"))
                    c3.metric("Hybrid length", row.get("hybrid_length"))

                    subseq_dp, hybrid_dp = extract_hybrid_from_raw_row(row.get("raw_row", ""))

                    st.markdown("**ASCII Duplex View**")
                    ascii_text = build_ascii_duplex(subseq_dp, hybrid_dp)
                    st.code(ascii_text, language="text")

                    st.markdown("**Graphical Duplex Plot (schematic)**")
                    make_graphical_duplex(subseq_dp, hybrid_dp)

                    st.markdown("**Raw CSV row**")
                    st.code(row.get("raw_row", ""), language="text")

            if export_csv:
                st.markdown("#### 📑 Export Results")
                st.download_button(
                    "Download predictions as CSV",
                    df.to_csv(index=False),
                    file_name="srna_mrna_predictions.csv",
                    mime="text/csv",
                )


# -----------------------------------------------------------------------------
# EXPLAIN TAB
# -----------------------------------------------------------------------------
with tab_explain:
    st.subheader(" Explain Top Interaction")

    if st.button(" Explain prediction for current input"):
        payload = {
            "srna_name": srna_name,
            "mrna_name": mrna_name,
            "srna": srna_seq,
            "mrna": mrna_seq,
            "max_hits": max_hits,
            "demo": demo_mode,
        }
        exp = api_post("/explain", payload)
        if not exp:
            st.stop()

        explain_block = exp.get("explain", {})
        top = exp.get("top_interaction")

        st.markdown("#### Explanation Summary")
        st.write(explain_block.get("why_rank1", ""))
        st.metric("Relative affinity", explain_block.get("relative_affinity", 0.0))

        st.markdown("#### Top Interaction Details")
        st.json(top)

        if top:
            subseq_dp, hybrid_dp = extract_hybrid_from_raw_row(top.get("raw_row", ""))
            st.markdown("#### ASCII Duplex for Top Interaction")
            ascii_text = build_ascii_duplex(subseq_dp, hybrid_dp)
            st.code(ascii_text, language="text")

            st.markdown("#### Graphical Duplex (Top Interaction)")
            make_graphical_duplex(subseq_dp, hybrid_dp)

        # --- RNAcofold / dot-bracket scaffold ---
        st.markdown("#### RNAcofold / dot-bracket arc view (if available)")
        cofold = exp.get("cofold")  # optional future backend field
        if cofold and isinstance(cofold, dict):
            seq_concat = cofold.get("sequence")
            struct = cofold.get("structure")
            st.text(f"Dot-bracket: {struct}")
            make_dotbracket_arc(seq_concat, struct)
        else:
            st.info(
                "Backend has not yet provided RNAcofold/dot-bracket information. "
                "Once the API returns a `cofold: {sequence, structure}` field in /explain, "
                "this panel will render arc plots from the dot-bracket."
            )

# -----------------------------------------------------------------------------
# BATCH TAB
# -----------------------------------------------------------------------------
with tab_batch:
    st.subheader(" Batch Mode")

    st.write(
        """
        Upload a CSV or Excel file with at least:

        - `srna_sequence` or `srna` (sRNA sequence, 5'→3')
        - `mrna_sequence` or `mrna` (mRNA sequence, 5'→3')
        - optionally: `srna_name`, `mrna_name`
        - optionally: `label_strength` (for validation / performance metrics)
        """
    )

    uploaded = st.file_uploader("Upload CSV or Excel", type=["csv", "xlsx", "xls"])

    if uploaded is not None:
        name = uploaded.name.lower()
        try:
            if name.endswith(".csv"):
                batch_df = pd.read_csv(uploaded)
            else:  # .xlsx or .xls
                batch_df = pd.read_excel(uploaded)
        except Exception as e:
            st.error(f"Failed to read uploaded file: {e}")
            batch_df = None


        if batch_df is not None:
            st.markdown("### Preview of uploaded data")
            st.dataframe(batch_df.head(), use_container_width=True)

            # We’ll try to be flexible with column names
            def get_col(df, options, required=False, desc=""):
                for c in options:
                    if c in df.columns:
                        return c
                if required:
                    st.error(f"Missing required column for {desc}. "
                             f"Looked for: {options}")
                    st.stop()
                return None

            col_srna_seq = get_col(
                batch_df,
                ["srna", "srna_sequence", "srna_seq"],
                required=True,
                desc="sRNA sequence",
            )
            col_mrna_seq = get_col(
                batch_df,
                ["mrna", "mrna_sequence", "mRNA_Target_Region_(Short, 5'-3')"],
                required=True,
                desc="mRNA sequence",
            )
            col_srna_name = get_col(
                batch_df,
                ["srna_name", "sRNA_Name", "sRNA"],
                required=False,
            )
            col_mrna_name = get_col(
                batch_df,
                ["mrna_name", "mRNA_Name", "Target_mRNA"],
                required=False,
            )
            col_label = get_col(
                batch_df,
                ["label_strength", "Label", "label"],
                required=False,
            )

            if st.button(" Run Batch Predictions"):
                requests_list = []
                for _, row in batch_df.iterrows():
                    srna_seq_row = str(row[col_srna_seq])
                    mrna_seq_row = str(row[col_mrna_seq])

                    srna_name_row = (
                        str(row[col_srna_name])
                        if col_srna_name
                        else "srna"
                    )
                    mrna_name_row = (
                        str(row[col_mrna_name])
                        if col_mrna_name
                        else "mrna"
                    )

                    requests_list.append(
                        {
                            "srna": srna_seq_row,
                            "mrna": mrna_seq_row,
                            "srna_name": srna_name_row,
                            "mrna_name": mrna_name_row,
                            "max_hits": max_hits,
                            "demo": demo_mode,
                        }
                    )

                out = api_post("/batch_predict", {"requests": requests_list})
                if not out:
                    st.stop()

                all_interactions = []
                for res in out.get("results", []):
                    for it in res.get("interactions", []):
                        all_interactions.append(it)

                if not all_interactions:
                    st.warning("No interactions returned for any pairs.")
                    st.stop()

                out_df = pd.DataFrame(all_interactions)
                st.markdown("### Combined batch interactions")
                st.dataframe(out_df, use_container_width=True)

                # ---- If label_strength exists, compute evaluation metrics ----
                if col_label is not None:
                    st.markdown("### 📈 Batch Performance vs. Ground Truth")

                    # Create a key for merging: (srna_name, mrna_name)
                    tmp_gold = batch_df.copy()
                    tmp_gold["__srna_name"] = (
                        tmp_gold[col_srna_name]
                        if col_srna_name
                        else "srna"
                    )
                    tmp_gold["__mrna_name"] = (
                        tmp_gold[col_mrna_name]
                        if col_mrna_name
                        else "mrna"
                    )

                    gold = tmp_gold[["__srna_name", "__mrna_name", col_label]].rename(
                        columns={
                            "__srna_name": "srna_name",
                            "__mrna_name": "mrna_name",
                            col_label: "label_strength",
                        }
                    )

                    merged = out_df.merge(
                        gold, on=["srna_name", "mrna_name"], how="left"
                    )

                    eval_df = merged.dropna(subset=["label_strength", "ml_score"]).copy()
                    st.write(f"Evaluable interactions: {eval_df.shape[0]}")

                    if eval_df.empty:
                        st.info(
                            "No overlap between predicted interactions and labeled rows "
                            "(check srna_name/mrna_name consistency)."
                        )
                    else:
                        # Basic metrics: correlation & RMSE
                        y_true = eval_df["label_strength"].astype(float)
                        y_score = eval_df["ml_score"].astype(float)

                        corr = y_true.corr(y_score)
                        rmse = np.sqrt(((y_true - y_score) ** 2).mean())
                        mae = np.abs(y_true - y_score).mean()

                        m1, m2, m3 = st.columns(3)
                        m1.metric("Pearson r", f"{corr:.3f}")
                        m2.metric("RMSE", f"{rmse:.3f}")
                        m3.metric("MAE", f"{mae:.3f}")

                        # Scatter plot label_strength vs ml_score
                        if HAS_PLOTLY:
                            import plotly.express as px

                            fig_scatter = px.scatter(
                                eval_df,
                                x="label_strength",
                                y="ml_score",
                                hover_data=["srna_name", "mrna_name", "deltaG"],
                                title="Label strength vs ML score (batch)",
                            )
                            st.plotly_chart(fig_scatter, use_container_width=True)
                        else:
                            st.info(
                                "Install `plotly` to see performance scatter plot "
                                "(pip install plotly)"
                            )

                        if export_csv:
                            st.download_button(
                                "Download merged predictions + labels as CSV",
                                merged.to_csv(index=False),
                                file_name="batch_predictions_with_labels.csv",
                                mime="text/csv",
                            )
                else:
                    st.info(
                        "No `label_strength` column detected, so only raw predictions "
                        "are shown. Add label_strength to your CSV for evaluation."
                    )

# -----------------------------------------------------------------------------
# VALIDATION TAB
# -----------------------------------------------------------------------------
with tab_validation:
    st.subheader(" Model Validation & Benchmarking")

    col_left, col_right = st.columns([1, 1])

    # -------------------------
    # LEFT: trigger validation & show raw JSON
    # -------------------------
    with col_left:
        st.markdown("#### Run batch validation on server")
        st.write(
            "This will run `notebooks/validation.py` on the VM, "
            "recompute metrics for RandomForest / XGBoost / IntaRNA baseline, "
            "and refresh the summary JSON exposed via `/metrics`."
        )

        refresh = st.button(" Run batch validation (recompute metrics)")

        try:
            if refresh:
                with st.spinner("Running validation on server..."):
                    resp = requests.get(
                        f"{API_BASE}/metrics", params={"refresh": "true"}
                    )
            else:
                resp = requests.get(f"{API_BASE}/metrics")
        except Exception as e:
            st.error(f"Backend metrics endpoint unavailable: {e}")
            st.stop()

    # Now we are outside col_left, but still inside tab_validation,
    # and `resp` is guaranteed to exist if we didn't hit the exception.

    if resp.ok:
        metrics = resp.json()
    else:
        st.error("Failed to fetch metrics from backend.")
        st.text(resp.text)
        st.stop()

    # Show dataset summary
    ds = metrics.get("dataset", {})
    st.markdown(
        f"**Samples used:** {ds.get('n_samples', '–')}  |  "
        f"Target: `{ds.get('target', 'label_strength')}`"
    )

    # Build comparison table
    rows = []
    for key, m in metrics.get("models", {}).items():
        if not m.get("available", False):
            continue
        reg = m.get("regression", {}) or {}
        binm = m.get("binary", {}) or {}
        rows.append(
            {
                "model_key": key,
                "model": m.get("name", key),
                "R²": reg.get("r2"),
                "RMSE": reg.get("rmse"),
                "MAE": reg.get("mae"),
                "ROC AUC": binm.get("roc_auc"),
                "Avg Precision": binm.get("average_precision"),
            }
        )

    if not rows:
        st.warning("No available models found in metrics JSON.")
        st.json(metrics)
        st.stop()

    df_metrics = pd.DataFrame(rows)

    # -------------------------
    # RIGHT: model table + plots
    # -------------------------
    with col_right:
        st.markdown("#### Model comparison table")
        st.dataframe(df_metrics.set_index("model"), use_container_width=True)

    # Simple bar plot comparing R² across models
    st.markdown("#### R² comparison")
    fig = go.Figure()
    fig.add_bar(
        x=df_metrics["model"],
        y=df_metrics["R²"],
    )
    fig.update_layout(
        xaxis_title="Model",
        yaxis_title="R²",
        template="plotly_white",
        height=400,
    )
    st.plotly_chart(fig, use_container_width=True)

    # Optional: second chart for ROC AUC
    st.markdown("#### ROC AUC comparison (where available)")
    fig2 = go.Figure()
    fig2.add_bar(
        x=df_metrics["model"],
        y=df_metrics["ROC AUC"],
    )
    fig2.update_layout(
        xaxis_title="Model",
        yaxis_title="ROC AUC",
        template="plotly_white",
        height=400,
    )
    st.plotly_chart(fig2, use_container_width=True)

    st.caption(
        f"Last metrics generation: {metrics.get('generated_at', 'n/a')}"
    )

    # --- Download per-interaction validation CSV ---
    st.markdown("#### ⬇️ Download per-interaction validation results")

    csv_path = Path("data") / "validation_with_predictions.csv"

    if csv_path.exists():
        try:
            df_val = pd.read_csv(csv_path)
            csv_bytes = df_val.to_csv(index=False).encode("utf-8")

            st.download_button(
                label=" Download validation_with_predictions.csv",
                data=csv_bytes,
                file_name="validation_with_predictions.csv",
                mime="text/csv",
            )

            with st.expander("Preview first rows"):
                st.dataframe(df_val.head(), use_container_width=True)
        except Exception as e:
            st.error(f"Failed to load validation_with_predictions.csv: {e}")
    else:
        st.info(
            "validation_with_predictions.csv not found yet. "
            "Click **' Run batch validation (recompute metrics)'** above to generate it."
        )

# -----------------------------------------------------------------------------
# FOOTER
# -----------------------------------------------------------------------------
st.write("---")
st.caption(
    "Backend: FastAPI + IntaRNA + ML model | Frontend: Streamlit | "
    "Extra: Batch CSV/XLSX, violin plots, duplex visualization, RNAcofold scaffold"
)

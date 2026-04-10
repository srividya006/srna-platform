import streamlit as st
import pandas as pd
import requests
import plotly.graph_objects as go

API_BASE = "http://localhost:8080"  # change to VM IP when hosted

# === Streamlit Global Config ===
st.set_page_config(
    page_title="sRNA–mRNA Interaction Dashboard",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- Styled Header ---
st.markdown("""
<style>
    .main-title {
        text-align: center;
        font-size: 32px;
        font-weight: bold;
        color: #2866cc;
    }
</style>
<div class="main-title">🧬 sRNA–mRNA Interaction Intelligence Platform</div>
""", unsafe_allow_html=True)


# ======================================================================================
# API Helpers
# ======================================================================================

def try_request(url, payload):
    try:
        res = requests.post(url, json=payload, timeout=60)
        res.raise_for_status()
        return res.json()
    except Exception as e:
        st.error(f"❌ API unreachable: {e}")
        return None


def run_predict(srna, mrna, srna_name, mrna_name, demo, max_hits):
    payload = {
        "srna": srna,
        "mrna": mrna,
        "srna_name": srna_name,
        "mrna_name": mrna_name,
        "max_hits": max_hits,
        "demo": demo,
    }
    return try_request(f"{API_BASE}/predict_intarna", payload)


def run_explain(srna, mrna, srna_name, mrna_name, demo, max_hits):
    payload = {
        "srna": srna,
        "mrna": mrna,
        "srna_name": srna_name,
        "mrna_name": mrna_name,
        "max_hits": max_hits,
        "demo": demo,
    }
    return try_request(f"{API_BASE}/explain", payload)


# ======================================================================================
# Sidebar Inputs
# ======================================================================================
st.sidebar.header("Sequence Input Panel")

srna_name = st.sidebar.text_input("sRNA Name", "MicF")
mrna_name = st.sidebar.text_input("mRNA Name", "ompF")

srna_seq = st.sidebar.text_area("sRNA Sequence (5'→3')", "", height=120)
mrna_seq = st.sidebar.text_area("mRNA Sequence (5'→3')", "", height=120)

demo_mode = st.sidebar.checkbox("Relax RNA constraints (demo mode)", True)
max_hits = st.sidebar.slider("Max Predictions", 1, 10, 5)

st.sidebar.write("---")
export_option = st.sidebar.checkbox("Enable CSV Export", True)

# ======================================================================================
# Tab Navigation
# ======================================================================================
tabs = st.tabs([" Home", "Predict", " Explain", " Batch Mode"])


# ======================================================================================
# HOME TAB
# ======================================================================================
with tabs[0]:
    st.subheader(" Platform Overview")
    st.write("""
    This dashboard integrates:
    - IntaRNA thermodynamic predictions  
    - ML-based ranking  
    - Sequence feature extraction  
    - Arc-based visualization  
    """)

    st.success("Backend connection expected at: **http://localhost:8080**")


# ======================================================================================
# PREDICT TAB
# ======================================================================================
with tabs[1]:
    st.header(" Run Prediction")

    if st.button(" Run IntaRNA + ML Ranking"):
        result = run_predict(srna_seq, mrna_seq, srna_name, mrna_name, demo_mode, max_hits)

        if result and "interactions" in result:
            interactions = result["interactions"]
            if interactions:
                df = pd.DataFrame(interactions)
                st.subheader(" Ranked Predictions Table")
                st.dataframe(df)

                # === ΔG Bar Chart ===
                fig = go.Figure()
                fig.add_trace(go.Bar(
                    x=df["rank"],
                    y=df["deltaG"],
                    text=df["ml_score"],
                    name="Energy",
                    marker_color="#1976D2"
                ))
                fig.update_layout(
                    title="ΔG Binding Energies",
                    xaxis_title="Rank",
                    yaxis_title="ΔG (kcal/mol)"
                )
                st.plotly_chart(fig, use_container_width=True)

                # === Gene Browser Interaction Arc Plot ===
                st.subheader(" Interaction Arc Plot")

                arc = go.Figure()
                for _, row in df.iterrows():
                    arc.add_shape(
                        type="line",
                        x0=row["srna_start"],
                        y0=0,
                        x1=row["mrna_start"],
                        y1=1,
                        line=dict(color="red", width=2)
                    )
                arc.update_layout(
                    xaxis_title="sRNA Position",
                    yaxis_title="mRNA Position",
                    title="Pairing Region Arc Map"
                )
                st.plotly_chart(arc, use_container_width=True)

                if export_option:
                    st.download_button(
                        "⬇ Export CSV",
                        df.to_csv(index=False),
                        file_name="interaction_predictions.csv",
                        mime="text/csv"
                    )
            else:
                st.warning("No interactions detected by the backend.")


# ======================================================================================
# EXPLAIN TAB
# ======================================================================================
with tabs[2]:
    st.header(" Explain Model Decision")

    if st.button(" Explain Top Interaction"):
        exp = run_explain(srna_seq, mrna_seq, srna_name, mrna_name, demo_mode, max_hits)

        if exp:
            st.subheader(" Why this interaction was ranked highest")
            st.write(exp["explain"]["why_rank1"])
            st.metric("Relative Affinity", exp["explain"]["relative_affinity"])
            st.json(exp["top_interaction"])


# ======================================================================================
# BATCH MODE TAB (Prototype UI)
# ======================================================================================
with tabs[3]:
    st.header(" Batch Processing")

    st.write(" Coming soon: upload CSV and get ranked outputs!")

    uploaded = st.file_uploader("Upload batch CSV (srna,mrna pairs)")

    if uploaded:
        df = pd.read_csv(uploaded)
        st.dataframe(df)
        st.info("Batch compute button will be implemented next phase.")


# ======================================================================================
# Footer
# ======================================================================================
st.write("---")
st.caption("Built with FastAPI + Streamlit + IntaRNA + ML ranking model")

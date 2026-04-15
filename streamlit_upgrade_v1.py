import streamlit as st
import requests
from requests.exceptions import ConnectionError
import pandas as pd
import altair as alt

# ======================
# CONFIG
# ======================
API_BASE = "http://srna-platform-backend.onrender.com"

st.set_page_config(
    page_title="sRNA-mRNA Interaction Platform",
    layout="wide",
    page_icon="🧬"
)

# ======================
# UI Styling
# ======================
CARD_STYLE = """
<div style="
    background-color:#f8f9fa;
    padding:18px;
    border-radius:12px;
    border:1px solid #ddd;
    box-shadow:1px 1px 4px rgba(0,0,0,0.1);
    margin-bottom:12px;">
{content}
</div>
"""

def card(text):
    st.markdown(CARD_STYLE.format(content=text), unsafe_allow_html=True)

# ======================
# API Helpers
# ======================
def run_predict(srna, mrna, srna_name, mrna_name, demo=False):
    payload = {
        "srna": srna,
        "mrna": mrna,
        "srna_name": srna_name,
        "mrna_name": mrna_name,
        "max_hits": 5,
        "demo": demo
    }
    try:
        res = requests.post(f"{API_BASE}/predict_intarna", json=payload, timeout=30)
        res.raise_for_status()
        return res.json()
    except ConnectionError:
        st.error(f"❌ Could not reach backend API at "
                 f"{API_BASE}. Check deployment!")
        return None
    except Exception as e:
        st.error(f"❌ API error: {e}")
        return None

def run_explain(srna, mrna, srna_name, mrna_name):
    payload = {
        "srna": srna,
        "mrna": mrna,
        "srna_name": srna_name,
        "mrna_name": mrna_name
    }
    res = requests.post(f"{API_BASE}/explain", json=payload)
    return res.json()

# ======================
# TABS
# ======================
tab_home, tab_predict, tab_explain, tab_batch = st.tabs(
    ["🏠 Home", "⚡ Predict", "🔍 Explain", "📦 Batch"]
)

# ======================
# HOME TAB
# ======================
with tab_home:
    st.title("🧬 sRNA-mRNA Cloud Prediction Platform")
    card("""
    🔹 Predict RNA-RNA interactions  
    🔹 Retrieve thermodynamics  
    🔹 ML-ranked interaction scoring  
    🔹 Explainability engine  
    🔹 CSV export  
    """)

    st.markdown("### What this UI supports")
    st.success("""
    ✔ IntaRNA backend  
    ✔ GC/seed feature computation  
    ✔ ML-based re-ranking  
    ✔ Visualization front-end  
    """)

# ======================
# PREDICT TAB
# ======================
with tab_predict:
    st.header("⚡ Run Prediction")

    col1, col2 = st.columns(2)
    with col1:
        srna_name = st.text_input("sRNA Name", "MicF")
        srna = st.text_area("sRNA Sequence", "")
    with col2:
        mrna_name = st.text_input("mRNA Name", "ompF")
        mrna = st.text_area("mRNA Sequence", "")

    demo = st.checkbox("Demo Mode (force interactions)", True)

    if st.button("Run IntaRNA Prediction"):
        result = run_predict(srna, mrna, srna_name, mrna_name, demo)
        
        if result and "interactions" in result:
            interactions = result["interactions"]
            if interactions:
                df = pd.DataFrame(interactions)

                st.subheader("📌 Ranked Predictions")
                st.table(df)

                # ========= ΔG Bar Plot ===========
                st.subheader("📊 ΔG Energy Chart")
                chart = alt.Chart(df).mark_bar().encode(
                    x="rank:O",
                    y="deltaG:Q",
                    color="rank:O",
                    tooltip=["rank", "deltaG", "ml_score"]
                )
                st.altair_chart(chart, use_container_width=True)

                # ========= CSV Export ============
                st.download_button(
                    label="📑 Download CSV",
                    data=df.to_csv(index=False),
                    file_name="prediction_results.csv",
                    mime="text/csv"
                )

# ======================
# EXPLAIN TAB
# ======================
with tab_explain:
    st.header("🔍 Explain Model Decisions")

    col1, col2 = st.columns(2)
    with col1:
        srna_name = st.text_input("Explain: sRNA name", "MicF")
        srna = st.text_area("Explain: sRNA seq")
    with col2:
        mrna_name = st.text_input("Explain: mRNA name", "ompF")
        mrna = st.text_area("Explain: mRNA seq")

    if st.button("Explain Top Interaction"):
        exp = run_explain(srna, mrna, srna_name, mrna_name)

        st.subheader("📌 Explanation Summary")
        st.write(exp.get("explain", {}))

        if exp.get("top_interaction"):
            t = exp["top_interaction"]
            st.info(
                f"Top ΔG = {t['deltaG']} | ML score = {t.get('ml_score')}"
            )

# ======================
# BATCH TAB
# ======================
with tab_batch:
    st.header("📦 Batch Processing")

    uploaded = st.file_uploader("Upload CSV with srna,mrna columns")

    if uploaded:
        batch_df = pd.read_csv(uploaded)
        st.write(batch_df.head())

        if st.button("Run Batch Prediction"):
            reqs = []
            for _, row in batch_df.iterrows():
                reqs.append({
                    "srna": row["srna"],
                    "mrna": row["mrna"],
                    "srna_name": row.get("srna_name", "srna"),
                    "mrna_name": row.get("mrna_name", "mrna"),
                    "max_hits": 3
                })

            res = requests.post(f"{API_BASE}/batch_predict", json={"requests": reqs})
            output = res.json()

            st.success(f"Processed {output['count']} pairs")

            all_results = []
            for r in output["results"]:
                for i in r.get("interactions", []):
                    all_results.append(i)

            if all_results:
                df_out = pd.DataFrame(all_results)
                st.write(df_out)

                st.download_button(
                    label="📑 Download Combined Output CSV",
                    data=df_out.to_csv(index=False),
                    file_name="batch_results.csv",
                    mime="text/csv"
                )

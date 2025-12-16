import streamlit as st
import requests
import pandas as pd

API_BASE = "http://localhost:8080"  # Change to your external IP when deployed


st.set_page_config(
    page_title="sRNA–mRNA Interaction Platform",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("🧬 sRNA–mRNA Interaction Prediction Platform")
st.write("Cloud-enabled RNA Interaction Intelligence System")

# Sidebar Inputs
st.sidebar.header("Sequence Inputs")

srna_name = st.sidebar.text_input("sRNA Name", "MicF")
mrna_name = st.sidebar.text_input("mRNA Name", "ompF")

srna_seq = st.sidebar.text_area(
    "sRNA Sequence (5'→3')",
    "AUGCAGGUUCAUCAACACCUUAGGAAACCCUGC",
    height=120,
)

mrna_seq = st.sidebar.text_area(
    "mRNA Target Sequence (5'→3')",
    "AGACACATAAAGACACCAAACUCUCAUCA",
    height=120,
)

run_demo = st.sidebar.checkbox("Relax RNA constraints (demo mode)", value=True)
max_hits = st.sidebar.slider("Max Predictions", 1, 10, 5)

st.sidebar.write("---")
run_rnafold = st.sidebar.checkbox("Compute secondary structure (RNAfold)", value=False)

col1, col2 = st.columns([1, 2])


# ===========================================================
# Run RNAfold
# ===========================================================
with col1:
    if run_rnafold and st.button("Run RNAfold"):
        resp = requests.post(
            API_BASE + "/rnafold",
            json={"name": srna_name, "sequence": srna_seq},
        )
        if resp.ok:
            fold = resp.json()
            st.subheader("RNAfold Result")
            st.write(f"Structure: `{fold['structure']}`")
            st.write(f"MFE: **{fold['mfe']} kcal/mol**")
        else:
            st.error(f"RNAfold Error: {resp.text}")


# ===========================================================
# Prediction Section
# ===========================================================
with col2:
    st.subheader("Predict Interaction (IntaRNA + ML Ranking)")

    if st.button("⚡ Predict Interaction"):
        payload = {
            "srna_name": srna_name,
            "mrna_name": mrna_name,
            "srna": srna_seq,
            "mrna": mrna_seq,
            "max_hits": max_hits,
            "demo": run_demo,
        }

        resp = requests.post(API_BASE + "/predict_intarna", json=payload)

        if not resp.ok:
            st.error("❌ API Error")
            st.write(resp.text)
        else:
            data = resp.json()

            st.write(f"**GC sRNA:** {round(data['gc_content_srna'],2)}")
            st.write(f"**GC mRNA:** {round(data['gc_content_mrna'],2)}")
            st.write(f"**Seed Features:** {data['seed_features']}")

            interactions = data["interactions"]
            if interactions:
                df = pd.DataFrame(interactions)
                st.dataframe(df)

                best = interactions[0]
                st.success(
                    f"Top Interacting Region: {best['srna_start']}→{best['srna_end']} "
                    f"(ΔG={best['deltaG']}, ML Score={best.get('ml_score', '-'):})"
                )
            else:
                st.warning("⚠ No interactions returned")


# ===========================================================
# Explainability
# ===========================================================
st.write("---")
st.header("🔍 Explain Top Interaction Prediction")

if st.button("Explain Model Decision"):
    payload = {
        "srna_name": srna_name,
        "mrna_name": mrna_name,
        "srna": srna_seq,
        "mrna": mrna_seq,
        "max_hits": max_hits,
        "demo": run_demo,
    }

    expl = requests.post(API_BASE + "/explain", json=payload)

    if expl.ok:
        exp = expl.json()

        st.write("### Explanation")

        st.write(exp["explain"]["why_rank1"])
        st.write(f"Relative affinity: **{exp['explain']['relative_affinity']}**")
        st.write("Top predicted interaction:")
        st.json(exp["top_interaction"])
    else:
        st.error("Explainability Failed")
        st.write(expl.text)


# ===========================================================
# Footer / Credits
# ===========================================================
st.write("---")
st.info("Built using FastAPI + IntaRNA + Streamlit + ML ranking")



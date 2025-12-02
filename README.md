<<<<<<< HEAD
# srna-platform
srna-mrna interaction prediction
=======
# Cloud-Native Hybrid Prediction Platform for sRNA–mRNA Interactions  
*With Special Focus on Mycobacterium tuberculosis*

---

## 📌 Project Overview

This project develops a **cloud-native prediction system** for identifying potential regulatory interactions between **bacterial small RNAs (sRNAs) and target mRNAs**.  
It integrates:

- Thermodynamic RNA structure prediction (RNAfold)
- RNA–RNA interaction modelling (IntaRNA)
- Feature extraction and scoring
- Web-based access through FastAPI

The initial biological application focuses on **Mycobacterium tuberculosis**, but the framework is **generalizable to other bacterial genomes**.

---

## 🚀 Features

✔ Predict sRNA–mRNA interactions via IntaRNA  
✔ Compute RNA structural profiles using RNAfold  
✔ Lightweight REST API interface  
✔ Cloud-ready deployment architecture  
✔ Modular pipelines for future ML integration  
✔ Extendable for other pathogens

---

## 📂 Project Structure

srna_project/
│── src/
│ ├── app.py # FastAPI backend
│ ├── utils.py # Helper functions (future extension)
│ └── models/ # ML models or scoring logic (future)
│
├── data/ # Input datasets (ignored from git)
├── results/ # Output predictions (ignored from git)
├── tests/ # Unit tests (future)
│
├── environment.yml # Conda environment definition
├── README.md # Project documentation
└── .gitignore
## 🔧 Installation

### 1️⃣ Clone repository

```bash
git clone <your_repo_url>
cd srna_project

2️⃣ Create conda environment
conda env create -f environment.yml
conda activate srna

3️⃣ Install required tools

Ensure IntaRNA and RNAfold are installed and callable:

RNAfold --version
IntaRNA --version

▶️ Usage
Start the API server:
uvicorn src.app:app --host 0.0.0.0 --port 8080

Example API call (curl)
curl -X POST "http://<server_ip>:8080/predict_intarna" \
     -H "Content-Type: application/json" \
     -d '{"srna": "AUGCUAUGCUA", "mrna": "GGGGAUACGAUAGCUAGCUA"}'

Roadmap / Milestones
-	Phase	Status
1. Environment Setup (conda, RNAfold, IntaRNA)	✔ Completed	
2. API Prototype for simple predictions	✔ Completed	
3. Add structured scoring + feature extraction	🔄 In Progress	
4. Add visualization layer	⏳	
5. Extend framework to other pathogens	⏳	
6. Cloud deployment (Cloud Run / Docker)	⏳

Test Execution

You may validate tool installation using:

echo -e ">test\nAUGCUAUGCUA" | RNAfold


and

IntaRNA -q test_srna.fa -t test_mrna.fa

Architecture Diagram
User → FastAPI → Prediction Engine → (RNAfold / IntaRNA) → Results JSON

Contributions welcome for:

Frontend UI

ML ranking layer

Other bacterial genomes

✔ License

Academic / Research Use Only.
Please cite appropriately if using in research output.
>>>>>>> 167edfe (Add README.md project documentation)

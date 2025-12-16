FINAL_README.md
sRNA–mRNA Interaction Prediction Platform

IntaRNA + Machine Learning + Dockerized Demo

1. Project Overview

This project implements an end-to-end platform for predicting bacterial sRNA–mRNA interactions by integrating:

IntaRNA (v3.4.1) for thermodynamic interaction prediction

Machine Learning re-ranking using extracted sequence and structural features

FastAPI for backend inference and batch processing

Streamlit for an interactive web-based frontend

Docker & Docker Compose for reproducible deployment and long-term demoability

The system supports single prediction, batch prediction, explainability, validation, and stress testing workflows.

2. System Architecture
+------------------+        HTTP        +----------------------+
|  Streamlit UI    |  <--------------> |  FastAPI Backend     |
|  (Port 8501)     |                  |  (Port 8080)         |
+------------------+                  +----------------------+
                                              |
                                              |
                                      +------------------+
                                      |   IntaRNA        |
                                      |   ML Re-ranker   |
                                      +------------------+


The frontend collects user input and visualizes results

The backend runs IntaRNA, extracts features, applies ML scoring, and returns ranked interactions

Containers communicate over an internal Docker network

3. Features Implemented
i. Prediction

Single sRNA–mRNA interaction prediction

Ranked interactions with ΔG and ML score

GC content computation

Seed detection heuristic

ii.Batch Mode

Accepts JSON batch input (/batch_predict)

Supports stress testing with 50–100+ interaction pairs

Measures latency and failure rates

iii.Explainability

Displays:

ΔG binding energy

Hybrid alignment

ML confidence score

Seed match presence

iv.Validation

Comparison against literature-validated datasets

Metrics produced:

Top-1 accuracy

Precision@k

Recall@k

Results stored as CSV and JSON for reproducibility

4. Technology Stack
Component	Technology
Interaction Tool	IntaRNA 3.4.1
ML Model	Random Forest
Backend API	FastAPI + Uvicorn
Frontend UI	Streamlit
Containerization	Docker, Docker Compose
Language	Python 3.10
5. Prerequisites

To run the demo locally:

Docker ≥ 20.x

Docker Compose (v2 or docker-compose)

Linux / macOS / Windows (via Docker Desktop)

No Conda or Python installation is required on the host system.

6. Running the Project (Demo Mode)

Step 1: Clone the Repository
git clone https://github.com/<your-username>/srna-platform.git
cd srna-platform

Step 2: Build and Start Containers
docker-compose up --build


This will:

Build the backend image with IntaRNA and ML dependencies

Start the FastAPI service

Start the Streamlit frontend

Step 3: Access the Services

Service     	URL
Streamlit UI	http://localhost:8501

FastAPI Docs	http://localhost:8080/docs

7. Demonstration Workflow
A. Interactive Prediction (UI)

Open Streamlit UI

Enter sRNA and mRNA sequences

Run prediction

Inspect:

Ranked interactions

ΔG values

ML confidence

Seed features

B. Batch Prediction (API)

Upload or POST batch JSON to /batch_predict

Review aggregated results and latency

C. Validation

Validation scripts and outputs are stored in:

data/
notebooks/


Includes comparison with known sRNA–target pairs

8. Dockerized Deployment Rationale

Ensures reproducibility across systems

Allows the project to be demonstrated after cloud free tiers expire

Eliminates dependency conflicts (Boost, ViennaRNA, IntaRNA)

Enables examiners or reviewers to run the system with a single command

9. Repository Structure
srna-platform/
├── Dockerfile
├── docker-compose.yml
├── README.md
├── FINAL_README.md
├── src/                # FastAPI backend
├── streamlit_app/      # Streamlit UI
├── data/               # Datasets, validation outputs
├── models/             # Trained ML model
├── notebooks/          # Validation and analysis
└── tools/              # Utility scripts

10. Stopping the Application
docker-compose down

11. Notes for Evaluators

The system is fully self-contained

No external services are required

Both UI and API can be tested independently

All results are reproducible using the provided datasets and containers

12. Author & Project Status

Author: Rama Srividya Chintalapalle

Project Type: Individual MCA Research Project

Status: Phase 5 complete (Testing, Validation & Optimization)

Deployment: Dockerized, demo-ready

13. License

This project is intended for academic and research demonstration purposes.

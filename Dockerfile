FROM mambaorg/micromamba:1.5.8

# Set workdir
WORKDIR /app

# Copy env first (better layer caching)
COPY environment.yml .

# Create conda environment
RUN micromamba create -y -n srna -f environment.yml && \
    micromamba clean --all --yes

# Activate environment
ENV PATH=/opt/conda/envs/srna/bin:$PATH
ENV CONDA_DEFAULT_ENV=srna

# Copy project files
COPY src/ src/
COPY data/ data/
COPY models/ models/
COPY requirements.txt .

# Expose API port
EXPOSE 8080

# Start FastAPI
CMD ["uvicorn", "src.app:app", "--host", "0.0.0.0", "--port", "8080"]

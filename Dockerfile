FROM mambaorg/micromamba:1.5.8

WORKDIR /app

# Copy environment first for caching
COPY environment.yml .

# Create environment with all necessary scientific and visualization tools
# We add plotly, pandas, and scikit-learn to be safe
RUN micromamba create -y -n srna -f environment.yml && \
    micromamba install -y -n srna -c conda-forge setuptools plotly pandas scikit-learn && \
    micromamba clean --all --yes

# Activate environment
ENV PATH=/opt/conda/envs/srna/bin:$PATH
ENV CONDA_DEFAULT_ENV=srna

# Copy the entire project structure
COPY . .

# Final pip safety check for any missed dependencies in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt || true
RUN pip install --no-cache-dir plotly

# Expose ports for both services
EXPOSE PORT=10000

# The CMD is usually handled by docker-compose, so we leave it generic here
CMD ["uvicorn", "src.app:app", "--host", "0.0.0.0", "--port", "10000"]



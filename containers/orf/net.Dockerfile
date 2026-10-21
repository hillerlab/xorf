# ==============================================================================
# net.Dockerfile - Container for 'orf nets' subcommand
# ==============================================================================

# ---------- Build Stage ----------
# bookworm like the runtime stage: bullseye is past its end of life and its security
# archive answers 404 for the packages installed below
FROM rust:1.93-slim-bookworm AS builder

# Build dependencies for Rust project
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    libssl-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

COPY modules/orf/Cargo.toml modules/orf/Cargo.lock ./
COPY modules/orf/src/ ./src/

RUN cargo build --release && \
    strip target/release/orf

# ---------- Python Dependencies Stage ----------
FROM python:3.10-slim-bookworm AS py_builder

# Build tools for Python wheels that need compilation
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Copy only pyproject.toml first for better layer caching
COPY modules/nets/pyproject.toml ./

# Install dependencies in a virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install Python dependencies
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir \
    torch==1.12.1 \
    numpy==1.23.5 \
    pandas==2.0.3 \
    tqdm==4.62.3 \
    transformers==4.36.0 \
    huggingface_hub>=0.20.3 \
    biopython>=1.78 \
    scikit-learn>=0.24.0

# ---------- Runtime Stage ----------
FROM python:3.10-slim-bookworm

# Install minimal runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    wget \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/target/release/orf /usr/local/bin/orf
COPY --from=py_builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy the tools
COPY modules/nets/netstart2 /opt/netstart2
COPY modules/nets/TRANSAID /opt/transaid

# Downloading models for netstart2
RUN mkdir -p /opt/netstart2/src/model_config/models \
    && wget -O /opt/netstart2/src/model_config/models/netstart_model1.pth https://huggingface.co/linesandvad/netstart2_models/resolve/main/netstart_model1.pth?download=true \
    && wget -O /opt/netstart2/src/model_config/models/netstart_model2.pth https://huggingface.co/linesandvad/netstart2_models/resolve/main/netstart_model2.pth?download=true \
    && wget -O /opt/netstart2/src/model_config/models/netstart_model3.pth https://huggingface.co/linesandvad/netstart2_models/resolve/main/netstart_model3.pth?download=true \
    && wget -O /opt/netstart2/src/model_config/models/netstart_model4.pth https://huggingface.co/linesandvad/netstart2_models/resolve/main/netstart_model4.pth?download=true

# Download and persist tokenizer assets at build time to avoid runtime Hub requests.
RUN mkdir -p /opt/netstart2/src/model_config/pretrained_models/tokenizers/facebook_esm2_t6_8M_UR50D \
    && python -c "from transformers import AutoTokenizer; AutoTokenizer.from_pretrained('facebook/esm2_t6_8M_UR50D', do_lower_case=False).save_pretrained('/opt/netstart2/src/model_config/pretrained_models/tokenizers/facebook_esm2_t6_8M_UR50D')"

# Make scripts executable and expose convenience shims
RUN chmod +x /opt/netstart2/netstart2.py && \
    chmod +x /opt/transaid/transaid.py && \
    chmod +x /usr/local/bin/orf && \
    ln -s /opt/netstart2/netstart2.py /usr/local/bin/netstart2 && \
    ln -s /opt/transaid/transaid.py /usr/local/bin/transaid

# Add to PATH so 'predict' command works
ENV PATH="/opt/netstart2:$PATH"
ENV PATH="/opt/transaid:$PATH"
ENV HF_HUB_OFFLINE=1
ENV TRANSFORMERS_OFFLINE=1

# Set up non-root user
RUN useradd -m -u 1000 netuser && \
    chown -R netuser:netuser /opt/netstart2/src/model_config/models && \
    chown -R netuser:netuser /opt/netstart2/src/model_config/pretrained_models/tokenizers && \
    chown -R netuser:netuser /opt/netstart2/src/model_config/hyperparameters && \
    chown -R netuser:netuser /opt/transaid/model

USER netuser

WORKDIR /data

# Test that dependencies are available
RUN python -c "import pandas; import numpy; import sklearn; from transformers import AutoTokenizer, AutoModel; from huggingface_hub import hf_hub_download; print('All dependencies OK')"

# Create test fasta file to verify the tools are working
RUN echo ">test\nACGTCAATGAGATGACACAGTAGATGATGG" > /data/test.fasta

# Verify the script is accessible
RUN netstart2 --help
RUN transaid --help

# Test tools are working
RUN netstart2 -in /data/test.fasta -compute_device cpu -o chordata -out /data/test_netstart2_out
RUN transaid --input /data/test.fasta --gpu -1 --output /data/test_transaid_out

# Check that the output files are not empty
RUN test -s /data/test_netstart2_out.csv
RUN test -s /data/test_transaid_out.csv

# Check versioning
RUN netstart2 --version
RUN transaid --version

# Clean up test fasta file
RUN rm -rf /data/test_transaid_out.csv /data/test_netstart2_out.csv

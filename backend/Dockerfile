FROM node:20-alpine AS frontend-build
WORKDIR /frontend
COPY frontend/package*.json ./
RUN npm install
COPY frontend/ ./
RUN npm run build

FROM mambaorg/micromamba:1.5.8-bookworm-slim

USER root

# Install system dependencies needed for graphics/weasyprint and general tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libpango-1.0-0 \
    libpangoft2-1.0-0 \
    libharfbuzz0b \
    libffi-dev \
    libjpeg-dev \
    libopenjp2-7 \
    libgomp1 \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

WORKDIR /workspace

# Copy environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /workspace/environment.yml

# Create environment from environment.yml
RUN micromamba create -y -n docklab -f /workspace/environment.yml && \
    micromamba clean --all --yes

# Set environment path in shell
ENV PATH="/opt/conda/envs/docklab/bin:$PATH"
ENV CONDA_DEFAULT_ENV="docklab"

# Copy source code
COPY --chown=$MAMBA_USER:$MAMBA_USER . /workspace
COPY --from=frontend-build --chown=$MAMBA_USER:$MAMBA_USER /frontend/dist /workspace/frontend/dist

# Default entrypoint activates micromamba env
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["sh", "-c", "uvicorn backend.app.main:app --host 0.0.0.0 --port ${PORT:-8000}"]


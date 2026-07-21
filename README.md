# Genomic Pipeline

This repository contains a Nextflow DSL2 pipeline for genomic analysis, including FastQC, alignment, variant calling, annotation, and MultiQC reporting.

## Requirements

- Nextflow installed and available on your PATH
- Docker installed and running
- Optional: a Python virtual environment if you use one for local tooling

## Recommended setup

If you are using a Python virtual environment, activate it before running Nextflow:

```bash
# macOS / Linux
source /path/to/venv/bin/activate

# Windows (PowerShell)
# .\path\to\venv\Scripts\Activate.ps1
```

If you are not using a Python virtualenv, you can skip this step.

## Run the pipeline

From the repository root:

```bash
nextflow run main.nf -resume -with-docker
```

This command resumes any previous run and executes pipeline processes inside Docker containers.

## Files

- `main.nf` - Nextflow pipeline definition
- `nextflow.config` - pipeline configuration
- `.gitignore` - files and folders ignored by Git

## Notes

- Input data is expected under `data/`
- Results are written to `results/`
- Nextflow working files are stored in `work/`
- If you need to restart cleanly, remove the `work/` directory or use the Nextflow `-resume` option carefully

# AlphaOmega Genomic Analysis Pipeline

AlphaOmega is a local, privacy-first pipeline that takes a raw human genome (e.g., a massive WGS BAM file) and generates interactive clinical-grade reports. It annotates your variants against SNPedia, ClinVar, and PharmCAT without ever sending your genetic data to the cloud.

---

## 🚀 Getting Started

Ensure you have Docker and Docker Compose installed.

### 1. Build the Environment
Always run this when code or dependencies in `requirements.txt` are updated:
```bash
docker compose build
```

### 2. Initial Setup
Download the necessary knowledge bases (SNPedia, ClinVar, basic reference genomes):
```bash
docker compose run --rm alphaomega setup
```
*Note: This will populate the `./db` and `./reference` folders on your host machine.*

---

## 🧬 Core Commands

### Analyze a Genome
This command performs variant calling (running `bcftools mpileup`), annotates the called variants with rsIDs and clinical data, and generates a standalone interactive HTML report.

```bash
docker compose run --rm alphaomega analyze --bam data/your_genome.bam --profile "YourName"
```

**Options:**
- `--force`: If your profile already has a cached raw VCF, AlphaOmega will skip the lengthy variant calling step and just re-run annotation and reporting. Add `--force` if you want to explicitly overwrite the cache and restart variant calling from scratch.

### Compare & Find Filiation
Compare multiple analyzed profiles to discover shared genetic traits and evaluate relatedness/filiation.

```bash
docker compose run --rm alphaomega compare --profiles "Profile1,Profile2" --relatedness
```
*What this does:*
- **Shared vs Unique Variants:** Lists variants shared between the designated profiles.
- **Relatedness Analysis:** Calculates IBS0/IBS1/IBS2 ratios, checks Mendelian errors, and predicts the exact biological relationship (e.g., Parent-Child, Full Siblings) along with a confidence score.

---

## 🛠️ Handling Corner Scenarios & Troubleshooting

During development, several edge cases were identified when processing 100GB+ Whole Genome Sequencing files. Here is how to handle them:

### 1. Chromosome Naming Mismatches
**The Problem:** Your BAM file uses Ensembl/UCSC standard naming (e.g., `chr1`, `chr2`), but the reference index expects raw numbers (e.g., `1`, `2`). This causes an immediate pipeline crash because regions cannot be aligned.
**The Fix:** Use `samtools reheader` to strip the `chr` prefix and create a new BAM file.
```bash
samtools view -H data/original.bam | sed -e 's/SN:chr/SN:/' | samtools reheader - data/original.bam > data/fixed.bam
samtools index data/fixed.bam
```
*You can then throw away `original.bam` to save space and run the pipeline on `fixed.bam`.*

### 2. Missing rsIDs (`.`) In VCF Output
**The Problem:** `bcftools mpileup` outputs variant coordinates but does not stamp them with standard `rsIDs` (it outputs `.` in the ID column). The downstream annotator expects rsIDs and will produce an empty report with 0 findings.
**The Fix:** AlphaOmega has an automated injection fallback. It downloads the **build-matched ClinVar VCF** (e.g., GRCh37), extracts the `rs`-prefixed IDs from the ClinVar `INFO/RS` field, creates a specialized `POS->rsID` index, and automatically annotates your raw VCF by matching position and chromosome before passing it to the report generator.

### 3. Exhausted Docker Disk Space (WSL2)
**The Problem:** Docker running on Windows (WSL2) operates within a virtual disk (ext4.vhdx). Processing massive BAMs and generating uncompressed temporary VCF files can quickly bloat this virtual disk. When the disk limit is hit, the Docker daemon crashes with `no space left on device` or `File write failed`.
**The Fix:** 
- Keep your data on the Windows host and selectively mount it (which `.docker-compose.yml` is already configured to do).
- Periodically clear stopped containers and dangling images: `docker system prune`
- If you use Docker Desktop on Windows, you may need to use `diskpart` to compact the WSL2 `ext4.vhdx` to reclaim unallocated space after deleting large files.

### 4. Variant Calling Timeouts
**The Problem:** `bcftools mpileup` can take upwards of **4 to 8 hours** on a 110GB+ BAM file. The Python `subprocess.run` default timeouts will prematurely kill the process, destroying the output cache.
**The Fix:** The pipeline is tuned to allow long execution limits (`timeout=28800` / 8 hours for VCF execution) and features a background threading monitor that prints progress logs (e.g., "Still calling variants... 119 mins elapsed, 66.0 MB written") every 60 seconds without risking pipe deadlocks on `stderr`.

---

## 📊 The Interactive HTML Report
The pipeline outputs an interactive report at `profiles/<ProfileName>/report.html`.
Because the frontend is powered by vanilla JavaScript directly inside the HTML file, **no backend server is required**. You can double-click the file to open it in your browser offline.

The report includes a **Filter Bar** to sift through millions of variants:
- **Domain:** Filter variants by physiological impact (e.g., Metabolism, Disease Risk, Pharmacogenomics).
- **Significance:** Filter by clinical annotations (Pathogenic, Benign, Uncertain).
- **Notable Findings Only:** Instantly toggle the view to only show findings with a high impact (SNPedia Magnitude &ge; 2 or Pathogenic clinical significance).

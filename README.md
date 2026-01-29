# genoSubSampler

A tool for subsampling paired-end FASTQ files to meet file size constraints. Useful for uploading sequencing data to platforms with file size limits.

An interactive version is available at [DrDx.Me](http://www.drdx.me/).

## Features

- **Automatic subsampling** - Reduces FASTQ files to meet specified size limits while maintaining paired-end integrity
- **Better compression** - Uses 7zip instead of pigz/gzip for ~5% better compression ratios
- **M.tb filtering** - Optional filtering for *Mycobacterium tuberculosis* aligned sequences using bowtie2
- **Detailed summary** - Generates a CSV report with read counts and retention statistics
- **Reproducible** - Uses fixed seed for deterministic subsampling

## Installation

### Quick Install

```bash
git clone https://github.com/SemiQuant/genoSubSampler.git
cd genoSubSampler
chmod +x install.sh subsample.sh
./install.sh
```

### Install with M.tb Filtering Support

```bash
./install.sh --with-filter
```

This installs additional dependencies (bowtie2, samtools, bedtools, htslib) and builds a bowtie2 index from the included *M. tuberculosis* H37Rv reference genome.

### Manual Installation

#### Core Dependencies

| Tool | Purpose |
|------|---------|
| [seqtk](https://github.com/lh3/seqtk) | Fast FASTQ/FASTA manipulation |
| [7zip](https://www.7-zip.org/) | Compression (better ratios than gzip) |
| bc | Arithmetic calculations |

#### Additional Dependencies (for M.tb filtering)

| Tool | Purpose |
|------|---------|
| [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) | Sequence alignment |
| [samtools](http://www.htslib.org/) | BAM file processing |
| [bedtools](https://bedtools.readthedocs.io/) | BAM to FASTQ conversion |
| [htslib](http://www.htslib.org/) | bgzip compression |

#### Using micromamba/conda

```bash
# Core dependencies
micromamba install -c bioconda -c conda-forge seqtk p7zip bc

# With M.tb filtering
micromamba install -c bioconda -c conda-forge seqtk p7zip bc bowtie2 samtools bedtools htslib
```

#### Using Homebrew (macOS)

```bash
brew install seqtk p7zip bc
# For filtering: brew install bowtie2 samtools bedtools htslib
```

#### Using apt (Debian/Ubuntu)

```bash
sudo apt-get install seqtk p7zip-full bc
# For filtering: sudo apt-get install bowtie2 samtools bedtools tabix
```

## Usage

### Basic Usage

```bash
./subsample.sh -d /path/to/fastq/files
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-d, --dir` | Directory containing FASTQ files (required) | - |
| `-m, --max-size` | Maximum file size in MB | 98 |
| `-r, --reads` | Initial number of reads for subsampling | 2,000,000 |
| `-t, --threads` | Number of threads for bowtie2 | 4 |
| `-f, --filter` | Filter for M.tb sequences | disabled |
| `--ref` | Bowtie2 reference index path | (built-in) |
| `-h, --help` | Show help message | - |
| `-v, --version` | Show version | - |

### Examples

#### Subsample files to default 98MB limit

```bash
./subsample.sh -d /data/sequencing/run001
```

#### Subsample with custom size limit (50MB)

```bash
./subsample.sh -d /data/sequencing/run001 -m 50
```

#### Filter for M.tb sequences and subsample

```bash
./subsample.sh -d /data/sequencing/run001 -f -t 8
```

#### Use custom reference for filtering

```bash
./subsample.sh -d /data/sequencing/run001 -f --ref /path/to/bowtie2_index
```

## Input Requirements

- Paired-end FASTQ files with naming convention: `*_R1_001.fastq.gz` and `*_R2_001.fastq.gz`
- Files must be gzip compressed (`.fastq.gz`)

## Output

### Subsampled Files

- Original files are replaced with subsampled versions
- Naming convention: `*_R1_001_subsampled.fastq.gz`

### Summary Report

A `subsampling_summary.txt` file is generated with the following columns:

| Column | Description |
|--------|-------------|
| Sample_ID | Sample identifier |
| Original_R1_Reads | Original R1 read count |
| Original_R2_Reads | Original R2 read count |
| Subsampled_R1_Reads | Final R1 read count |
| Subsampled_R2_Reads | Final R2 read count |
| R1_Retention_% | Percentage of R1 reads retained |
| R2_Retention_% | Percentage of R2 reads retained |
| File_Size_R1_MB | Final R1 file size (MB) |
| File_Size_R2_MB | Final R2 file size (MB) |

Example output:

```
Sample_ID          Original_R1  Original_R2  Subsampled_R1  Subsampled_R2  R1_%    R2_%    R1_MB  R2_MB
Sample001          5000000      5000000      1850000        1850000        37.00   37.00   95     96
Sample002          2000000      2000000      2000000        2000000        100.00  100.00  45     46
```

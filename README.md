# ATAC-seq Analysis Pipeline
A Snakemake workflow for processing ATAC-seq data from raw reads to peak calling, with extensive QC metrics and analysis. This pipeline handles paired-end Illumina ATAC-seq data.

## Overview

This pipeline performs the following steps:

1. **Quality control and preprocessing**
   - FastQC for raw read QC
   - Fastp for adapter trimming and quality filtering

2. **Alignment**
   - Alignment with HISAT2
   - Filtering, sorting and indexing with Samtools

3. **Post-processing**
   - Deduplication with Picard
   - Blacklist region filtering with BEDTools

4. **Peak calling and analysis**
   - Peak calling with MACS2
   - IDR analysis for replicate consistency

5. **QC metrics and visualization**
   - Fragment size distribution
   - TSS enrichment
   - Fingerprint analysis
   - GC bias assessment
   - PCA/correlation analysis
   - Library complexity metrics (NRF, PBC1, PBC2)
   - Fraction of Reads in Peaks (FRiP)
   - Mitochondrial content analysis

## Dependencies

The pipeline uses conda environments to manage dependencies. The main requirements are:

- Snakemake ≥7.0.0
- Python ≥3.8
- FastQC
- Fastp
- HISAT2
- Samtools
- Picard
- BEDTools
- MACS2
- DeepTools
- IDR

## Installation

```bash
# Clone this repository
git clone https://github.com/yourusername/atac-seq-pipeline.git
cd atac-seq-pipeline

# Create conda environment with Snakemake
conda env create -f environment.yaml
conda activate atac_seq
```

## Usage

### 1. Prepare your data

Place raw FASTQ files in the `data/` directory following this naming convention:
```
data/{sample}_R1_001.fastq.gz
data/{sample}_R2_001.fastq.gz
```

### 2. Create sample sheet

Create a CSV file (e.g., `samples.csv`) with sample information:
```csv
sample_id,group
sample1,control
sample2,control
sample3,treatment
sample4,treatment
```

### 3. Configure the pipeline

Edit `ref/config.yaml` to define parameters and reference files:
```yaml
# Sample information
samples_table: "samples.csv"  # CSV with sample information

# Reference data
hisat2_index: "path/to/index"
blacklist: "path/to/blacklist.bed"
gtf_file: "path/to/gencode.v36.annotation.gtf"
genome_2bit: "path/to/hg38.2bit"
```

### 4. Run the pipeline

```bash
# Dry run to check workflow
snakemake -n

# Run locally
snakemake --use-conda --cores 24

# Run on a cluster (SLURM example)
snakemake --use-conda --cluster "sbatch -p {params.partition} -c {threads} -t {params.time}" --jobs 100
```
flowchart TD
    %% Main workflow steps
    rawdata[/"Raw FASTQ Files
    {sample}_R1_001.fastq.gz
    {sample}_R2_001.fastq.gz"/]
    
    %% Quality Control and Preprocessing
    subgraph "Quality Control & Preprocessing"
        fastqc["FastQC
        Quality Assessment"]
        fastp["Fastp
        Adapter Trimming
        Quality Filtering"]
    end
    
    %% Alignment and Filtering
    subgraph "Alignment & Filtering"
        hisat2["HISAT2 Alignment
        No Spliced Alignment
        --no-discordant --no-mixed"]
        samtools["Samtools
        Filter Proper Pairs
        Sort & Index"]
        picarddedup["Picard MarkDuplicates
        Remove PCR Duplicates"]
        blfilter["Bedtools
        Blacklist Region Filtering"]
    end
    
    %% Peak Calling
    subgraph "Peak Calling"
        macs2["MACS2
        Peak Calling
        --nomodel --shift -75"]
        idr["IDR Analysis
        Replicate Consistency"]
    end
    
    %% Signal Generation
    subgraph "Signal Generation"
        bigwig["DeepTools bamCoverage
        BigWig Generation"]
        bedgraph["DeepTools bamCoverage
        Bedgraph Generation"]
    end
    
    %% QC Metrics
    subgraph "QC Metrics"
        frip["FRiP Score
        Fraction of Reads in Peaks"]
        libcomplex["Library Complexity
        NRF, PBC1, PBC2"]
        fragsize["Fragment Size Analysis
        bamPEFragmentSize"]
        fingerprint["Fingerprint Analysis
        plotFingerprint"]
        gcbias["GC Bias Assessment
        computeGCBias"]
        mtcontent["MT Content Analysis
        Mitochondrial Contamination"]
        tss["TSS Enrichment Analysis
        computeMatrix & plotHeatmap"]
        correlation["Sample Correlation
        multiBamSummary & plotCorrelation"]
        pca["Principal Component Analysis
        plotPCA"]
        blstats["Blacklist Filtering Statistics"]
    end
    
    %% Connections
    rawdata --> fastqc
    rawdata --> fastp
    fastp --> hisat2
    hisat2 --> samtools
    samtools --> picarddedup
    picarddedup --> blfilter
    blfilter --> macs2
    blfilter --> bigwig
    blfilter --> bedgraph
    
    %% QC Connections
    macs2 --> idr
    blfilter --> frip
    macs2 --> frip
    hisat2 --> libcomplex
    blfilter --> fragsize
    blfilter --> fingerprint
    blfilter --> gcbias
    blfilter --> mtcontent
    bigwig --> tss
    blfilter --> correlation
    correlation --> pca
    picarddedup --> blstats
    blfilter --> blstats
    
    %% Style
    classDef processing fill:#e1f5fe,stroke:#01579b,stroke-width:2px,color:#01579b
    classDef qc fill:#f3e5f5,stroke:#4a148c,stroke-width:2px,color:#4a148c
    classDef data fill:#fffde7,stroke:#f57f17,stroke-width:2px,color:#f57f17
    
    class fastqc,fastp,hisat2,samtools,picarddedup,blfilter,macs2,bigwig,bedgraph processing
    class frip,libcomplex,fragsize,fingerprint,gcbias,mtcontent,tss,correlation,pca,blstats,idr qc
    class rawdata data

## Output Structure

```
results/
├── fastqc/                    # FastQC reports
├── fastp/                     # Trimmed reads and reports
├── aligned/                   # HISAT2 alignment files
├── filtered/                  # Filtered BAM files 
├── dedup/                     # Deduplicated BAM files
├── blacklist_filtered/        # Blacklist-filtered BAM files
├── peaks/                     # MACS2 peak calls
├── bigwig/                    # BigWig files for visualization
├── bedgraph/                  # Bedgraph files
├── idr/                       # IDR analysis results
├── deeptools/                 # DeepTools analysis results
│   ├── fragmentSize.png       # Fragment size distribution
│   ├── ATACseq_fingerprint.pdf # Fingerprint plot
│   ├── gc_content.pdf         # GC bias assessment
│   ├── Heatmap_TSS.pdf        # TSS heatmap
│   └── Profile_TSS.pdf        # TSS profile plot
├── FRiP/                      # Fraction of reads in peaks
├── library_complexity/        # Library complexity metrics
└── qc/                        # QC summary reports
```

## Quality Control Metrics

| Metric | Description | Good | Acceptable | Poor |
|--------|-------------|------|------------|------|
| TSS Enrichment | Signal enrichment at transcription start sites | >10 | 5-10 | <5 |
| FRiP | Fraction of reads in peaks | >0.3 | 0.2-0.3 | <0.2 |
| NRF | Non-redundant fraction | >0.9 | 0.8-0.9 | <0.8 |
| PBC1 | PCR bottleneck coefficient 1 | >0.9 | 0.7-0.9 | <0.7 |
| PBC2 | PCR bottleneck coefficient 2 | >3 | 1-3 | <1 |
| MT Content | Mitochondrial read % | <20% | 20-30% | >30% |
| Fragment Size | Nucleosomal pattern | Clear pattern | Weak pattern | No pattern |
| IDR | Irreproducible Discovery Rate | >80% | 70-80% | <70% |

## Troubleshooting

### Common Issues

1. **Low alignment rate**: Check adapter trimming and reference genome compatibility
2. **High duplication rate**: Indicates low library complexity; optimize input material
3. **Poor TSS enrichment**: Issues with sample quality or tagmentation
4. **High mitochondrial content**: Consider nuclear isolation protocols
5. **Failed IDR**: Low reproducibility between replicates; check protocol consistency
6. **GC bias**: Consider GC bias correction methods for downstream analysis
7. **Peak number too low/high**: Adjust MACS2 parameters (q-value/p-value cutoffs)
8. **Blacklist filtering removes too many reads**: Check blacklist regions compatibility

### Logs

Error logs for each step are stored in the `logs/` directory, organized by rule name:
```
logs/
├── fastqc/
├── fastp/
├── hisat2/
├── samtools/
├── dedup/
├── blacklist_filter/
├── macs2/
└── ...
```

## Citation

If you use this pipeline in your research, please cite:

```
Your Name et al. (20XX). ATAC-seq analysis pipeline for comprehensive chromatin accessibility profiling.
```

📄 License
This project is licensed under the MIT License - see the LICENSE file for details.

📧 Contact
For questions or feedback, please open an issue on the GitHub repository or contact the author.

Last updated: April 30, 2025
Created by: gynecoloji

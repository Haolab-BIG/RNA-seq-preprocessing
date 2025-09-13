# RNA-seq Pipeline

This unified RNA-seq pipeline processes raw FASTQ files through to differential-expression results, using Singularity for reproducibility and supporting batch analysis of multiple samples.

## Workflow Diagram

<img width="1750" height="447" alt="pipeline" src="https://github.com/user-attachments/assets/a9d49331-4be1-4093-b955-a97db6be8083" />

## Features

  - **Single Command Execution**: Executes the entire workflow—from FASTQ input through sequential per-sample processing to final differential-expression analysis—with a single command.
  - **Two Modes**: Supports both alignment-based (`align`) and alignment-free (`quant`) modes.
      - **Alignment-based (`align`) mode**: This is the traditional and most comprehensive approach. It maps each sequencing read to a reference genome, providing precise information about where each read came from.

          - **When to use**: Offers the most detailed, base-level view of the transcriptome, enabling downstream analyses beyond expression quantification. It generates BAM and BIGWIG files that can be visualized in genome browsers like IGV.

      - **Alignment-free (`quant`) mode**: A lightweight approach that bypasses full genome alignment. It directly estimates transcript abundance by decomposing reads into k-mers and matching them to a prebuilt transcriptome index.

          - **When to use**: Choose this mode if your primary goal is transcript-level or gene-level differential expression analysis. It is significantly faster and requires less computational resources.
  - **Reproducible**: All software is encapsulated within a Singularity container (`RNA.sif`).

## Requirements

1.  **Recommended System Configuration**:

      * 8-core CPU
      * 80 GB RAM

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
			libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```

3.  **Pipeline Files**:

      * `run_pipeline.sh`
      * `RNA.sif` (The Singularity container)

4.  **Reference Data**: A directory containing all necessary reference files (e.g., STAR indices, Salmon index, GTF annotation, etc.).

## Setup

### 1\. Prepare the Sample Sheet

This is the most critical input file. Create a CSV file named `samplesheet.csv`.

  - `sample`: A unique identifier for the sample (e.g., `Control_Rep1`). This name will be used for output subdirectories.
  - `condition`: The experimental group for the sample (e.g., `Control`, `Treated`). This is used for the DESeq2 design.
  - `fastq1_path`: The **absolute path** to the Read 1 FASTQ file.
  - `fastq2_path`: The **absolute path** to the Read 2 FASTQ file.

**Note on Sequencing Type:**
This pipeline supports both paired-end (PE) and single-end (SE) sequencing data. The example below shows the format for paired-end data. If you have single-end data, simply remove the `fastq2_path` column.

**Example `samplesheet.csv`:**

```csv
sample,condition,fastq1_path,fastq2_path
Control_Rep1,Control,/path/to/data/Control_Rep1_R1.fastq.gz,/path/to/data/Control_Rep1_R2.fastq.gz
Control_Rep2,Control,/path/to/data/Control_Rep2_R1.fastq.gz,/path/to/data/Control_Rep2_R2.fastq.gz
Treated_Rep1,Treated,/path/to/data/Treated_Rep1_R1.fastq.gz,/path/to/data/Treated_Rep1_R2.fastq.gz
Treated_Rep2,Treated,/path/to/data/Treated_Rep2_R1.fastq.gz,/path/to/data/Treated_Rep2_R2.fastq.gz
```

### 2\. Prepare the Reference Genome

The pipeline requires several pre-built reference files. Below are the steps to generate them for the human hg38 genome using GENCODE annotations.

#### Create Reference Directory

Create a dedicated directory for all reference data:

```bash
mkdir -p reference_data
cd reference_data
```

#### Common Reference Files for Both Modes

Both `align` and `quant` modes require the following base files:

**Download Reference Files:**

```bash
# Download Genome FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz

# Download GTF Annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.annotation.gtf.gz

# Unzip the files
gunzip GRCh38.primary_assembly.genome.fa.gz
gunzip gencode.v46.primary_assembly.annotation.gtf.gz

# Rename annotation file for consistency
mv gencode.v46.primary_assembly.annotation.gtf annotation.gtf
```

#### Reference Files for Align Mode

The `align` mode requires STAR indices for both genome alignment and rRNA depletion.

##### Create rRNA Reference:

```bash
# Create a list of rRNA intervals from the GTF for STAR's rRNA index
grep 'gene_type "rRNA"' annotation.gtf > gencode.v46.rRNA.gtf
```

##### Build STAR Indices:

```bash
# Build the main Genome Index
mkdir -p star_genome_index
singularity exec RNA.sif STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ./star_genome_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile annotation.gtf \
     --sjdbOverhang 149

# Build the rRNA Index
mkdir -p star_rrna_index
singularity exec RNA.sif STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ./star_rrna_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v46.rRNA.gtf \
     --sjdbOverhang 149
```

The final reference structure for `align` mode:

```
reference_data/
├── annotation.gtf
├── GRCh38.primary_assembly.genome.fa
├── gencode.v46.rRNA.gtf
├── star_genome_index/
└── star_rrna_index/
```

#### Reference Files for Quant Mode

The `quant` mode requires STAR rRNA index, Salmon index and transcript-to-gene mapping.

##### Download Transcriptome FASTA:

```bash
# Download transcriptome FASTA file from GENCODE
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.transcripts.fa.gz
gunzip gencode.v46.transcripts.fa.gz
```

##### Build Salmon Index:

```bash
# Build the Salmon Index
mkdir -p salmon_index
singularity exec RNA.sif salmon index -t gencode.v46.transcripts.fa -i ./salmon_index -k 31
```

##### Build STAR rRNA Index:

```bash
mkdir -p star_rrna_index
singularity exec RNA.sif STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir ./star_rrna_index \
     --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile gencode.v46.rRNA.gtf \
     --sjdbOverhang 149
```

##### Create Transcript-to-Gene Map:

```bash
# Create the transcript-to-gene map (t2g.tsv)
# This command extracts transcript and gene IDs from the transcriptome FASTA header
grep "^>" gencode.v46.transcripts.fa \
  | cut -d'|' -f1,2 \
  | sed 's/>//g' \
  | sed 's/\.[0-9]*//g' \
  | tr '|' '\t' > t2g.tsv
```

The final reference structure for `quant` mode:

```
reference_data/
├── annotation.gtf
├── gencode.v46.transcripts.fa
├── t2g.tsv
├── salmon_index/
└── star_rrna_index/
```

## Running

Execute the pipeline using a single command, providing the sample sheet, output directory, reference directory, and desired mode.

### Command Parameters

  - `-s`: Path to the sample sheet CSV file (required)
  - `-o`: Output directory path where results will be saved (required)
  - `-r`: Reference data directory containing indices and annotation files (required)
  - `-m`: Analysis mode, either `align` or `quant` (required)
  - `-c`: Path to the RNA.sif Singularity container file (required)
  - `-L`: Control condition name for differential expression analysis (optional,If unset, DESeq2 uses alphabetical order)
  - `-t`: Number of threads to use for processing (optional, default: 8)
  - `-g`: Path to transcript-to-gene mapping file (required only for `quant` mode)

### Example Commands

**Example for `align` mode:**

```bash
./run_pipeline.sh \
  -s ./samplesheet.csv \
  -o ./project_results \
  -r ./reference_data \
  -m align \
  -c ./RNA.sif \
  -L Control \
  -t 8
```

**Example for `quant` mode:**

Note the additional `-g` flag, which is required for this mode to provide the transcript-to-gene mapping file.

```bash
./run_pipeline.sh \
  -s ./samplesheet.csv \
  -o ./project_results \
  -r ./reference_data \
  -m quant \
  -c ./RNA.sif \
  -g ./reference_data/t2g.tsv \
  -L Control \
  -t 8
```

## Output Structure and Interpretation

After the pipeline completes, the output directory will contain several files and directories. Below is a detailed explanation of what each file is and how it can be used.

-----

### Align Mode Output

```
./project_results/
├── Control_Rep1/
│   ├── Control_Rep1.dedup.bam         # Final processed BAM file
│   ├── Control_Rep1.dedup.bam.bai     # BAM index file
│   └── Control_Rep1.bw                # BigWig signal track
├── Control_Rep2/
│   ├── Control_Rep2.dedup.bam         # Final processed BAM file
│   ├── Control_Rep2.dedup.bam.bai     # BAM index file
│   └── Control_Rep2.bw                # BigWig signal track
├── Treated_Rep1/
│   ├── Treated_Rep1.dedup.bam         # Final processed BAM file
│   ├── Treated_Rep1.dedup.bam.bai     # BAM index file
│   └── Treated_Rep1.bw                # BigWig signal track
├── Treated_Rep2/
│   ├── Treated_Rep2.dedup.bam         # Final processed BAM file
│   ├── Treated_Rep2.dedup.bam.bai     # BAM index file
│   └── Treated_Rep2.bw                # BigWig signal track
├── multiqc_report/
│   └── multiqc_report.html            # Aggregated QC report for all samples
├── deg_results.txt                    # Differential expression gene list from DESeq2
└── normalized_counts.txt              # Normalized counts matrix from DESeq2
```

#### Per-Sample Files (`Control_Rep1/`)

  - **`*.dedup.bam`**

      - **Content**: This is the main alignment file in Binary Alignment Map (BAM) format. It contains all the sequencing reads and their mapping coordinates on the reference genome. This version has had duplicate reads (PCR duplicates) removed.
      - **Application**: It's the primary evidence for read alignment and can be used for detailed inspection in genome browsers or for downstream analyses.

  - **`*.dedup.bam.bai`**

      - **Content**: This is the index file for the BAM file.
      - **Application**: It allows for fast random access to the BAM file, which is essential for visualization software (like IGV) to quickly load and display alignments for a specific genomic region without reading the entire file.

  - **`*.bw`**

      - **Content**: A BigWig file that represents the RNA-seq signal coverage across the genome. It shows the read density (how many reads cover each position) in a compressed format.
      - **Application**: Primarily used for visualization. You can load this file into a genome browser (e.g., IGV, UCSC Genome Browser) to see a "signal track" that shows gene expression levels visually across chromosomes. Highly expressed genes will appear as peaks.


<img width="2386" height="534" alt="CleanShot 2025-09-13 at 15 31 03@2x" src="https://github.com/user-attachments/assets/e7fa1554-dfd0-47fb-b6ee-d16c50dba478" />


#### Aggregate Result Files

  - **`deg_results.txt`**

      - **Content**: A tab-separated text file containing the results of the differential expression analysis from DESeq2. Each row corresponds to a gene, and columns typically include:
          - `baseMean`: Average normalized count across all samples.
          - `log2FoldChange`: The logarithm (base 2) of the fold change between the 'Treated' and 'Control' conditions. A positive value means the gene is upregulated in the 'Treated' group; a negative value means it is downregulated.
          - `lfcSE`: The standard error of the `log2FoldChange` estimate.
          - `stat`: The Wald statistic.
          - `pvalue`: The raw p-value for the statistical test.
          - `padj`: The p-value adjusted for multiple testing (e.g., using Benjamini-Hochberg correction).
      - **Application**: This is the final result file. You can filter this file based on `padj` (e.g., `padj < 0.05`) and `log2FoldChange` thresholds to obtain a list of statistically significant differentially expressed genes.

  - **`normalized_counts.txt`**

      - **Content**: A tab-separated text file containing a matrix of normalized expression counts. Rows represent genes, and columns represent samples. These counts are adjusted for differences in sequencing depth between libraries, making them comparable across samples.
          - `row names`: Gene name.
          - `Control_Rep1`: The normalized count for Control_Rep1.
          - `Control_Rep2`: The normalized count for Control_Rep2.
          - `Treated_Rep1`: The normalized count for Treated_Rep1.
          - `Treated_Rep2`: The normalized count for Treated_Rep2.
      - **Application**: This matrix is essential for downstream analyses and visualizations beyond simple DEG lists. It can be used as input for generating heatmaps, performing principal component analysis (PCA) to check for sample clustering, or conducting gene set enrichment analysis.

  - **`multiqc_report`**: Open `multiqc_report.html` in a web browser to explore all sections interactively.

  	- **Application**: This is the first file you should check to assess the overall quality of your sequencing data and the alignment process. It helps identify problematic samples (e.g., low alignment rate, high duplication) early on.

    	- **General Statistics**: A combined table summarizing important metrics for each sample:

      <img width="2242" height="918" alt="CleanShot 2025-09-13 at 15 41 20@2x" src="https://github.com/user-attachments/assets/92cc2049-51aa-4ac6-b4f8-4c531637904a" />

    	- **FastQC**: Quality-control metrics on raw and trimmed reads, including  
      'Sequence Counts', 'Sequence Quality Histograms', 'Per Sequence Quality Scores',  
      'Per Base Sequence Content', 'Per Sequence GC Content', 'Per Base N Content',  
      'Sequence Length Distribution', 'Sequence Duplication Levels',  
      'Overrepresented sequences by sample', 'Top overrepresented sequences', 'Adapter Content'.

      	- **Sequence Quality Histograms**: The mean quality value across each base position in the read.  

        <img width="2222" height="616" alt="CleanShot 2025-09-13 at 16 12 48@2x" src="https://github.com/user-attachments/assets/45253b25-d3b1-4baa-8907-3f75fc564c25" />

      	- **Adapter Content**: The cumulative percentage count of the proportion of your library which has seen each of the adapter sequences at each position.  

        <img width="2434" height="700" alt="CleanShot 2025-09-13 at 16 21 40@2x" src="https://github.com/user-attachments/assets/eb0aeb6c-28b6-408c-a0b7-1c0eb251b7fd" />

    	- **Cutadapt**: Reports the number of reads and bases trimmed for adapters and quality:

      <img width="2238" height="1272" alt="CleanShot 2025-09-13 at 15 38 41@2x" src="https://github.com/user-attachments/assets/b5428467-3bd7-4f86-bd68-ab236f9793a6" />

    	- **STAR**: Alignment statistics such as total reads, uniquely mapped reads, and multi-mapping rates:

      <img width="1900" height="1650" alt="CleanShot 2025-09-13 at 15 37 58@2x" src="https://github.com/user-attachments/assets/4c71ed96-0137-4a29-82c9-9e0c6b7f52c6" />

    	- **featureCounts**: Gene-level quantification results, including total counts and assignment rates:

      <img width="1904" height="1026" alt="CleanShot 2025-09-13 at 15 37 09@2x" src="https://github.com/user-attachments/assets/a94c1a74-1d74-422c-95ba-2397a273b10c" />
	


-----

### Quant Mode Output

The output is more concise as it does not generate per-sample alignment files.

```
./project_results/
├── multiqc_report/
│   └── multiqc_report.html            # Aggregated QC report for all samples
├── deg_results.txt                    # Differential expression gene list from DESeq2
└── normalized_counts.txt              # Normalized counts matrix from DESeq2
```

  - **`deg_results.txt`**: Same as in `align` mode. It lists the differentially expressed genes based on transcript-level quantifications that have been summarized to the gene level.
  - **`normalized_counts.txt`**: Same as in `align` mode. This matrix of normalized counts is derived from Salmon's quantifications and is ready for downstream visualization and analysis like PCA and heatmap generation.
  - **`multiqc_report`** : Open multiqc_report.html in a web browser to explore all sections interactively. (Same as in `align` mode, but without STAR alignment statistics, instead, Salmon statistics are available.)
 
       - **Salmon**: Fragment length distribution – Shows the estimated insert-size (fragment length) distribution of the sequencing library:
	  This reflects the typical distance between read pairs after library preparation and is important for accurate transcript abundance estimation. Abnormal distributions may indicate library prep issues or unexpected sample characteristics.
		<img width="1772" height="1292" alt="CleanShot 2025-09-13 at 15 36 04@2x" src="https://github.com/user-attachments/assets/087233e4-380a-4a1d-a12c-de70db7c06c8" />


## Video Tutorials

### Align Mode Tutorial

Watch this video tutorial to see a complete walkthrough of running the pipeline in `align` mode:



https://github.com/user-attachments/assets/2009351c-23d9-4d23-b56d-2b6f6538d318



### Quant Mode Tutorial

Watch this video tutorial to see a complete walkthrough of running the pipeline in `quant` mode:



https://github.com/user-attachments/assets/774937f6-2990-4176-a879-141dd138d99d



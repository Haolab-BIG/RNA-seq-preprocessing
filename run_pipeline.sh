#!/bin/bash

# --- Stricter Shell Safety ---
# Exit immediately if a command exits with a non-zero status.
set -e
# Fail a pipeline if any command fails, not just the last one.
set -o pipefail
# Set a sane default for the Internal Field Separator.
IFS=$'\n\t'

# Author: ZHANG HE
#
# Description:
#   This script executes a complete RNA-seq workflow with a single command. It
#   incorporates several improvements for robustness and flexibility, including
#   shell practices, and more robust file handling.

# --- 1. Default Parameters & Help Function ---
THREADS=8
MODE=""
REF_LEVEL=""
# MODIFICATION: Set a default for the SIF path, can be overridden by -c
SIF_PATH="$(dirname "$0")/RNA.sif"


HELP_MSG="Usage: $0 -s <sample_sheet.csv> -o <out_dir> -r <ref_dir> -m <'align'|'quant'> [OPTIONS]

Required:
  -s  Sample sheet (CSV). Format: sample,condition,fastq1_path,fastq2_path
  -o  Main output directory.
  -r  Reference data directory.
  -m  Workflow mode: 'align' (alignment-based) or 'quant' (alignment-free).

Optional:
  -g  Transcript-to-gene map file (Required when -m is 'quant').
  -c  Path to the Singularity container (RNA.sif). (Default: same directory as script)
  -t  Number of threads to use (Default: 8).
  -L  Reference level for DESeq2 comparison (e.g., 'Control'). If unset, DESeq2 uses alphabetical order.
  -h  Display this help message.
"

# --- 2. Parse Command-line Arguments ---
T2G_MAP=""
# MODIFICATION: Added 'c:' to getopts string
while getopts "s:o:r:m:g:c:t:L:h" opt; do
  case ${opt} in
    s ) SAMPLE_SHEET=$(realpath "${OPTARG}") ;;
    o ) OUT_DIR=$(realpath "${OPTARG}") ;;
    r ) REF_DIR=$(realpath "${OPTARG}") ;;
    m ) MODE=${OPTARG} ;;
    g ) T2G_MAP=$(realpath "${OPTARG}") ;;
    # MODIFICATION: Handle the new -c argument for the SIF path
    c ) SIF_PATH=$(realpath "${OPTARG}") ;;
    t ) THREADS=${OPTARG} ;;
    L ) REF_LEVEL=${OPTARG} ;;
    h ) echo "${HELP_MSG}"; exit 0 ;;
    \? ) echo "Invalid option: -${OPTARG}" >&2; echo "${HELP_MSG}"; exit 1 ;;
  esac
done

# Check for mandatory arguments
if [ -z "${SAMPLE_SHEET}" ] || [ -z "${OUT_DIR}" ] || [ -z "${REF_DIR}" ] || [ -z "${MODE}" ]; then
    echo "Error: Missing mandatory arguments." >&2; echo "${HELP_MSG}"; exit 1
fi
if [[ "${MODE}" != "align" && "${MODE}" != "quant" ]]; then
    echo "Error: Mode (-m) must be either 'align' or 'quant'." >&2; exit 1
fi
if [[ "${MODE}" == "quant" && -z "${T2G_MAP}" ]]; then
    echo "Error: In 'quant' mode, a transcript-to-gene map must be provided with -g." >&2; exit 1
fi
if [ ! -f "${SAMPLE_SHEET}" ]; then
    echo "Error: Sample sheet not found: ${SAMPLE_SHEET}" >&2; exit 1
fi

# --- 3. Setup Environment and Directories ---
# MODIFICATION: The SIF_PATH is now validated here, whether it's the default or user-provided.
if [ ! -f "${SIF_PATH}" ]; then
    echo "Error: Singularity container not found at: ${SIF_PATH}" >&2
    exit 1
fi
mkdir -p "${OUT_DIR}"

# Define the main Singularity execution command template
SINGULARITY_BASE_CMD="singularity exec --cleanenv -B ${OUT_DIR}:/output -B ${REF_DIR}:/reference"

echo "================================================="
echo "====== RNA-seq Automated Pipeline Started ======="
echo "================================================="
echo "Mode: ${MODE}"
echo "Sample Sheet: ${SAMPLE_SHEET}"
echo "Output Directory: ${OUT_DIR}"
echo "Threads: ${THREADS}"
echo "Singularity Container: ${SIF_PATH}"
if [ -n "${REF_LEVEL}" ]; then echo "DESeq2 Ref Level: ${REF_LEVEL}"; fi
echo "================================================="

# --- 4. Loop Through and Process Each Sample ---
# Use process substitution and 'tr' to robustly handle CSV files with CRLF line endings
while IFS=, read -r SAMPLE_NAME CONDITION FQ1 FQ2; do
    # Trim whitespace from all fields read from the sample sheet
    SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '[:space:]')
    CONDITION=$(echo "${CONDITION}" | tr -d '[:space:]')
    FQ1=$(echo "${FQ1}" | tr -d '[:space:]')
    FQ2=$(echo "${FQ2}" | tr -d '[:space:]')

    echo -e "\n\n======== Starting to process sample: [${SAMPLE_NAME}] ========"
    
    if [ ! -f "${FQ1}" ] || [ ! -f "${FQ2}" ]; then
        echo "Error: FASTQ files for sample ${SAMPLE_NAME} not found. Skipping." >&2
        continue
    fi

    # --- 4.1. Set up directories and dynamic binds ---
    FQ_DIR=$(dirname "${FQ1}")
    SAMPLE_OUT_DIR="/output/${SAMPLE_NAME}" # Path inside the container
    mkdir -p "${OUT_DIR}/${SAMPLE_NAME}"   # Create directory on the host

    SINGULARITY_CMD="singularity exec --cleanenv -B "${OUT_DIR}:/output" -B "${REF_DIR}:/reference" -B "${FQ_DIR}:/reads" "${SIF_PATH}""

    # --- 4.2. Shared Pipeline Steps (QC, Trim, rRNA depletion) ---
    echo "======== [${SAMPLE_NAME}] Step 1 & 2: Trimming and Raw QC ========"
    # Use --cores, the correct flag for parallel execution in Trim Galore
    eval "$SINGULARITY_CMD fastqc -t "${THREADS}" "${FQ1}" "${FQ2}" -o ${SAMPLE_OUT_DIR}"
    eval "$SINGULARITY_CMD trim_galore --paired --fastqc --cores "${THREADS}" \
        -o ${SAMPLE_OUT_DIR} \
        /reads/$(basename "${FQ1}") /reads/$(basename "${FQ2}")"

    # Robustly determine trimmed file names
    FQ1_BASENAME=$(basename "${FQ1}")
    FQ2_BASENAME=$(basename "${FQ2}")
    FQ1_TRIMMED_BASENAME=$(echo "${FQ1_BASENAME}" | sed -E 's/(\.fq|\.fastq)(\.gz)?$/_val_1.fq.gz/')
    FQ2_TRIMMED_BASENAME=$(echo "${FQ2_BASENAME}" | sed -E 's/(\.fq|\.fastq)(\.gz)?$/_val_2.fq.gz/')
    TRIMMED_FQ1="${SAMPLE_OUT_DIR}/${FQ1_TRIMMED_BASENAME}"
    TRIMMED_FQ2="${SAMPLE_OUT_DIR}/${FQ2_TRIMMED_BASENAME}"

    echo "======== [${SAMPLE_NAME}] Step 3: rRNA Depletion using STAR ========"
    # Removed non-standard --outSAMattributes flag
    eval "$SINGULARITY_CMD STAR \
        --runThreadN "${THREADS}" --genomeDir /reference/star_rrna_index \
        --readFilesIn "${TRIMMED_FQ1}" "${TRIMMED_FQ2}" --readFilesCommand zcat \
        --outStd BAM_Unsorted --outSAMtype BAM Unsorted \
        --outFilterMultimapNmax 20 --outFilterScoreMin 10 --outReadsUnmapped Fastx \
        --outFileNamePrefix ${SAMPLE_OUT_DIR}/rrna_depletion. \
        > /dev/null 2>&1"

    NON_RRNA_FQ1="${SAMPLE_OUT_DIR}/rrna_depletion.Unmapped.out.mate1"
    NON_RRNA_FQ2="${SAMPLE_OUT_DIR}/rrna_depletion.Unmapped.out.mate2"

    echo "======== [${SAMPLE_NAME}] Step 4: Post-depletion FastQ QC ========"
    eval "$SINGULARITY_CMD fastqc -t "${THREADS}" "${NON_RRNA_FQ1}" "${NON_RRNA_FQ2}" -o ${SAMPLE_OUT_DIR}"

    # Compress unmapped reads to save space
    echo "--> Compressing unmapped reads..."
    eval "$SINGULARITY_CMD gzip "${NON_RRNA_FQ1}""
    eval "$SINGULARITY_CMD gzip "${NON_RRNA_FQ2}""
    NON_RRNA_FQ1="${NON_RRNA_FQ1}.gz"
    NON_RRNA_FQ2="${NON_RRNA_FQ2}.gz"

    # --- 4.3. Mode-Specific Branched Steps ---
    if [ "${MODE}" == "align" ]; then
        echo "======== [${SAMPLE_NAME}] EXECUTING ALIGNMENT-BASED PATH ========"
        STAR_OUT_PREFIX="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}."
        
        echo "--> Step 5a: Genome Alignment (STAR)"
        eval "$SINGULARITY_CMD STAR \
            --runThreadN "${THREADS}" --genomeDir /reference/star_genome_index \
            --readFilesIn "${NON_RRNA_FQ1}" "${NON_RRNA_FQ2}" --readFilesCommand zcat \
            --outFileNamePrefix "${STAR_OUT_PREFIX}" --outSAMtype BAM SortedByCoordinate"
        
        SORTED_BAM="${STAR_OUT_PREFIX}Aligned.sortedByCoord.out.bam"
        UNIQUE_BAM="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.unique.bam"
        echo "--> Step 5b: Filter for Uniquely Mapped Reads and Sort"
        UNIQUE_SORTED_BAM="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.unique.sorted.bam"
        
        eval "$SINGULARITY_CMD samtools view -q 255 -b "${SORTED_BAM}"" | \
            eval "$SINGULARITY_CMD samtools sort -n -@ ${THREADS} -o "${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.namesort.bam" -"

        eval "$SINGULARITY_CMD samtools fixmate -m "${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.namesort.bam" "${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.fixmate.bam""
        eval "$SINGULARITY_CMD samtools sort -@ ${THREADS} -o ${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.fixmate.sorted.bam ${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.fixmate.bam"
        echo "--> Step 5c: Remove PCR Duplicates"
        DEDUP_BAM="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.dedup.bam"
        eval "$SINGULARITY_CMD samtools markdup -@ ${THREADS} -r "${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.fixmate.sorted.bam" "${DEDUP_BAM}""
        eval "$SINGULARITY_CMD samtools index "${DEDUP_BAM}""

        BIGWIG_FILE="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.bw"
        echo "--> Step 5d: Generate BigWig"
        eval "$SINGULARITY_CMD bamCoverage -b "${DEDUP_BAM}" -o "${BIGWIG_FILE}" -p "${THREADS}" --normalizeUsing CPM"
        
        COUNTS_FILE="${SAMPLE_OUT_DIR}/${SAMPLE_NAME}.counts.txt"
        echo "--> Step 5e: Gene Counting (featureCounts)"
        eval "$SINGULARITY_CMD featureCounts \
            -T "${THREADS}" -p -s 0 \
            -a /reference/annotation.gtf \
            -o "${COUNTS_FILE}" \
            "${DEDUP_BAM}""

    elif [ "${MODE}" == "quant" ]; then
        echo "======== [${SAMPLE_NAME}] EXECUTING ALIGNMENT-FREE PATH ========"
        echo "--> Step 5q: Salmon Quantification"
        # Salmon natively handles gzipped input
        eval "$SINGULARITY_CMD salmon quant \
            -i /reference/salmon_index -l A -1 "${NON_RRNA_FQ1}" -2 "${NON_RRNA_FQ2}" \
            --validateMappings --gcBias --threads "${THREADS}" \
            -o "${SAMPLE_OUT_DIR}""
    fi
    echo "======== [${SAMPLE_NAME}] Processing finished. ========"
done < <(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r')

echo -e "\n\n==============================================="
echo "======= All samples processed. Starting final analysis... ======="
echo "==============================================="

# --- 5. Embed and Execute DESeq2 R Script ---
echo "======== Step 5: Differential Expression Analysis (DESeq2) ========"
DESEQ_R_SCRIPT="${OUT_DIR}/deseq_analysis_runner.R"
DESEQ_SAMPLESHEET="${OUT_DIR}/deseq_samplesheet.csv"
OUTPUT_PREFIX="/output/final_results"

echo "sample,condition" > "${DESEQ_SAMPLESHEET}"
# Re-parse the sample sheet to create the simple version for DESeq2
while IFS=, read -r SAMPLE_NAME CONDITION _; do
    SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '[:space:]')
    CONDITION=$(echo "${CONDITION}" | tr -d '[:space:]')
    echo "${SAMPLE_NAME},${CONDITION}" >> "${DESEQ_SAMPLESHEET}"
done < <(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r')

cat > "${DESEQ_R_SCRIPT}" << 'EOF'
#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tximport))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript deseq_runner.R <mode> <sample_sheet> <data_path> <output_prefix> [t2g_path] [ref_level]", call. = FALSE)
}

MODE <- args[1]
SAMPLE_SHEET_PATH <- args[2]
DATA_PATH <- args[3]
OUTPUT_PREFIX <- args[4]
T2G_PATH <- if (length(args) >= 5) args[5] else ""
REF_LEVEL <- if (length(args) >= 6) args[6] else ""

message("Reading sample sheet from: ", SAMPLE_SHEET_PATH)
sample_table <- read.csv(SAMPLE_SHEET_PATH, header = TRUE, row.names = 1)

if (MODE == "featurecounts") {
    message("Mode: featurecounts. Finding and merging individual count files from: ", DATA_PATH)
    count_files <- file.path(DATA_PATH, rownames(sample_table), paste0(rownames(sample_table), ".counts.txt"))
    names(count_files) <- rownames(sample_table)
    if (!all(file.exists(count_files))) {
        stop("Error: Missing count files: ", paste(names(count_files)[!file.exists(count_files)], collapse=", "))
    }
    list_of_dfs <- lapply(names(count_files), function(sample_name) {
        df <- read.table(count_files[[sample_name]], header = TRUE, sep = "\t", comment.char = "#")
        df <- df[, c(1, 7)]
        colnames(df) <- c("Geneid", sample_name)
        return(df)
    })
    count_df <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), list_of_dfs)
    rownames(count_df) <- count_df$Geneid
    count_df$Geneid <- NULL
    count_df[is.na(count_df)] <- 0
    count_data <- as.matrix(count_df)

} else if (MODE == "salmon") {
    message("Mode: salmon. Importing quant files from: ", DATA_PATH)
    if (T2G_PATH == "") stop("Transcript-to-gene map (t2g_path) is required for Salmon mode.")
    files <- file.path(DATA_PATH, rownames(sample_table), "quant.sf")
    names(files) <- rownames(sample_table)
    
    fix_ids <- function(f) {
        q <- read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        q$Name <- sub("\\|.*", "", q$Name)   # 去掉 | 后的所有部分
        q$Name <- sub("\\..*", "", q$Name)   # 去掉版本号 (.2)
        tmpfile <- tempfile(fileext = ".sf")
        write.table(q, tmpfile, quote = FALSE, row.names = FALSE, sep = "\t")
        return(tmpfile)
    }
    fixed_files <- sapply(files, fix_ids, USE.NAMES = TRUE)

    t2g <- read.table(T2G_PATH, header = FALSE, col.names = c("TXNAME", "GENEID"))
    txi <- tximport(fixed_files, type = "salmon", tx2gene = t2g,
                    ignoreTxVersion = TRUE, ignoreAfterBar = TRUE)
    count_data <- txi$counts
}

count_data <- count_data[, rownames(sample_table)]
dds <- DESeqDataSetFromMatrix(countData = round(count_data), colData = sample_table, design = ~ condition)

if (REF_LEVEL != "" && REF_LEVEL %in% dds$condition) {
    message("Setting reference level for condition to: ", REF_LEVEL)
    dds$condition <- relevel(dds$condition, ref = REF_LEVEL)
} else if (REF_LEVEL != "") {
    warning("Provided reference level '", REF_LEVEL, "' not found in conditions. Using alphabetical default.")
}

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
message("Running DESeq2 analysis...")
dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj, na.last = TRUE),]

message("Writing results to output files prefixed with: ", OUTPUT_PREFIX)
write.table(as.data.frame(counts(dds, normalized=TRUE)),
            file = paste0(OUTPUT_PREFIX, "_normalized_counts.txt"),
            sep = "\t", quote = FALSE)
write.table(res_df,
            file = paste0(OUTPUT_PREFIX, "_DEG_results.txt"),
            sep = "\t", quote = FALSE)
message("Analysis complete.")
EOF

# Build and execute the Singularity command for DESeq2
if [ "${MODE}" == "align" ]; then
    DESEQ2_MODE="featurecounts"
    eval "$SINGULARITY_BASE_CMD "${SIF_PATH}" Rscript /output/deseq_analysis_runner.R \
        "${DESEQ2_MODE}" \
        /output/deseq_samplesheet.csv \
        /output/ \
        "${OUTPUT_PREFIX}" \
        "unused" \
        "${REF_LEVEL}""

elif [ "${MODE}" == "quant" ]; then
    DESEQ2_MODE="salmon"
    T2G_DIR=$(dirname "${T2G_MAP}")
    T2G_BASENAME=$(basename "${T2G_MAP}")
    # Add a specific bind for the t2g map's directory
    eval "$SINGULARITY_BASE_CMD -B "${T2G_DIR}":/t2g_dir "${SIF_PATH}" Rscript /output/deseq_analysis_runner.R \
        "${DESEQ2_MODE}" \
        /output/deseq_samplesheet.csv \
        /output/ \
        "${OUTPUT_PREFIX}" \
        "/t2g_dir/${T2G_BASENAME}" \
        "${REF_LEVEL}""
fi

mv "${OUT_DIR}/final_results_DEG_results.txt" "${OUT_DIR}/deg_results.txt"
mv "${OUT_DIR}/final_results_normalized_counts.txt" "${OUT_DIR}/normalized_counts.txt"
echo "Differential expression results saved to: ${OUT_DIR}/deg_results.txt"

# --- 6. Generate the Final MultiQC Report ---
echo "======== Step 6: Generating Final MultiQC Report ========"
mkdir -p "${OUT_DIR}/multiqc_report"
eval "$SINGULARITY_BASE_CMD "${SIF_PATH}" multiqc /output -o /output/multiqc_report --force"


# --- 7. MODIFICATION: Final Cleanup of Sample Directories ---
echo "======== Step 7: Cleaning up intermediate files ========"
while IFS=, read -r SAMPLE_NAME CONDITION FQ1 FQ2; do
    SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '[:space:]')
    SAMPLE_DIR="${OUT_DIR}/${SAMPLE_NAME}"

    if [ ! -d "${SAMPLE_DIR}" ]; then
        continue
    fi

    echo "--> Cleaning directory: ${SAMPLE_DIR}"
    if [ "${MODE}" == "align" ]; then
        # In align mode, keep only the final BAM, its index, and the BigWig file.
        # The `find` command deletes all files that do NOT match the kept patterns.
        find "${SAMPLE_DIR}" -type f \
            ! -name "*.dedup.bam" \
            ! -name "*.dedup.bam.bai" \
            ! -name "*.bw" \
            -delete
        rm -f "${OUT_DIR}/deseq_analysis_runner.R"
        rm -f "${OUT_DIR}/deseq_samplesheet.csv"
    elif [ "${MODE}" == "quant" ]; then
        # In quant mode, no bam/bw files are generated that need to be kept.
        # The entire sample-specific directory is removed.
        rm -rf "${SAMPLE_DIR}"
        rm -f "${OUT_DIR}/deseq_analysis_runner.R"
        rm -f "${OUT_DIR}/deseq_samplesheet.csv"
    fi
done < <(tail -n +2 "${SAMPLE_SHEET}" | tr -d '\r')
echo "======== Cleanup complete ========"


echo -e "\n\n==============================================="
echo "====== RNA-seq Pipeline Completed Successfully! ======"
echo "==============================================="
echo "All results are located in: ${OUT_DIR}"
echo "View the interactive QC report at: ${OUT_DIR}/multiqc_report/multiqc_report.html"
echo "==============================================="
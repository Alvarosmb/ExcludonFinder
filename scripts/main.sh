#!/bin/bash
# Function to check if command exists
command_exists() {
  command -v "$1" >/dev/null 2>&1
}

# Check for required commands
for cmd in conda bwa-mem2 samtools minimap2; do
  if ! command_exists "$cmd"; then
    echo "Error: Required command '$cmd' not found" >&2
    exit 1
  fi
done



# Set default values for options
threshold=0.5
N_threads=8
alignment_tool="bwa-mem2"
Paired="FALSE"

# Variables to store input files
fasta_input=""
fastq_input1=""
fastq_input2=""
gff_input=""
use_minimap2=false
run_control_script=false

while getopts "f:1:2:g:t:j:lC" opt; do  # Changed -r to -t and -t to -j
  case ${opt} in
    f ) fasta_input=$(realpath "${OPTARG}");;
    1 ) fastq_input1=$(realpath "${OPTARG}");;
    2 ) fastq_input2=$(realpath "${OPTARG}");;
    g ) gff_input=$(realpath "${OPTARG}");;
    t ) threshold=$OPTARG;;   
    j ) N_threads=$OPTARG;;   
    l ) use_minimap2=true;;
    C ) run_control_script=true;;
    \? ) echo "Invalid option: -$OPTARG" 1>&2; exit 1;;
    : ) echo "Option -$OPTARG requires an argument." 1>&2; exit 1;;
  esac
done
# Check that mandatory input files are specified
if [ -z "$fasta_input" ] || [ -z "$gff_input" ]; then
  echo "Usage: $0 -f <fasta_input> -1 <fastq_input_1> [-2 <fastq_input_2>] -g <gff_input> [-t coverage threshol] [-l]" >&2
  exit 1
fi

# Check if input files exist and are readable
for file in "$fasta_input" "$gff_input" "$fastq_input1"; do
  if [ ! -r "$file" ]; then
    echo "Error: Cannot read file: $file" >&2
    exit 1
  fi
done

if [ -n "$fastq_input2" ] && [ ! -r "$fastq_input2" ]; then
  echo "Error: Cannot read file: $fastq_input2" >&2
  exit 1
fi

# Handle gzipped FASTQ files
process_fastq() {
  local input=$1
  if [[ $input == *.gz ]]; then
    gunzip -c "$input" > "${input%.gz}"
    echo "${input%.gz}"
  else
    echo "$input"
  fi
}

# Process input files
fastq_input1_processed=$(process_fastq "$fastq_input1")
if [ -n "$fastq_input2" ]; then
  fastq_input2_processed=$(process_fastq "$fastq_input2")
  Paired="TRUE"
  echo "PAIRED_END" >&2
fi

# Replace ID for gene_id if necessary
sed -i 's/ID/gene_id/g' "$gff_input"

# Set up alignment based on paired/single end
if [ -z "$fastq_input2" ]; then
  filename=$(basename "$fastq_input1")
  sample=${filename%.fastq*}
  input_fastq_option="$fastq_input1_processed"
else
  filename1=$(basename "$fastq_input1")
  sample=${filename1%_R1.fastq*}
  input_fastq_option="$fastq_input1_processed $fastq_input2_processed"
fi

# Create output directories
mkdir -p output/{alignment,coverage_data,index}

# Perform alignment
echo "STARTING_ALIGNMENT" >&2
if $use_minimap2; then
  echo "Data type: Long reads detected" >&2
  # Modified minimap2 command with additional error checking
  if ! minimap2 -t "$N_threads" -a -p 0.99 -k14 --MD -uf "$fasta_input" "$fastq_input1_processed" > "output/alignment/${sample}_temp.sam" 2>/dev/null; then
    echo "Error: minimap2 alignment failed" >&2
    exit 1
  fi
  
  # Convert SAM to sorted BAM
  if ! samtools view -@ "$N_threads" -bh "output/alignment/${sample}_temp.sam" 2>/dev/null | \
       samtools sort -@ "$N_threads" -o "output/alignment/${sample}_sorted.bam" 2>/dev/null; then
    echo "Error: SAM to BAM conversion failed" >&2
    exit 1
  fi
  
  # Clean up temporary SAM file
  rm "output/alignment/${sample}_temp.sam"
else
  # Create index files in the index directory
  index_prefix="output/index/$(basename "$fasta_input")"
  bwa-mem2 index -p "$index_prefix" "$fasta_input" >/dev/null 2>&1
  bwa-mem2 mem -t "$N_threads" "$index_prefix" $input_fastq_option 2>/dev/null | 
    samtools view -@ "$N_threads" -bh - 2>/dev/null | 
    samtools sort -@ "$N_threads" -o "output/alignment/${sample}_sorted.bam" 2>/dev/null
fi

# Index the final BAM file
if ! samtools index "output/alignment/${sample}_sorted.bam"; then
  echo "Error: BAM indexing failed" >&2
  exit 1
fi

# Get alignment stats
mapped_percent=$(samtools flagstat "output/alignment/${sample}_sorted.bam" | grep "mapped" | head -n1 | awk -F'[(%]' '{print $2}')
paired_percent=$(samtools flagstat "output/alignment/${sample}_sorted.bam" | grep "properly paired" | awk -F'[(%]' '{print $2}')

echo "MAPPED_READS:$mapped_percent" >&2
echo "PAIRED_READS:$paired_percent" >&2

bam="output/alignment/${sample}_sorted.bam"

if $run_control_script; then
    # Create a temporary file for both stdout and stderr
    temp_file=$(mktemp)
    
    # Run featureCounts and redirect both stdout and stderr to temp file
    if $use_minimap2; then
        featureCounts -a "$gff_input" -F GFF3 -t CDS -T "$N_threads" -L "$bam" -o "output/coverage_data/${sample}_counts.txt" > "$temp_file" 2>&1
    elif [ "$Paired" = "TRUE" ]; then
        featureCounts -a "$gff_input" -F GFF3 -t CDS -T "$N_threads" -p -B -C "$bam" -o "output/coverage_data/${sample}_counts.txt" > "$temp_file" 2>&1
    else
        featureCounts -a "$gff_input" -F GFF3 -t CDS -T "$N_threads" "$bam" -o "output/coverage_data/${sample}_counts.txt" > "$temp_file" 2>&1
    fi
    
    # Extract the alignment rate more robustly
    if ! alignment_rate=$(grep "Assignment rate" "$temp_file" | grep -o '[0-9]\+\.[0-9]\+'); then
        # If the first grep fails, try an alternative pattern
        alignment_rate=$(grep "Successfully assigned" "$temp_file" | grep -o '[0-9]\+\.[0-9]\+')
    fi
    
    # Check if we got a valid alignment rate
    if [ -n "$alignment_rate" ]; then
        if (( $(echo "$alignment_rate < 80" | bc -l) )); then
            echo "WARNING: Low percentage of mapped reads to CDSs detected (${alignment_rate}%)"
        else
            echo "PASS: Good rate of mapped reads to CDSs  (${alignment_rate}%)"
        fi
    else
        echo "WARNING: Could not determine alignment rate"
    fi
    
    # Clean up temp file
    rm -f "$temp_file"
fi

echo "CALCULATING_DEPTH" >&2
# Process reads based on strand
if [ "$Paired" = "TRUE" ]; then
    # Paired-end processing
    samtools view -h -f 0x40 -F 0x10 "$bam" > "output/alignment/${sample}_plus_strand_1.sam"
    samtools view -h -f 0x80 -f 0x10 "$bam" > "output/alignment/${sample}_plus_strand_2.sam"
    samtools merge -f "output/alignment/${sample}_plus_strand_merged.bam" \
        "output/alignment/${sample}_plus_strand_1.sam" \
        "output/alignment/${sample}_plus_strand_2.sam"

    samtools view -h -f 0x40 -f 0x10 "$bam" > "output/alignment/${sample}_minus_strand_1.sam"
    samtools view -h -f 0x80 -F 0x10 "$bam" > "output/alignment/${sample}_minus_strand_2.sam"
    samtools merge -f "output/alignment/${sample}_minus_strand_merged.bam" \
        "output/alignment/${sample}_minus_strand_1.sam" \
        "output/alignment/${sample}_minus_strand_2.sam"
else
    # Single-end processing
    samtools view -bh -F 0x14 "$bam" > "output/alignment/${sample}_plus_strand_merged.bam"
    samtools view -bh -f 0x10 -F 0x4 "$bam" > "output/alignment/${sample}_minus_strand_merged.bam"
fi

# Sort and index strands
for strand in plus minus; do
    samtools sort -@ "$N_threads" "output/alignment/${sample}_${strand}_strand_merged.bam" \
        -o "output/alignment/${sample}_${strand}_strand_sorted.bam" 2>/dev/null 
    samtools index "output/alignment/${sample}_${strand}_strand_sorted.bam"
    samtools depth -a "output/alignment/${sample}_${strand}_strand_sorted.bam" \
        > "output/coverage_data/${sample}_${strand}_depth.txt" 2>/dev/null
done

echo "DEPTH COMPLETED" >&2

echo "STARTING_ANNOTATION" >&2
conda deactivate

# Run Excludon annotation
Rscript scripts/Tus_annotation.R "$fasta_input" "$gff_input" "$bam" "$sample" "$threshold" "$N_threads"

# Cleanup temporary files
rm  -r "output/alignment/"
rm  -r "output/coverage_data"


echo "### The analysis is finished ###"

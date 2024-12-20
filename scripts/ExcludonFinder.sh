#!/bin/bash

print_header(){
  cat << 'EOF'
   ______          _           _               ______ _           _                        
  |  ____|        | |         | |             |  ____(_)         | |           
  | |__  __  _____| |_   _  __| | ___  _ __   | |__   _ _ __   __| | ___ _ __   
  |  __| \ \/ / __| | | | |/ _` |/ _ \| '_ \  |  __| | | '_ \ / _` |/ _ \ '__|   
  | |____ >  < (__| | |_| | (_| | (_) | | | | | |    | | | | | (_| |  __/ |     
  |______/_/\_\___|_|\__,_|\__,_|\___/|_| |_| |_|    |_|_| |_|\__,_|\___|_|      
     --------------------------------------------
               --------------------------------------------
                                           --------------------------------------------
                                                      --------------------------------------------
                                        V1.0.0                                                                                                                                                          
EOF
}

# Function to repeat a character n times
repeat_char() {
    local char="$1"
    local count="$2"
    printf "%${count}s" | tr " " "$char"
}

# Function to print formatted section headers
print_section() {
    local title="$1"
    local width=78
    local title_len=${#title}
    local padding=$(( (width - title_len) / 2 - 2 ))
    
    echo "//==============================================================================\\\\"
    echo "||                                                                              ||"
    printf "||%*s%s%*s||\n" $padding "" "$title" $(( width - title_len - padding  )) ""
    echo "||                                                                              ||"
}

# Function to print formatted status
print_status() {
    local message="$1"
    printf "||  %-76s||\n" "${message:0:76}"
}

# Function to print section footer
print_footer() {
    echo "||                                                                              ||"
    echo "\\\==============================================================================//"
    echo
}

# Function to print empty line
print_empty_line() {
    echo "||                                                                              ||"
}

# Function to parse command line arguments with improved error handling
parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            -f) shift; fasta_file="$1" ;;
            -g) shift; gff_file="$1" ;;
            -1) shift; fastq_r1="$1" ;;
            -2) shift; fastq_r2="$1" ;;
            -t) shift; threshold="$1" ;;
            -j) shift; threads="$1" ;;
            -l) long_reads="yes" ;;
            -C) qc="yes" ;;
            *) echo "Error: Unknown option: $1"; exit 1 ;;
        esac
        shift
    done
}

# Function to check if a file exists with proper error formatting
check_file() {
    local file="$1"
    local desc="$2"
    if [[ ! -f "$file" ]]; then
        print_status "Error: $desc file not found: $file"
        return 1
    fi
    return 0
}

# Function to check all requirements with consistent error reporting
check_requirements() {
    local errors=0

    # Check Scripts directory
    if [[ ! -d "scripts" ]]; then
        print_status "Error: scripts directory not found"
        return 1
    fi

    # Check mandatory input files
    check_file "$fasta_file" "FASTA" || ((errors++))
    check_file "$gff_file" "GFF" || ((errors++))
    check_file "$fastq_r1" "FASTQ R1" || ((errors++))
    
    # Only check R2 if it was provided
    if [[ -n "$fastq_r2" ]]; then
        check_file "$fastq_r2" "FASTQ R2" || ((errors++))
    fi

    # Create output directory structure
    mkdir -p output/{alignment,coverage_data} 2>/dev/null

    return $errors
}

move_cursor_up() {
    echo -en "\033[${1}A"
}

# Function to clear lines below current position
clear_lines_below() {
    local n=$1
    for ((i=0; i<n; i++)); do
        echo -en "\033[1B\033[2K"
    done
    move_cursor_up "$n"
}

main() {
    # Clear screen and print header
    clear
    print_header
    
    # Parse command line arguments
    parse_args "$@"
    
    # Check requirements
    if ! check_requirements; then
        print_section "Error"
        print_status "Fatal: Missing required files or directories"
        print_footer
        exit 1
    fi
    
    # Print configuration with consistent formatting
    print_section "Configuration"
    print_status "FASTA file  : $fasta_file"
    print_status "GFF file    : $gff_file"
    print_status "FASTQ R1    : $fastq_r1"
    [[ -n "$fastq_r2" ]] && print_status "FASTQ R2    : $fastq_r2"
    print_status "Threads     : ${threads:-8}"
    print_status "Threshold   : ${threshold:-0.5}"
    print_status "Data type   : $([[ -n "$fastq_r2" ]] && echo "Paired-end" || echo "Single-end")"
    [[ "$long_reads" == "yes" ]] && print_status "Long reads  : Yes"
    [[ "$qc" == "yes" ]] && print_status "QC enabled  : Yes"
    print_footer
    
    # Start progress monitoring with improved formatting
    print_section "Running pipeline"
    
    # Run the main script with formatted output processing
    bash scripts/main.sh "$@" 2>&1 | while IFS= read -r line; do
        case "$line" in
            "PAIRED_END")
                print_status "Data type: Paired-end reads detected"
                print_empty_line
                ;;
            "STARTING_ALIGNMENT")
                print_status "Performing Alignment..."
                ;;
            "CALCULATING_DEPTH")
                print_status "Calculating Coverage Depth For Each Nucleotide..."
                print_empty_line
                cat << 'EOF'
||                  ATGGCTAGCTTACGT                                             ||
||                  ..............                                              ||
||                  ........ .....                                              ||
||                  . .... . .. ..                                              ||
||                  . .... . .. ..                                              ||
EOF
                print_empty_line
                print_empty_line
                ;;
            "DEPTH COMPLETED")
                print_status "Completed"
                print_section "Detecting Overlapping Transcription"
                print_empty_line
                ;;
            MAPPED_READS:*)
                percentage=${line#MAPPED_READS:}
                print_status "    Mapped reads: ${percentage}%"
                ;;
            PAIRED_READS:*)
                if [[ -n "$fastq_r2" ]]; then
                    percentage=${line#PAIRED_READS:}
                    print_status "    Properly paired: ${percentage}%"
                    print_empty_line
                fi
                ;;
            *"PASS:"*|*"WARNING:"*)
                if [[ "$line" =~ "WARNING:" ]]; then
                    message="${line}"
                    if [[ "$message" =~ "Low percentage of mapped reads" ]]; then
                        percentage=${message##*\(}
                        print_status "    QC Status: WARNING: Low percentage of mapped reads to CDSs"
                        print_status "              detected (${percentage}"
                    else
                        print_status "    QC Status: $message"
                    fi
                else
                    print_status "    QC Status: $line"
                fi
                print_empty_line
                ;;
            *'[1] "Annotating TUs of plus genes"'*|*"Annotating TUs of plus genes"*)
                print_status "Annotating TUs of plus genes ==>"
                print_empty_line
                ;;
            *'[1] "Annotating TUs of minus genes"'*|*"Annotating TUs of minus genes"*)
                print_status "Annotating TUs of minus genes <=="
                print_empty_line
                ;;
            *"Number of Excludons:"*)
                print_empty_line
                print_status "Number of Excludons"
                ;;
            *"Convergent: "*)
                print_status "            Convergent: ${line#*Convergent: }"
                ;;
            *"Divergent: "*)
                print_status "            Divergent: ${line#*Divergent: }"
                ;;
            *"The analysis is finished"*)
                print_footer
                print_section "Analysis Complete"
                print_status "All processes finished successfully"
                print_footer
                ;;
            *"ERROR:"*)
                print_section "Error"
                print_status "${line#ERROR: }"
                print_footer
                ;;
            "#####"*)
                ;;
            *)
                if [[ -n "$line" && ! "$line" =~ ^[[:space:]]*$ && 
                    ! "$line" =~ ^\[1\] && ! "$line" =~ ^"STARTING_ANNOTATION"$ ]]; then
                    print_status "$line"
                fi
                ;;
        esac
    done

    exit ${PIPESTATUS[0]}
}

# Execute main function with all arguments
main "$@"

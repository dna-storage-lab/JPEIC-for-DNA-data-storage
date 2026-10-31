#!/bin/bash
#==============================================================================
# Script Name: rapid_process_4000.sh
# Description: This script performs read-by-read decoding for 248nt single read
#              readout analysis on 4k pool data. It processes sequencing reads
#              individually to decode molecular barcodes and extract relevant
#              sequence information for downstream analysis.
# 
# Purpose:     The 248nt single read readout pipeline - executes
#              read-by-read decoding algorithm (version 3) to identify and
#              classify individual sequencing reads from the 4k pool dataset.
#
# Author:      [Qi Ge]
# Created:     [2025-7-4]
# Version:     3.0
#
# Usage:       ./rapid_process_4000.sh
#
#==============================================================================
echo "Begin Processing time :" `date`
ulimit -s 655350
# Record script start time
start_time=$(date +%s)

ideal_reference=./configure_file

DEFAULT_TARGET_DIR_BASE="../ONT_FASTQ"
DATE_PATTERN=""
TARGET_SUBDIR=""

MAX_FILES=20000
CHECK_INTERVAL=3
COMPLETE_CHECK_TIMES=2
COMPLETE_CHECK_INTERVAL=2
STACK_SIZE=655350

EXPECTED_READS=100

THREAD_COUNT=20

output_result_path="../Processing_result"
if [ ! -d "$output_result_path" ]; then
	mkdir -p "$output_result_path"
fi

# Function to compile all necessary programs
compile_programs() 
{
    echo -e "\033[32m=========== Compiling programs... ===========\033[0m"
    
    #step1
    g++ -o primer_alignment_nanopore_online primer_alignment_nanopore_online.c edlib.cpp -lpthread -lm
    check_status $? "g++ primer_alignment_nanopore_online compilation" "/dev/null"

    #step2
    gcc -c alloc.c -o alloc.o
    check_status $? "gcc alloc.c compilation" "/dev/null"

    gcc -c check.c -o check.o
    check_status $? "gcc check.c compilation" "/dev/null"

    gcc -c intio.c -o intio.o
    check_status $? "gcc intio.c compilation" "/dev/null"

    gcc -c open.c -o open.o
    check_status $? "gcc open.c compilation" "/dev/null"

    gcc -c mod2sparse.c -o mod2sparse.o
    check_status $? "gcc mod2sparse.c compilation" "/dev/null"

    gcc -c mod2dense.c -o mod2dense.o
    check_status $? "gcc mod2dense.c compilation" "/dev/null"

    gcc -c mod2convert.c -o mod2convert.o
    check_status $? "gcc mod2convert.c compilation" "/dev/null"

    g++ -c ldpc_decode.cpp -o ldpc_decode.o
    check_status $? "g++ ldpc_decode.cpp compilation" "/dev/null"

    gcc -c edlib.cpp -o edlib.o
    check_status $? "gcc edlib.cpp compilation" "/dev/null"

    g++  indel_correcting_and_LDPC_decoding.o  alloc.o check.o intio.o open.o mod2sparse.o mod2dense.o mod2convert.o ldpc_decode.o f_b_algorithm.o edlib.o -o indel_correcting_and_LDPC_decoding -std=c++11 -lpthread -lm
    check_status $? "g++ indel_correcting_and_LDPC_decoding linking" "/dev/null"

    #step3
    gcc -o index_identification_and_majority_voting_online index_identification_and_majority_voting_online.c -lm
    check_status $? "gcc index_identification_and_majority_voting_online compilation" "/dev/null"

    #step4
    g++ -o decode_RS decode_RS.cpp -lm
    check_status $? "gcc decode_RS compilation" "/dev/null"
    g++ -o verification verification.cpp -lm
    check_status $? "gcc verification compilation" "/dev/null"

    echo -e "\033[32m[SUCCESS] All programs compiled successfully\033[0m"
}

# Error handling function
check_status() 
{
    local exit_code=$1
    local program_name="$2"
    local log_file="$3"
    
    if [ $exit_code -ne 0 ]; then
        echo -e "\033[31m[ERROR] $program_name failed with exit code $exit_code\033[0m"
        echo "[$(date)] ERROR: $program_name failed with exit code $exit_code" >> "$log_file"
        echo "Stopping execution due to error in $program_name"
        exit $exit_code
    # else
    #     echo -e "\033[32m[SUCCESS] $program_name completed successfully\033[0m"
    #     echo "[$(date)] SUCCESS: $program_name completed successfully" >> "$log_file"
    fi
}

# Check if file exists and is not empty
check_file_exists() 
{
    local file_path="$1"
    local program_name="$2"
    local log_file="$3"
    
    if [ ! -f "$file_path" ]; then
        echo -e "\033[31m[ERROR] Expected output file $file_path does not exist after $program_name\033[0m"
        echo "[$(date)] ERROR: Expected output file $file_path does not exist after $program_name" >> "$log_file"
        exit 1
    elif [ ! -s "$file_path" ]; then
        echo -e "\033[33m[WARNING] Output file $file_path is empty after $program_name\033[0m"
        echo "[$(date)] WARNING: Output file $file_path is empty after $program_name" >> "$log_file"
    else
        echo -e "\033[32m[SUCCESS] Output file $file_path created successfully\033[0m"
    fi
}

# Function to check if primer decoding has been completed
check_primer_completed() 
{
    # Check completion status of single primer
    if [[ "$primer_completed" == "true" ]]; then
        return 0  # Primer completed
    fi
    return 1  # Primer not completed
}

# Function to create primer completion status file
save_primer_completion_status() {
    local output_dir="$1"
    
    local completion_flag_file="${output_dir}/primer_completed.flag"
    echo "completed" > "$completion_flag_file"
    echo "[$(date)] Primer marked as completed" >> "${output_dir}/execution_error.log"
}

# Function to load primer completion status from file
load_primer_completion_status() {
    local output_dir="$1"
    
    # Reset primer status to not completed
    primer_completed=false
    
    # Check primer completion status file
    local completion_flag_file="${output_dir}/primer_completed.flag"
    if [[ -f "$completion_flag_file" ]]; then
        local content=$(tr -d '[:space:]' < "$completion_flag_file")
        if [[ "$content" == "completed" ]]; then
            primer_completed=true
            echo -e "\033[32m[INFO] Primer was previously completed\033[0m"
        fi
    fi
}

# Check if file is complete by monitoring size stability
is_file_complete()
{
    local filename="$1"
    local check_times=${2:-$COMPLETE_CHECK_TIMES}
    local interval=${3:-$COMPLETE_CHECK_INTERVAL}

    local prev_size=-1
    local current_size

    for ((i=0; i<check_times; i++)); do
        if [ ! -f "$filename" ]; then
            return 1
        fi

        current_size=$(stat -c%s "$filename" 2>/dev/null) || return 1

        if [ "$prev_size" -eq "$current_size" ]; then
            return 0
        fi

        prev_size=$current_size

        if [ $i -lt $((check_times-1)) ]; then
            sleep $interval
        fi
    done

    return 1
}

is_fastq_file_complete_by_lines() 
{
    local filename="$1"
    local check_times=${2:-$COMPLETE_CHECK_TIMES}
    local interval=${3:-$COMPLETE_CHECK_INTERVAL}
    local expected_reads="${4:-$EXPECTED_READS}"

    local expected_lines=$((expected_reads * 4))
    local line_count

    for ((i=0; i<check_times; i++)); 
    do
        if [ ! -f "$filename" ]; then
            return 1
        fi

        line_count=$(wc -l < "$filename")

        if [ "$line_count" -eq "$expected_lines" ]; then
            return 0
        else
            echo "[INFO] $filename currently has $line_count lines (expecting $expected_lines), waiting..."
            sleep "$interval"
        fi
    done

    echo "[ERROR] $filename did not reach expected line count after $check_times checks."
    return 1
}

# Wait for file existence and completeness
wait_for_file()
{
    local target_dir="$1"
    local file_index="$2"
    
    # Find the file with pattern *_${file_index}.fastq
    local pattern="${target_dir}/*_${file_index}.fastq"
    local filename=""
    
    echo "Waiting for file with pattern: $pattern"
    
    # Wait for file to exist
    while [ -z "$filename" ]; do
        filename=$(ls $pattern 2>/dev/null | head -1)
        if [ -z "$filename" ]; then
            echo "File matching pattern $pattern does not exist, waiting..."
            sleep $CHECK_INTERVAL
        fi
    done
    
    # Set global variable for use in main function
    current_file="$filename"
    
    echo "File $filename exists, checking completeness..."

    # while ! is_file_complete "$filename"; do
    #     echo "File $filename is still being written, waiting for completion..."
    #     sleep $CHECK_INTERVAL
    # done
    while ! is_fastq_file_complete_by_lines "$filename"; do
        echo "File $filename is still being written, waiting for completion..."
        sleep $CHECK_INTERVAL
    done

    echo "File $filename is complete"
}
# Function to search for directory with date pattern
find_target_directory() 
{
    local base_dir="$1"
    local date_pattern="$2"
    local subdir="$3"
    
    echo "Searching for directory with pattern: ${base_dir}/${date_pattern}/${subdir}" >&2  # Output to stderr
    
    while true; do
        if [ ! -d "$base_dir" ]; then
            echo "Base directory $base_dir does not exist, waiting..." >&2
            sleep $CHECK_INTERVAL
            continue
        fi
        
        local found_dirs=($(find "$base_dir" -maxdepth 1 -type d -name "$date_pattern" 2>/dev/null))
        
        if [ ${#found_dirs[@]} -eq 0 ]; then
            echo "No directories found matching pattern $date_pattern in $base_dir, waiting..." >&2
            sleep $CHECK_INTERVAL
            continue
        fi
        
        for dir in "${found_dirs[@]}"; do
            local full_path="${dir}/${subdir}"
            if [ -d "$full_path" ]; then
                echo "Found target directory: $full_path" >&2
                printf "%s" "$full_path"
                return 0
            fi
        done
        
        echo "Found directories matching date pattern, but none contain $subdir subdirectory, waiting..." >&2
        sleep $CHECK_INTERVAL
    done
}

# Main loop
main()
{
    
    # First compile all necessary programs
    compile_programs
    
    # Find the actual target directory
    echo "Starting directory search..."
    local target_dir=${DEFAULT_TARGET_DIR_BASE}
    
    if [ -z "$target_dir" ]; then
        echo "Error: Could not find target directory"
        exit 1
    fi
    
    echo "Using target directory: $target_dir"

    
    local max_files=${1:-$MAX_FILES}
    local file_index=0
    local thread_setting=${THREAD_COUNT}
    local expected_reads=${EXPECTED_READS}

    local output_path=${output_result_path}
    result_path=${output_path}/result
    if [ ! -d "$result_path" ]; then
        mkdir -p "$result_path"
    fi
    analysis_path=${result_path}/analysis_result
    if [ ! -d "$analysis_path" ]; then
        mkdir -p "$analysis_path"
    fi
    erasure_record_path=${result_path}/erasure_record
    if [ ! -d "$erasure_record_path" ]; then
        mkdir -p "$erasure_record_path"
    fi
    Temporary=${result_path}/Temporary
    if [ ! -d "$Temporary" ]; then
        mkdir -p "$Temporary"
    fi
    decoding_path=${result_path}/decoding_result
    if [ ! -d "$decoding_path" ]; then
        mkdir -p "$decoding_path"
    fi
    pic_recovery_path=${result_path}/picture_recovery
    if [ ! -d "$pic_recovery_path" ]; then
        mkdir -p "$pic_recovery_path"
    fi
    # Create error log file
    error_log_file=${result_path}/execution_error.log
    echo "=== Execution Log Started at $(date) ===" > $error_log_file

    time_log_file=${result_path}/group_processing_time_fast.log
    echo "Group,Start Time,End Time,Processing Time(s),Cumulative Time(s),Date" > $time_log_file

    #======================================= main function ========================================
    loop_start_time=$(date +%s)
    
    while [ ! -d "$target_dir" ]; do
        echo -e "\033[31mWaiting for directory '${target_dir}' to be created...\033[0m"
        sleep $CHECK_INTERVAL
    done
    echo -e "\033[32mDirectory '${target_dir}' exists.\033[0m"

    echo "Starting to detect FASTQ files..."
    echo "Target directory: ${target_dir}"
   

    while [ $file_index -lt $max_files ];
    do

        echo -e "\n================================================="
        echo "Current file_index: $file_index"
        
        wait_for_file "$target_dir" "$file_index"
        echo -e "======\033[32m[$file_index] Starting to process: $(basename "$current_file")======\033[0m"
        

        primer_completed=false
        load_primer_completion_status "$result_path"

        group_start_time=$(date +%s)
        # Check if primer has been completed before starting each group
        if check_primer_completed; then
            echo -e "\033[32m[SUCCESS] Primer has been successfully decoded! Stopping processing.\033[0m"
            echo "[$(date)] SUCCESS: Primer completed - stopping at group $s" >> "$error_log_file"
            break
        fi
        #====== step1 primer_alignment_nanopore_online ========================
        echo -e "\033[32m==== Step1: Primer alignment and identification ====\033[0m"
        primer_start_time=$(date +%s)
        fastq_file=${current_file}
        primer_list=${ideal_reference}/primer_list.txt
        primer_num=2
        primer_len=20
        temp_output_prefix=${result_path}
        min_length=248
        reads_primer_match_counts=${result_path}/primer_match_counts_${file_index}.txt
        result_output_prefix=${result_path}
        thread=${thread_setting}
        valid_reads_record=${result_path}/valid_reads_record_${file_index}.txt
        payload_size=208
        valid_bais_size=10
        file_index_number=${file_index}
        
        ./primer_alignment_nanopore_online $fastq_file \
                                $primer_list \
                                $primer_num \
                                $primer_len \
                                $temp_output_prefix \
                                $min_length \
                                $reads_primer_match_counts \
                                $result_output_prefix \
                                $thread \
                                $valid_reads_record \
                                $payload_size \
                                $valid_bais_size \
                                $file_index_number
        check_status $? "primer_alignment_nanopore_online for file $file_index" "$error_log_file"
        
        primer_end_time=$(date +%s)
        primer_time=$((primer_end_time - primer_start_time))
        echo -e "[FASTQ: ${file_index}] primer_alignment completed in ${primer_time} seconds"

        #====== step2 fb_correct_ldpc_208nt_online ========================
        echo -e "\033[32m====Step2: Forward-backward ecc and LDPC decoding ====\033[0m"
        fb_start_time=$(date +%s)
        Reads_without_primer=${result_path}/primer_pair_1_${file_index}.fasta
        
        # Calculate primer-identified reads count and coverage
        if [[ -f "$Reads_without_primer" ]]; then
            primer_reads_count=$(grep -c "^>" "$Reads_without_primer" 2>/dev/null || echo 0)
            coverage=$(echo "scale=4; $primer_reads_count / 4000" | bc)
            cumulative_reads=$((cumulative_reads + primer_reads_count))
            cumulative_coverage=$(echo "scale=2; $cumulative_reads / 4000" | bc | sed 's/^\./0./')
            echo "Cumulative Coverage: $cumulative_reads reads / 4000 total (${cumulative_coverage}x)"
        else
            echo "Warning: $Reads_without_primer not found, cannot calculate coverage"
        fi
        # check_file_exists "$Reads_without_primer" "primer_alignment_nanopore_online primer_pair_1" "$error_log_file"
        
        water_file=${ideal_reference}/SequenceLimitRunLength
        LDPC_output=${result_path}/ldpc_decoded_bits_${file_index}.txt
        output_prefix=${result_path}
        summary_error=${analysis_path}/fb_ldpc_summary_error_${file_index}.txt
        ldpc_failed_read=${result_path}/LDPC_failed_reads_${file_index}.fasta
        p_ins=0.005
        p_del=0.008
        p_sub=0.005
        thread=${thread_setting}
        ./indel_correcting_and_LDPC_decoding "$Reads_without_primer" \
                    "$water_file" \
                    "$LDPC_output" \
                    "$output_prefix" \
                    "$summary_error" \
                    "$ldpc_failed_read" \
                    "$p_ins" \
                    "$p_del" \
                    "$p_sub" \
                    "$thread" 

        check_status $? "indel_correcting_and_LDPC_decoding for file $file_index" "$error_log_file"


        fb_end_time=$(date +%s)
        fb_time=$((fb_end_time - fb_start_time))
        echo "[${file_index}] Indel error correction completed in ${fb_time} seconds"
        
        #=============== step3 address_parse_and_majority_voting ========================
        echo -e "\033[32m==== Step3: Index identification and majority voting ====\033[0m"
        # Parsing address to determine whether to decode
        address_start_time=$(date +%s)

        if [ $file_index -eq 0 ]; then
            prev_s=$file_index
        else
            prev_s=$((file_index-1))
        fi
        decoded_filename=${result_path}/ldpc_decoded_bits_${file_index}.txt
        ifdecoding_flag=${result_path}/ifdecoding_flag_${file_index}.txt
        old_found_codeword_file=${result_path}/found_codeword_${prev_s}.txt
        address_info=${erasure_record_path}/address_info_${file_index}.txt
        total_erasure_record=${erasure_record_path}/total_erasure_record_${file_index}.txt
        new_found_codeword_file=${result_path}/found_codeword_${file_index}.txt
        result_for_rs_decoding=${result_path}/result_for_rs_decoding_${file_index}.txt

        ./index_identification_and_majority_voting_online \
                                "$decoded_filename" \
                                "$ifdecoding_flag" \
                                "$old_found_codeword_file" \
                                "$address_info" \
                                "$total_erasure_record" \
                                "$new_found_codeword_file" \
                                "$result_for_rs_decoding"

        check_status $? "index_identification_and_majority_voting_online for file $file_index" "$error_log_file"

        address_end_time=$(date +%s)
        address_time=$((address_end_time - address_start_time))
        echo -e "[FASTQ: ${file_index}] index parsing and majority voting completed in ${address_time} seconds"

        echo -e "\033[32m==== Step4: RS decoding ====\033[0m"
        #=============== step4 decode_RS ========================
        # echo "=========== Checking if decoding is possible... ==========="  
        ifdecoding_flag_file=${result_path}/ifdecoding_flag_${file_index}.txt
        if [[ ! -f "$ifdecoding_flag_file" ]]; then
            echo "Error: ${ifdecoding_flag_file} does not exist"
        else
            content=$(tr -d '[:space:]' < "$ifdecoding_flag_file")

            if [[ "$content" == "yes" ]]; then
                echo -e "\033[32mDecoding is possible. Starting RS decoding and recovery\033[0m"

                result_for_rs_decoding=${result_path}/result_for_rs_decoding_${file_index}.txt
                rs_decoding_for_recovery=${decoding_path}/rs_decoding_for_recovery_${file_index}.txt
                scramble_sequence=${ideal_reference}/scramble_sequence.txt
               # check_file_exists "$scramble_sequence" "decode_RS" "$error_log_file"

                success_decoding_flag=${decoding_path}/success_decoding_${file_index}.txt

                ./decode_RS \
                    "$result_for_rs_decoding" \
                    "$rs_decoding_for_recovery" \
                    "$scramble_sequence" \
                    "$success_decoding_flag"

                check_status $? "decode_RS for file $file_index" "$error_log_file"
                # check_file_exists "$rs_decoding_for_recovery" "decode_RS" "$error_log_file"
                # check_file_exists "$success_decoding_flag" "decode_RS" "$error_log_file"

                # echo "=========== Checking if RS decoding is success... ==========="

                content_rs=$(tr -d '[:space:]' < "$success_decoding_flag")

                if [[ "$content_rs" == "yes" ]]; then
                    echo -e "\033[32m===== [FASTQ: $file_index] RS decoding successful. =====\033[0m"

                    primer_completed=true
                    save_primer_completion_status "$result_path"

                    inputStream=${rs_decoding_for_recovery}
                    filename=${pic_recovery_path}/recovered_pic3.jpg
                    pictureSize=70033

                    # Picture recovery verification
                    ./verification "$inputStream" "$filename" "$pictureSize"
                    check_status $? "verification for file $file_index_number" "$error_log_file"

                    # Open picture
                    xdg-open "$filename" 2>/dev/null
                    sleep 0.5s
                    # Set picture position to bottom right corner of screen
                    wmctrl -r "$(basename "$filename")" -e 0,1320,680,690,517

                    echo "[FASTQ: $file_index] Completed RS decoding and recovery"
                else
                    echo -e "\033[33m======[FASTQ: $file_index] Decoding failed for this file. Moving to the next file.======\033[0m"
                fi
            else
                echo -e "\033[33mDecoding is not possible. Skipping RS decoding and recovery.\033[0m"
            fi
        fi

        # Calculate total time for this group
        group_end_time=$(date +%s)
        group_total_time=$((group_end_time - group_start_time))
        cumulative_time=$((group_end_time - loop_start_time))
        echo "======[FASTQ: $file_index] Processing completed — Time: ${group_total_time} seconds======"
        echo "======[${cumulative_coverage}x] Cumulative Time: ${cumulative_time} seconds======"
        # Write time information to log file
        echo "$file_index,$(date -d @$group_start_time +'%H:%M:%S'),$(date -d @$group_end_time +'%H:%M:%S'),$group_total_time,$cumulative_time,$(date +'%Y-%m-%d')" >> $time_log_file
        
        if [[ "$primer_completed" == "true" ]]; then
            echo -e "\033[33m======Readout already completed!!! Used $((file_index+1)) files======\033[0m"
            echo -e "\033[33m======Sequencing coverage used for data readout: ${cumulative_coverage}x======\033[0m"
            end_time=$(date +%s)
            runtime=$((end_time - start_time))
            hours=$((runtime / 3600))
            minutes=$(((runtime % 3600) / 60))
            seconds=$((runtime % 60))
            echo "----------------------------------------"
            echo "Total runtime: ${hours}h ${minutes}m ${seconds}s (${runtime}s)"
            echo "End Processing time: $(date)"
            echo "----------------------------------------"
            # Write total time to log file
            echo "Total runtime: ${hours}h ${minutes}m ${seconds}s (${runtime}s)" >> $time_log_file
            break
        # else
        #     echo "Decoding failed, continuing to process the next file"
            
        fi

        file_index=$((file_index + 1))

    done
}
    main "$@"


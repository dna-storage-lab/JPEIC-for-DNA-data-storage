#!/bin/bash
#==============================================================================
# Script Name: run_video_readout_NGS.sh
# Description: Offline decoding pipeline for NGS paired-end assembled FASTQ.
#
# Original Author: [Qi Ge]
# Modified:        [2026-01-26]
# Version:         1.0 (NGS Offline)
#
# Usage:           ./run_video_readout_NGS.sh
#==============================================================================

echo "Begin Processing time : $(date '+%F %T')"
ulimit -s 655350

# ================== Input Configuration ==================
# Hardcoded input file path
INPUT_FASTQ="./Assembled_fastq/assembled_test.fastq"

if [ ! -f "$INPUT_FASTQ" ]; then
    echo -e "\033[31m[ERROR] Input file not found: $INPUT_FASTQ\033[0m"
    echo "Please ensure 'Assembled_fastq/assembled_test.fastq' exists."
    exit 1
fi

# ================== Basic Configuration ==================
THREAD_SETTING=72
output_result_path=./NGS_recovery_results
configuration_files=./configuration_files

# Create Output Directories
output_dir=${output_result_path}/result
bch_result_path=${output_dir}/bch_result
erasure_record_path=${output_dir}/erasure_record

mkdir -p "$output_dir" "$bch_result_path" "$erasure_record_path"

error_log_file=${output_dir}/execution_error.log
time_log_file=${output_dir}/processing_time.log

# ================== Helper Functions ==================
check_status()
{
    local exit_code=$1
    local program_name="$2"
    local log_file="$3"
    
    if [ $exit_code -ne 0 ]; then
        echo -e "\033[31m[ERROR] $program_name failed (exit code $exit_code)\033[0m"
        [ -n "$log_file" ] && echo "[$(date)] ERROR: $program_name failed with exit code $exit_code" >> "$log_file"
        exit $exit_code
    fi
}

check_file_exists() 
{
    local file_path="$1"
    local program_name="$2"
    local log_file="$3"
    
    if [ ! -f "$file_path" ]; then
        echo -e "\033[31m[ERROR] Output file $file_path missing after $program_name\033[0m"
        [ -n "$log_file" ] && echo "[$(date)] ERROR: Missing $file_path after $program_name" >> "$log_file"
        exit 1
    elif [ ! -s "$file_path" ]; then
        echo -e "\033[33m[WARN] Output file $file_path is empty after $program_name\033[0m"
        [ -n "$log_file" ] && echo "[$(date)] WARNING: Empty $file_path after $program_name" >> "$log_file"
    fi
}

# ================== Main Process ==================

echo
echo "============================================================"
echo "[Init] NGS Readout Pipeline started at $(date '+%F %T')"
echo "       Input File: $INPUT_FASTQ"
echo "============================================================"

start_time=$(date +%s)

# ================== Step 1: Primer alignment ==================
echo "[Step 1] Primer alignment..."
primer_start_time=$(date +%s)

# Configuration for Step 1
primer_list=${configuration_files}/primer_list_10nt.txt
primer_num=2
payload_size=130
oligo_len=170
valid_bais_size=5

# Outputs for Step 1
reads_primer_match_counts=${output_dir}/primer_match_counts.txt
result_output=${output_dir}/primer_alignment_result.fasta
valid_reads_record=${output_dir}/valid_reads_record.txt
length_record_file=${output_dir}/reads_length_record.txt

# Run Command
./ngs_primer_align "$INPUT_FASTQ" \
                   "$primer_list" \
                   "$primer_num" \
                   "$payload_size" \
                   "$reads_primer_match_counts" \
                   "$result_output" \
                   "$valid_reads_record" \
                   "$oligo_len" \
                   "$valid_bais_size" \
                   "$length_record_file" \
                   "$THREAD_SETTING"

check_status $? "ngs_primer_align" "$error_log_file"
check_file_exists "$result_output" "ngs_primer_align" "$error_log_file"

primer_end_time=$(date +%s)
echo "[OK] Step 1 finished ($((primer_end_time - primer_start_time))s)"

# ================== Step 2: FB + BCH decoding ==================
echo "[Step 2] FB + BCH decoding..."
fb_start_time=$(date +%s)

# Inputs/Configs for Decoding
# NOTE: Using Step 1 output directly as input
selectReads=${result_output}
water_file=${configuration_files}/SequenceLimitRunLength
p_ins=0.005
p_del=0.008
p_sub=0.005

# Outputs for Decoding
hanmming_file=${bch_result_path}/fb_error_record.txt
BCH_output=${bch_result_path}/BCH_decoded_bits.txt
summary_error=${bch_result_path}/fb_bch_summary_error.txt

# Run Command
./indel_correction_BCH_decoding "$selectReads" \
                                "$water_file" \
                                "$hanmming_file" \
                                "$BCH_output" \
                                "$summary_error" \
                                "$p_ins" \
                                "$p_del" \
                                "$p_sub" \
                                "$THREAD_SETTING"

check_status $? "indel_correction_BCH_decoding" "$error_log_file"
check_file_exists "$BCH_output" "indel_correction_BCH_decoding" "$error_log_file"

fb_end_time=$(date +%s)
echo "[OK] Step 2 finished ($((fb_end_time - fb_start_time))s)"

# ================== Step 3: Majority voting and Deinterleave ==================
echo "[Step 3] Majority voting and Deinterleave..."
deinterleave_start_time=$(date +%s)

# Using BCH output from Step 2 directly
merged_bch_output=${BCH_output}
interleave_order_file=${configuration_files}/rs4000_250group_diagonal_interleave_order
deinterleaved_output_dir=${output_dir}/codeword_for_rs_decoding
mkdir -p "$deinterleaved_output_dir"

flag_file_for_rs=${deinterleaved_output_dir}/RS_decoding_flag.txt
deinterleaved_erasure_record=${deinterleaved_output_dir}/deinterleaved_erasure_record.txt
bch_summary_file=${output_dir}/bch_dropout_summary.txt

./majority_voting_deinterleave \
    "$merged_bch_output" \
    "$interleave_order_file" \
    "$deinterleaved_output_dir" \
    "$flag_file_for_rs" \
    "$deinterleaved_erasure_record" \
    "$bch_summary_file"

check_status $? "majority_voting_deinterleave" "$error_log_file"
check_file_exists "$flag_file_for_rs" "majority_voting_deinterleave" "$error_log_file"

deinterleave_end_time=$(date +%s)
echo "[OK] Majority voting and Deinterleave finished ($((deinterleave_end_time - deinterleave_start_time))s)"

# ================== Step 4: RS Decoding & Video Recovery ==================
# Check if majority voting and deinterleaving considers the data sufficient for RS decoding
content_deinter=$(tr -d '[:space:]' < "$flag_file_for_rs")

if [[ "$content_deinter" == "yes" ]]; then
    echo "[Step 4] RS decoding & video recovery..."
    rsdecoding_start_time=$(date +%s)

    rs_decoding_for_recovery=${output_dir}/rs_decoding_for_recovery
    mkdir -p "$rs_decoding_for_recovery"

    status_file=${rs_decoding_for_recovery}/status_file.txt
    success_decoding_flag=${rs_decoding_for_recovery}/success_decoding.txt
    bits_for_video_recovery=${rs_decoding_for_recovery}/bits_for_video_recovery.txt

    # Run RS Decoding
    ./rs_decode_parallel \
        "$deinterleaved_output_dir" \
        "$rs_decoding_for_recovery" \
        "$THREAD_SETTING" \
        "$status_file" \
        "$success_decoding_flag" \
        "$bits_for_video_recovery"

    check_status $? "rs_decode_parallel" "$error_log_file"

    rsdecoding_end_time=$(date +%s)
    echo "[OK] RS decoding step finished ($((rsdecoding_end_time - rsdecoding_start_time))s)"

    # Check if RS decoding was actually successful
    content_rs=$(tr -d '[:space:]' < "$success_decoding_flag")
    
    if [[ "$content_rs" == "yes" ]]; then
        echo "[Step 5] Recovering Video..."
        
        recovery_video_path=${output_dir}/video_recovery
        mkdir -p "$recovery_video_path"
        
        inputStream=${bits_for_video_recovery}
        scramble_code=${configuration_files}/watermark_448800bit.txt
        video_recovery=${recovery_video_path}/recovered_video.mp4
        
        # Run Video Recovery
        ./recover_video_from_bitstream "$scramble_code" "$inputStream" "$video_recovery"
        check_status $? "recover_video_from_bitstream" "$error_log_file"

        echo -e "\033[32m[SUCCESS] Video recovered at: $video_recovery\033[0m"
        
        # Optional: Open video if in GUI environment
        if command -v xdg-open &> /dev/null; then
             xdg-open "$video_recovery" 2>/dev/null
        fi
    else
        echo -e "\033[31m[FAIL] RS decoding failed. Too many erasures or errors.\033[0m"
    fi
else
    echo -e "\033[33m[WARN] Data insufficient or erasure rate too high after deinterleaving. Skipping RS decoding.\033[0m"
fi

# ================== Summary ==================
end_time=$(date +%s)
total_time=$((end_time - start_time))

echo
echo "============================================================"
echo "[Done] Pipeline completed in ${total_time}s ($(date '+%F %T'))"
echo "============================================================"

# Write time log
echo "Input_File,Start_Time,End_Time,Total_Time(s),Date" > "$time_log_file"
echo "$INPUT_FASTQ,$(date -d @$start_time +'%H:%M:%S'),$(date -d @$end_time +'%H:%M:%S'),$total_time,$(date +'%Y-%m-%d')" >> "$time_log_file"
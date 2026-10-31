#!/bin/bash
#==============================================================================
# Script Name: run_video_rapid_recovery.sh
# Description: Online decoding pipeline for 170nt single-read readout (1M pool)
#              [Two-chip version: process one FASTQ at a time from 2 chips]
# 
# Author:      [Qi Ge]
# Created:     [2025-12-14]
# Version:     2.0 (Two-chip, single fastq per group, dynamic scheduling)
#
# Usage:       ./run_video_rapid_recovery.sh
#
#==============================================================================

echo "Begin Processing time : $(date '+%F %T')"
ulimit -s 655350

# ================== Basic Configuration ==================
DEFAULT_TARGET_DIR_BASE="./ECS-Pool2_nanopore_data"

MAX_FILES=2000000
CHECK_INTERVAL=3

EXPECTED_READS=50000 # Target number of reads per fastq file
THREAD_SETTING=72

TOTAL_EXPECTED_SEQS=1000000  # Original 1 million sequences
output_result_path=./1M_real_time_results

[ ! -d "$output_result_path" ] && mkdir -p "$output_result_path"


LOG_FILE="${output_result_path}/terminal_full_record.log"
mkdir -p "$(dirname "$LOG_FILE")"
exec > >(tee -a "$LOG_FILE") 2>&1
# ======================================================

configuration_files=../NGS_Decoder_ECS_pool2/configuration_files

cumulative_reads=${cumulative_reads:-0}
cumulative_bases=${cumulative_bases:-0}

# Populated at runtime: fastq_pass directories for the two chips
CHIP1_FASTQ_DIR=""
CHIP2_FASTQ_DIR=""

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

# ================== Main Process ==================
main()
{
    echo
    echo "============================================================"
    echo "[Init] Rapid readout started at $(date '+%F %T')"
    echo "============================================================"

    local max_files=${1:-$MAX_FILES}
    local s=0

    local result_path=${output_result_path}
    rs_success=0    
    # Directory preparation
    output_dir=${result_path}/result
    Temporary=${output_dir}/Temporary
    reconstruction_path=${output_dir}/reconstructed_refs
    analysis_path=${output_dir}/analysis_result
    bch_result_path=${output_dir}/bch_result
    erasure_record_path=${output_dir}/erasure_record

    mkdir -p "$output_dir" "$Temporary" "$reconstruction_path" \
             "$analysis_path" "$bch_result_path" "$erasure_record_path"

    error_log_file=${output_dir}/execution_error.log
    time_log_file=${output_dir}/group_processing_time.log
    coverage_log="${output_dir}/coverage_report.txt"
    echo "Group,Coverage,Start Time,End Time,Processing Time(s),Cumulative Time(s),Date" > "$time_log_file"
    echo -e "File_Index\tCurrent_Reads\tCumulative_Reads\tCumulative_Bases\tCumulative_Coverage" > "$coverage_log"

    loop_start_time=0
    total_primer_time=0
    total_compare_time=0
    total_fb_time=0
    total_bch_comb_time=0
    total_deinterleave_time=0
    total_rsdecoding_time=0
    total_video_readout_time=0

    # ================== Main loop (process one FASTQ per group) ==================
    while [ $s -lt $max_files ]; do
        #group_start_time=$(date +%s)

        # ---------- Group header information ----------
        echo
        echo "============================================================"
        echo -e "\033[35m[Group $s] Start at $(date '+%H:%M:%S')\033[0m"
        echo "  Next index: $s"
        echo "============================================================"

        # -------- Select the current chip and fastq files to be processed --------
        fastq_file="${DEFAULT_TARGET_DIR_BASE}/sub_${s}.fastq"
        
        if [[ "$loop_start_time" -eq 0 ]]; then
            loop_start_time=$(date +%s)
        fi
        group_start_time=$(date +%s)
        
        # ================== Step 1: Primer alignment ==================
        echo "[Step 1] Primer alignment..."
        primer_start_time=$(date +%s)

        primer_list=${configuration_files}/primer_list_10nt.txt
        primer_num=2
        payload_size=130
        reads_primer_match_counts=${output_dir}/primer_match_counts_${s}.txt
        result_output=${output_dir}/primer_alignment_result_${s}.fasta
        valid_reads_record=${output_dir}/valid_reads_record_${s}.txt
        oligo_len=170
        valid_bais_size=5
        length_record_file=${output_dir}/reads_length_record_${s}.txt
        thread=${THREAD_SETTING}

        ./ont_primer_align "$fastq_file" \
                                "$primer_list" \
                                "$primer_num" \
                                "$payload_size" \
                                "$reads_primer_match_counts" \
                                "$result_output" \
                                "$valid_reads_record" \
                                "$oligo_len" \
                                "$valid_bais_size" \
                                "$length_record_file" \
                                "$thread"
        
        primer_end_time=$(date +%s)
        primer_time=$((primer_end_time - primer_start_time))
        total_primer_time=$((total_primer_time + primer_time))
        echo "[OK] Step 1 finished (${primer_time}s)"

        Reads_without_primer=${result_output}
        primer_reads_count=0
        total_bases_step=0

        if [[ -f "$Reads_without_primer" ]]; then
            read -r primer_reads_count total_bases_step <<< $(awk '/^>/ {seqs++} !/^>/ {bases += length($0)} END {print seqs+0, bases+0}' "$Reads_without_primer")
            primer_reads_count=${primer_reads_count:-0}
            total_bases_step=${total_bases_step:-0}

            cumulative_reads=$((cumulative_reads + primer_reads_count))
            cumulative_bases=$((cumulative_bases + total_bases_step))

            total_expected_bases=$(($TOTAL_EXPECTED_SEQS * 130))
            total_expected_bases=$(($TOTAL_EXPECTED_SEQS * 130))
            if [ "${total_expected_bases:-0}" -gt 0 ]; then
                cumulative_coverage=$(echo "scale=4; $cumulative_bases / $total_expected_bases" | bc | awk '{printf "%.4f", $0}')
            else
                cumulative_coverage="0.0000"
            fi

            echo -e "${s}\t${primer_reads_count}\t${cumulative_reads}\t${cumulative_bases}\t${cumulative_coverage}" >> "$coverage_log"
        else
            echo "  [WARN] $Reads_without_primer not found."
            cumulative_coverage=$(echo "scale=4; $cumulative_bases / $total_expected_bases" | bc | awk '{printf "%.4f", $0}')
        fi

        # ================== Step 2: Compare with regenerative references ==================
        echo "[Step 2] Filtering reads with regenerative references..."
        compare_start_time=$(date +%s)

        if [ $s -eq 0 ]; then
            prev_s=$s
            selectReads=${output_dir}/primer_alignment_result_${s}.fasta
        else
            prev_s=$((s-1))
            
            Reads_without_primer=${output_dir}/primer_alignment_result_${s}.fasta
            selectReads=${output_dir}/subfile_${s}_selectreads.fasta
            log=${erasure_record_path}/subfile_${s}_reads_select_log.txt

            if [ ! -s "$Reads_without_primer" ]; then
                echo "  [INFO] $Reads_without_primer is empty; skipping minimap2 and save_unmapped_reads."
                : > "$selectReads" 
            else
                found_oligopool=${reconstruction_path}/cumulative_refs_${prev_s}.fasta
                fulsam=${output_dir}/subfile_${s}.sam

                minimap2_start_time=$(date +%s)
                ./minimap2 -ax sr -t 60 "${found_oligopool}" "$Reads_without_primer" > "$fulsam" 2> /dev/null

                ./extract_unmapped_reads "${fulsam}" "${selectReads}" "${log}"
                minimap2_end_time=$(date +%s)
                minimap2_time=$((minimap2_end_time - minimap2_start_time))
                echo "[OK] minimap2 + filtering finished (${minimap2_time}s)"
            fi

            if [ ! -s "$selectReads" ]; then
                echo "  [WARN] $selectReads is empty (no new reads), continue to next group."
            fi
        fi

        compare_end_time=$(date +%s)
        compare_time=$((compare_end_time - compare_start_time))
        total_compare_time=$((total_compare_time + compare_time))
        echo "[OK] Step 2 finished (${compare_time}s)"

        # ================== Step 3: FB + BCH decoding ==================
        echo "[Step 3] FB + BCH decoding..."
        fb_start_time=$(date +%s)

        water_file=${configuration_files}/SequenceLimitRunLength
        hanmming_file=${bch_result_path}/fb_error_record_${s}.txt
        BCH_output=${bch_result_path}/BCH_decoded_bits_${s}.txt
        summary_error=${bch_result_path}/fb_bch_summary_error_${s}.txt
        recovered_refs=${reconstruction_path}/BCH_decoded_refs_${s}.fasta
        p_ins=0.005
        p_del=0.010
        p_sub=0.005
        thread=${THREAD_SETTING}
   
        ./indel_correction_BCH_decoding "$selectReads" \
                    "$water_file" \
                    "$hanmming_file" \
                    "$BCH_output" \
                    "$summary_error" \
                    "$recovered_refs" \
                    "$p_ins" \
                    "$p_del" \
                    "$p_sub" \
                    "$thread"

        fb_end_time=$(date +%s)
        fb_time=$((fb_end_time - fb_start_time))
        total_fb_time=$((total_fb_time + fb_time))
        echo "[OK] Step 3 finished (${fb_time}s)"

        # ================== Step 4: Updating BCH decoding results ==================
        echo "[Step 4] Updating BCH decoding results..."
        bch_comb_start_time=$(date +%s)
        echo "  -----------------------------------------"
        prev_s=$((s-1))
        old_bch_file=${bch_result_path}/Cumulative_BCH_decoded_bits_${prev_s}.txt
        new_bch_file=${bch_result_path}/BCH_decoded_bits_${s}.txt
        merged_bch_output=${bch_result_path}/Cumulative_BCH_decoded_bits_${s}.txt
        
        old_dna_file=${reconstruction_path}/cumulative_refs_${prev_s}.fasta
        new_dna_file=${reconstruction_path}/BCH_decoded_refs_${s}.fasta
        merged_dna_output=${reconstruction_path}/cumulative_refs_${s}.fasta
        flag_file=${bch_result_path}/RS_decoding_flag_${s}.txt
        rate_file=${bch_result_path}/dropout_rate_${s}.txt

        ./merge_bch_results "$old_bch_file" \
                                  "$new_bch_file" \
                                  "$merged_bch_output" \
                                  "$old_dna_file" \
                                  "$new_dna_file" \
                                  "$merged_dna_output" \
                                  "$flag_file" \
                                  "$rate_file"
        echo "  -----------------------------------------"
        bch_comb_end_time=$(date +%s)
        bch_comb_time=$((bch_comb_end_time - bch_comb_start_time))
        total_bch_comb_time=$((total_bch_comb_time + bch_comb_time))
        echo "[OK] Step 4 finished (${bch_comb_time}s)"

		rm -rf ${result_output}
		rm -rf ${selectReads}
		rm -rf ${fulsam}
		rm -rf ${old_bch_file}
		rm -rf ${new_bch_file}
		rm -rf ${old_dna_file}
		rm -rf ${new_dna_file}
		
        # ================== Step 5: Deinterleave & RS decoding ==================
        ifdecoding_flag_file=${bch_result_path}/RS_decoding_flag_${s}.txt
        if [[ ! -f "$ifdecoding_flag_file" ]]; then
            echo "[Info] RS_decoding_flag not found, skip RS for this group."
        else
            content=$(tr -d '[:space:]' < "$ifdecoding_flag_file")

            if [[ "$content" == "yes" ]]; then
                echo "[Step 5] Deinterleave..."
                deinterleave_start_time=$(date +%s)

                merged_bch_output=${bch_result_path}/Cumulative_BCH_decoded_bits_${s}.txt
                interleave_order_file=${configuration_files}/rs4000_250group_diagonal_interleave_order
                deinterleaved_output_dir=${output_dir}/codeword_for_rs_decoding
                mkdir -p "$deinterleaved_output_dir"

                flag_file_for_rs=${deinterleaved_output_dir}/RS_decoding_flag_${s}.txt
                deinterleaved_erasure_250groups=${deinterleaved_output_dir}/deinterleaved_erasure_record_${s}.txt
                
                ./deinterleave_rs_symbols \
                    "$merged_bch_output" \
                    "$interleave_order_file" \
                    "$deinterleaved_output_dir" \
                    "$flag_file_for_rs" \
                    "$deinterleaved_erasure_250groups"
                    
                deinterleave_end_time=$(date +%s)
                deinterleave_time=$((deinterleave_end_time - deinterleave_start_time))
                total_deinterleave_time=$((total_deinterleave_time + deinterleave_time))
                echo "[OK] Deinterleave finished (${deinterleave_time}s)"

                content_deinter=$(tr -d '[:space:]' < "$flag_file_for_rs")

                if [[ "$content_deinter" == "yes" ]]; then
                    echo "[Step 6] RS decoding & video recovery..."
                    rsdecoding_start_time=$(date +%s)

                    rs_decoding_for_recovery=${output_dir}/rs_decoding_for_recovery
                    mkdir -p "$rs_decoding_for_recovery"

                    thread=${THREAD_SETTING}
                    status_file=${rs_decoding_for_recovery}/status_file_${s}.txt
                    success_decoding_flag=${rs_decoding_for_recovery}/success_decoding_${s}.txt
                    recovery_video_path=${output_dir}/video_recovery
                    mkdir -p "$recovery_video_path"

                    bits_for_video_recovery=${recovery_video_path}/bits_for_video_recovery_${s}.txt
                    
                    ./rs_decode_parallel \
                        "$deinterleaved_output_dir" \
                        "$rs_decoding_for_recovery" \
                        "$thread" \
                        "$status_file" \
                        "$success_decoding_flag" \
                        "$bits_for_video_recovery"

                    rsdecoding_end_time=$(date +%s)
                    rsdecoding_time=$((rsdecoding_end_time - rsdecoding_start_time))
                    total_rsdecoding_time=$((total_rsdecoding_time + rsdecoding_time))
                    content_rs=$(tr -d '[:space:]' < "$success_decoding_flag")
                    if [[ "$content_rs" == "yes" ]]; then
                        echo "[OK] RS decoding successful (${rsdecoding_time}s)"

                        video_readout_start_time=$(date +%s)
                        inputStream=${bits_for_video_recovery}
                        scramble_code=${configuration_files}/watermark_448800bit.txt
                        video_recovery=${recovery_video_path}/recovered_video.mp4
                        
                        ./recover_video_from_bitstream "$scramble_code" "$inputStream"  "$video_recovery"
                        check_status $? "verification for file $s" "$error_log_file"

                        xdg-open "$video_recovery" 2>/dev/null
                        sleep 0.5s
                        wmctrl -r "$(basename "$video_recovery")" -e 0,1320,680,690,517

                        video_readout_end_time=$(date +%s)
                        video_readout_time=$((video_readout_end_time - video_readout_start_time))
                        total_video_readout_time=$((total_video_readout_time + video_readout_time))
                        echo "[OK] Video readout finished (${video_readout_time}s)"
                        echo -e "\033[32m[Group $s] RS decoding and recovery succeeded.\033[0m"
                        rs_success=1
                    else
                        echo "[WARN] Reed Solomon decoding failed for this group."
                    fi
                else
                    echo "[INFO] Deinterleaving completed, but RS erasure rate > 15%; skipping RS decoding."
                fi
            else
                echo "[INFO] Insufficient BCH results; RS decoding was not triggered."
            fi
        fi

        # ================== Summary ==================
        group_end_time=$(date +%s)
        group_total_time=$((group_end_time - group_start_time))
        cumulative_time=$((group_end_time - loop_start_time))
        
        echo "============================================================"
        printf "\033[32m[SUMMARY] Group %d \033[0m\n" \
            "$s" "$current_chip"
        printf "\033[32mNew reads: %d | Cum reads: %d | Cov: %s x\033[0m\n" \
            "$primer_reads_count" "$cumulative_reads" "${cumulative_coverage:-0.0000}"
        printf "\033[32mTime (this group): %ds | Total time: %ds\033[0m\n" \
            "$group_total_time" "$cumulative_time"
        echo "============================================================"

        echo "$s,${cumulative_coverage:-0.0000},$(date -d @$group_start_time +'%H:%M:%S'),$(date -d @$group_end_time +'%H:%M:%S'),$group_total_time,$cumulative_time,$(date +'%Y-%m-%d')" >> "$time_log_file"
        
        if [ "$rs_success" -eq 1 ]; then
            break
        fi

        s=$((s + 1))
    done

    loop_end_time=$(date +%s)
    loop_total_time=$((loop_end_time - loop_start_time))
    echo
    echo "============================================================"
    echo -e "\033[36m[Time Summary] Cumulative execution time per step:\033[0m"
    echo "------------------------------------------------------------"
    printf "  %-35s : %6ds\n" "Step 1 (Primer alignment)" "$total_primer_time"
    printf "  %-35s : %6ds\n" "Step 2 (Compare with refs)" "$total_compare_time"
    printf "  %-35s : %6ds\n" "Step 3 (FB + BCH decoding)" "$total_fb_time"
    printf "  %-35s : %6ds\n" "Step 4 (Updating BCH results)" "$total_bch_comb_time"
    printf "  %-35s : %6ds\n" "Step 5 (Deinterleave)" "$total_deinterleave_time"
    printf "  %-35s : %6ds\n" "Step 6a(RS decoding)" "$total_rsdecoding_time"
    printf "  %-35s : %6ds\n" "Step 6b(Video readout)" "$total_video_readout_time"
    echo "------------------------------------------------------------"
    echo -e "\033[32m[Done] All groups completed in ${loop_total_time}s ($(date '+%F %T'))\033[0m"
    echo "============================================================"
    
    echo "All,${cumulative_coverage:-0.0000},$loop_start_time,$loop_end_time,$loop_total_time,$loop_total_time,$(date +'%Y-%m-%d')" >> "$time_log_file"
}

main "$@"

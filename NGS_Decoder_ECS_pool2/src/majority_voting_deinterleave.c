#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>
#include <unistd.h>
#include <fcntl.h>

/* =====================================================================
 * Configuration Constants
 * ===================================================================== */
#define MAX_ADDR 1000000            // Maximum number of DNA strands
#define BIT_LENGTH 152              // Length of bit stream line
#define TARGET_START 66             // Address start bit index
#define TARGET_END 85               // Address end bit index
#define TARGET_LENGTH (TARGET_END - TARGET_START + 1)
#define RESULT_LENGTH 132           // Payload bits length
#define RS_CODE_NUM 11              // Number of RS symbols per strand
#define ROWS_PER_CODEWORD 4000      // Rows per interleaving matrix block
#define GROUP_SIZE 4000             // Same as ROWS_PER_CODEWORD
#define NUM_GROUPS 250              // Total groups (1M / 4000)
#define ERASURE_THRESHOLD 600       // Max allowed erasures per codeword
#define TOTAL_ELEMENTS (NUM_GROUPS * RS_CODE_NUM * ROWS_PER_CODEWORD)

/* =====================================================================
 * High-performance I/O Configuration
 * ===================================================================== */
#define OUT_BUF_SIZE (1024 * 1024)  // 1 MB output buffer
#define MAX_VAL 4096                // Max 12-bit symbol value + 1

/* =====================================================================
 * Global Variables & Lookup Tables
 * ===================================================================== */
// Lookup tables for fast integer-to-string conversion
char str_lut[MAX_VAL][8];
int  str_len[MAX_VAL];

const char *STR_ERASURE = "-1 ";
const int   LEN_ERASURE = 3;

// Data Matrices
int *raw_interleaved_data = NULL;   // Data before de-interleaving (Symbols)
int *deinterleaved_data = NULL;     // Data after de-interleaving
int *interleave_map = NULL;         // The mapping index

// Voting Structure
typedef struct {
    int vote_0[RESULT_LENGTH];
    int vote_1[RESULT_LENGTH];
    int total_votes;
} VoteResult;

/* =====================================================================
 * Helper Functions
 * ===================================================================== */

/**
 * Initialize lookup table for fast int->string conversion.
 */
void init_lookup_table() {
    for (int i = 0; i < MAX_VAL; i++) {
        sprintf(str_lut[i], "%d ", i);
        str_len[i] = strlen(str_lut[i]);
    }
}

/**
 * Recursive directory creation.
 */
void create_directory_recursive(const char *path) {
    char tmp[1024];
    char *p = NULL;
    size_t len;
    
    snprintf(tmp, sizeof(tmp), "%s", path);
    len = strlen(tmp);
    if (tmp[len - 1] == '/') tmp[len - 1] = 0;
    
    for (p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = 0;
            mkdir(tmp, 0775);
            *p = '/';
        }
    }
    mkdir(tmp, 0775);
}

/**
 * Binary string to Decimal converter.
 */
uint32_t bin_to_dec(const char *bin_str, int start, int end) {
    uint32_t result = 0;
    for (int i = start; i <= end; i++) {
        result = (result << 1) | (bin_str[i] - '0');
    }
    return result;
}

/**
 * Resets voting memory states.
 */
void reset_voting_state(VoteResult **found_votes, int *already_found, char **voted_bits) {
    for (uint32_t i = 0; i < MAX_ADDR; i++) {
        already_found[i] = 0;
        if (found_votes[i] != NULL) {
            memset(found_votes[i], 0, sizeof(VoteResult));
        }
        memset(voted_bits[i], 0, (RESULT_LENGTH + 1) * sizeof(char));
    }
}

/* =====================================================================
 * Main Program
 * ===================================================================== */
int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Usage: %s <input_bits_file> <interleave_order> <output_dir> <flag_file> <erasure_report> <summary_file>\n", argv[0]);
        printf("Example: %s BCH_decoded.txt interleave.bin ./output flag.txt erasures.txt summary.txt\n", argv[0]);
        return 1;
    }

    // --- Arguments ---
    char *input_file_path = argv[1];
    char *diagonal_interleave_order = argv[2];
    char *output_base_dir = argv[3];
    char *flag_file = argv[4];
    char *erasure_report_file = argv[5];
    char *summary_file_path = argv[6];

    // --- Initialization ---
    init_lookup_table();
    create_directory_recursive(output_base_dir);

    printf("--- DNA Storage: Voting, De-interleaving & Group Output ---\n");

    // --- Memory Allocation ---
    printf("Allocating memory...\n");
    
    // 1. Voting structures
    char **voted_bits = (char **)malloc(MAX_ADDR * sizeof(char *));
    VoteResult **found_votes = (VoteResult **)calloc(MAX_ADDR, sizeof(VoteResult *));
    int *already_found = (int *)calloc(MAX_ADDR, sizeof(int));
    
    for (uint32_t addr = 0; addr < MAX_ADDR; addr++) {
        voted_bits[addr] = (char *)calloc((RESULT_LENGTH + 1), sizeof(char));
    }

    // 2. Matrix buffers (Flattened 1D arrays for efficiency)
    raw_interleaved_data = (int *)malloc(TOTAL_ELEMENTS * sizeof(int));
    deinterleaved_data   = (int *)malloc(TOTAL_ELEMENTS * sizeof(int));
    interleave_map       = (int *)malloc(TOTAL_ELEMENTS * sizeof(int));

    if (!raw_interleaved_data || !deinterleaved_data || !interleave_map || !voted_bits) {
        perror("Memory allocation failed");
        return 1;
    }

    // Initialize raw data with Erasures (-1)
    for(int i=0; i<TOTAL_ELEMENTS; i++) raw_interleaved_data[i] = -1;

    // --- Step 1: Read Input & Perform Majority Voting ---
    printf("Step 1: Processing input file (Voting)...\n");
    FILE *file = fopen(input_file_path, "r");
    if (!file) { perror("Input file open error"); return 1; }

    char meta_line[256];
    char bit_stream[BIT_LENGTH + 16];
    char dummy_str[200];
    int ideal_index;
    double dummy_rate;

    while (1) {
        // Line 1: Metadata
        if (fgets(meta_line, sizeof(meta_line), file) == NULL) break;
        sscanf(meta_line, "%s %d %lf", dummy_str, &ideal_index, &dummy_rate);

        // Line 2: Bits
        if (fgets(bit_stream, sizeof(bit_stream), file) == NULL) break;

        // Parse Address
        uint32_t address = bin_to_dec(bit_stream, TARGET_START, TARGET_END);
        if (address >= MAX_ADDR) continue;

        if (found_votes[address] == NULL) {
            found_votes[address] = (VoteResult *)calloc(1, sizeof(VoteResult));
        }

        // Vote on Payload Bits
        int bit_idx = 0;
        for (int i = 0; i < BIT_LENGTH; i++) {
            if (i < TARGET_START || i > TARGET_END) {
                if (bit_stream[i] == '0') found_votes[address]->vote_0[bit_idx]++;
                else if (bit_stream[i] == '1') found_votes[address]->vote_1[bit_idx]++;
                bit_idx++;
            }
        }
        found_votes[address]->total_votes++;
        already_found[address] = 1;
    }
    fclose(file);

    // Calculate Dropout Stats
    int found_count = 0;
    for (int i = 0; i < MAX_ADDR; i++) { if (already_found[i]) found_count++; }
    double dropout_rate = (double)(MAX_ADDR - found_count) / MAX_ADDR * 100.0;
    printf(" -> Unique Addresses Found: %d\n", found_count);
    printf(" -> Dropout Rate: %.4f%%\n", dropout_rate);

    // --- Step 2: Determine Bits & Map to Symbols ---
    printf("Step 2: Mapping bits to RS symbols...\n");
    
    // Process every address (row)
    for (int addr = 0; addr < MAX_ADDR; addr++) {
        // 2a. Determine bits from votes
        if (already_found[addr] && found_votes[addr] != NULL) {
            for (int pos = 0; pos < RESULT_LENGTH; pos++) {
                if (found_votes[addr]->vote_1[pos] > found_votes[addr]->vote_0[pos]) voted_bits[addr][pos] = '1';
                else if (found_votes[addr]->vote_1[pos] < found_votes[addr]->vote_0[pos]) voted_bits[addr][pos] = '0';
                else voted_bits[addr][pos] = 'N'; // Ambiguous
            }
        } else {
            memset(voted_bits[addr], 'E', RESULT_LENGTH);
        }

        // 2b. Convert bits to RS Symbols (12 bits -> 1 int)
        for (int g2 = 0; g2 < RS_CODE_NUM; g2++) {
            int sym_val = 0;
            int base_idx = 12 * g2;
            int is_erasure = 0;

            for (int gg = 0; gg < 12; gg++) {
                char b = voted_bits[addr][base_idx + gg];
                if (b == 'E' || b == 'N') { is_erasure = 1; break; }
                if (b == '1') sym_val |= (1 << (11 - gg));
            }

            // Calculate linear index in the raw array
            long linear_idx = (long)addr * RS_CODE_NUM + g2;
            
            if (is_erasure) raw_interleaved_data[linear_idx] = -1;
            else raw_interleaved_data[linear_idx] = sym_val;
        }
    }

    // Free voting memory to save space
    for (int i=0; i<MAX_ADDR; i++) {
        if(found_votes[i]) free(found_votes[i]);
        free(voted_bits[i]);
    }
    free(found_votes);
    free(voted_bits);
    free(already_found);

    // --- Step 3: De-interleaving ---
    printf("Step 3: De-interleaving data...\n");

    FILE *f_order = fopen(diagonal_interleave_order, "rb");
    if (!f_order) { perror("Interleave order file error"); return 1; }
    size_t read_count = fread(interleave_map, sizeof(int), TOTAL_ELEMENTS, f_order);
    if (read_count != (size_t)TOTAL_ELEMENTS) {
        fprintf(stderr, "Map size mismatch! Expected %d, got %zu\n", TOTAL_ELEMENTS, read_count);
        return 1;
    }
    fclose(f_order);

    for (int i = 0; i < TOTAL_ELEMENTS; i++) {
        int target_idx = interleave_map[i];
        if (target_idx >= 0 && target_idx < TOTAL_ELEMENTS) {
            deinterleaved_data[target_idx] = raw_interleaved_data[i];
        } else {
            // If map points outside valid range (shouldn't happen with correct file)
             // Handle if necessary, or ignore
        }
    }

    // --- Step 4: Write Output Files (250 Groups) ---
    printf("Step 4: Writing 250 group files...\n");

    FILE *fp_report = fopen(erasure_report_file, "w");
    if (!fp_report) { perror("Report file create error"); return 1; }

    char *custom_buffer = (char*)malloc(OUT_BUF_SIZE);
    char out_path[512];
    int all_codewords_pass = 1;
    int failed_group_global = 0;

    for (int g = 0; g < NUM_GROUPS; g++) {
        snprintf(out_path, sizeof(out_path), "%s/group_%d.txt", output_base_dir, g);
        FILE *fp_out = fopen(out_path, "w");
        if (!fp_out) { perror("Output file create error"); exit(1); }

        int buf_offset = 0;

        // Iterate columns (RS codes) then rows
        for (int c = 0; c < RS_CODE_NUM; c++) {
            int erasure_count = 0;
            
            for (int r = 0; r < ROWS_PER_CODEWORD; r++) {
                // Calculate index in deinterleaved array
                // Structure: [Group 0 ... Group 249]
                // Within Group: [Col 0 (4000 rows), Col 1 (4000 rows)...]
                long index = (long)g * (RS_CODE_NUM * ROWS_PER_CODEWORD) 
                           + (long)c * ROWS_PER_CODEWORD + r;
                
                int val = deinterleaved_data[index];

                // Buffer write logic
                if (val == -1) {
                    erasure_count++;
                    if (buf_offset + LEN_ERASURE >= OUT_BUF_SIZE) {
                        fwrite(custom_buffer, 1, buf_offset, fp_out);
                        buf_offset = 0;
                    }
                    memcpy(custom_buffer + buf_offset, STR_ERASURE, LEN_ERASURE);
                    buf_offset += LEN_ERASURE;
                } else {
                    int len = str_len[val];
                    if (buf_offset + len >= OUT_BUF_SIZE) {
                        fwrite(custom_buffer, 1, buf_offset, fp_out);
                        buf_offset = 0;
                    }
                    memcpy(custom_buffer + buf_offset, str_lut[val], len);
                    buf_offset += len;
                }
            } // End row

            // Newline after each codeword (column)
            if (buf_offset + 1 >= OUT_BUF_SIZE) {
                fwrite(custom_buffer, 1, buf_offset, fp_out);
                buf_offset = 0;
            }
            custom_buffer[buf_offset++] = '\n';

            // Check erasure threshold
            fprintf(fp_report, "%d ", erasure_count);
            if (erasure_count >= ERASURE_THRESHOLD) {
                all_codewords_pass = 0;
                failed_group_global++;
            }
        } // End column

        if (buf_offset > 0) fwrite(custom_buffer, 1, buf_offset, fp_out);
        fclose(fp_out);
        fprintf(fp_report, "\n"); // Newline in report for next group
    }
    
    fclose(fp_report);
    if(custom_buffer) free(custom_buffer);

    // --- Step 5: Final Summaries ---
    printf("Step 5: Writing summaries...\n");

    // Pass/Fail Flag
    FILE *fp_flag = fopen(flag_file, "w");
    if (fp_flag) {
        fprintf(fp_flag, "%s\n", all_codewords_pass ? "yes" : "no");
        fclose(fp_flag);
    }

    // Summary Statistics File
    FILE *fw_sum = fopen(summary_file_path, "w");
    if (fw_sum) {
        double success_rate = (double)(NUM_GROUPS * RS_CODE_NUM - failed_group_global) 
                              / (NUM_GROUPS * RS_CODE_NUM) * 100.0;
        
        fprintf(fw_sum, "File:%s\tDropout:%.4f%%\tPassRate:%.4f%%\tFailedCW:%d\n", 
                input_file_path, dropout_rate, success_rate, failed_group_global);
        fclose(fw_sum);
    }

    // Cleanup
    free(raw_interleaved_data);
    free(deinterleaved_data);
    free(interleave_map);

    printf("Done. All outputs generated.\n");
    return 0;
}
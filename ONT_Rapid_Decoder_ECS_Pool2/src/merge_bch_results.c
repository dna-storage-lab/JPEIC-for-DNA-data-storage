#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>

#define MAX_ADDR 1000000 // 0-999,999
#define BUFFER_SIZE (256 * 1024) // 256KB buffer to reduce system calls

// Global flag array: used only to track which IDs have been recovered
uint8_t addr_flags[MAX_ADDR] = {0};

// Helper function: fast file append + bypass statistics (for BCH)
// Purpose: Append the contents of input_file to output_fp
// At the same time: parse the ">Address" lines and update global addr_flags
void append_file_and_count(const char *input_file, FILE *output_fp) 
{
    if (access(input_file, F_OK) == -1) return;

    FILE *fp_in = fopen(input_file, "r");
    if (!fp_in) {
        perror("Error opening input file");
        return;
    }

    // Set buffer
    char *io_buf = (char *)malloc(BUFFER_SIZE);
    if (io_buf) setvbuf(fp_in, io_buf, _IOFBF, BUFFER_SIZE);

    char line[1024]; 
    char last_char = '\n'; // Default assuming the previous character is a newline

    while (fgets(line, sizeof(line), fp_in)) 
    {
        // 1. Write to output file
        fputs(line, output_fp);

        // 2. Update statistics
        if (line[0] == '>') {
            int addr = atoi(line + 1);
            if (addr >= 0 && addr < MAX_ADDR) {
                addr_flags[addr] = 1;
            }
        }
        
        // 3. Record the last character of the current line
        size_t len = strlen(line);
        if (len > 0) {
            last_char = line[len - 1];
        }
    }

    // 4. Safety check: if the file does not end with a newline, force add one
    // This prevents file1's tail from merging with file2's head (e.g., >ID1ACTG>ID2...)
    if (last_char != '\n') {
        fputc('\n', output_fp);
    }
    fclose(fp_in);
    if (io_buf) free(io_buf);
}

// Helper function: pure file append (for DNA)
// Purpose: Append the contents of input_file to output_fp
void append_file_raw(const char *input_file, FILE *output_fp) 
{
    if (access(input_file, F_OK) == -1) return;

    FILE *fp_in = fopen(input_file, "rb"); // Binary read
    if (!fp_in) return;

    // Allocate large memory buffer for transport
    char *buffer = (char *)malloc(BUFFER_SIZE);
    if (!buffer) { fclose(fp_in); return; }

    size_t n;
    while ((n = fread(buffer, 1, BUFFER_SIZE, fp_in)) > 0) {
        fwrite(buffer, 1, n, output_fp);
    }
    fclose(fp_in);
    free(buffer); 
}

int main(int argc, char *argv[]) 
{
    // Parameter check [Changed: increased arg count to 9]
    if (argc != 9) {
        printf("Usage: %s <old_bch> <new_bch> <merged_bch> <old_dna> <new_dna> <merged_dna> <flag_file> <rate_file>\n", argv[0]);
        return 1;
    }

    char *old_bch_file = argv[1];
    char *new_bch_file = argv[2];
    char *merged_bch_output = argv[3];
    
    char *old_dna_file = argv[4];
    char *new_dna_file = argv[5];
    char *merged_dna_output = argv[6];
    
    char *flag_file = argv[7];
    char *rate_file = argv[8];
    // printf("Starting fast merge process...\n");

    // =========================================================
    // 1. Merge DNA files (pure copy, no logic)
    // =========================================================
    // Logic: clear merged_dna, then append old and new files in sequence
    FILE *f_dna_out = fopen(merged_dna_output, "wb"); // "wb" binary write, clear file
    if (f_dna_out) {
        // Append old file
        append_file_raw(old_dna_file, f_dna_out);
        // Append new file
        append_file_raw(new_dna_file, f_dna_out);
        
        fclose(f_dna_out);
        // printf("[DNA] Merged %s and %s -> %s\n", old_dna_file, new_dna_file, merged_dna_output);
    } else {
        perror("Failed to open merged DNA output");
    }

    // =========================================================
    // 2. Merge BCH files and calculate erasure rate
    // =========================================================
    FILE *f_bch_out = fopen(merged_bch_output, "w"); // "w" text write
    if (f_bch_out) {
        // Set output buffer to reduce disk write frequency
        char *out_buf = (char*)malloc(BUFFER_SIZE);
        if(out_buf) setvbuf(f_bch_out, out_buf, _IOFBF, BUFFER_SIZE);

        // A. Process old file (Copy & Count)
        append_file_and_count(old_bch_file, f_bch_out);

        // B. Process new file (Copy & Count)
        append_file_and_count(new_bch_file, f_bch_out);
        
        fclose(f_bch_out);
        if(out_buf) free(out_buf);
        
    } else {
        perror("Failed to open merged BCH output");
        return 1;
    }

    // =========================================================
    // 3. Calculate statistics
    // =========================================================
    int total_collected = 0;
    // Traverse the flag array and count how many positions are 1
    for (int i = 0; i < MAX_ADDR; i++) {
        if (addr_flags[i]) {
            total_collected++;
        }
    }

    int total_erasure = MAX_ADDR - total_collected;
    double erasure_rate = (double)total_erasure / MAX_ADDR;

    printf("  Total Unique Index Collected: %d\n", total_collected);
    printf("  Total Erasures:         %d\n", total_erasure);
    printf("  Erasure Rate:           %.4f%%\n", erasure_rate * 100.0);

    // =========================================================
    // 4. Output flag file (Yes/No)
    // =========================================================
    FILE *f_flag = fopen(flag_file, "w");
    if (f_flag) {
        // Decision criteria: Erasure rate < 15%
        if (erasure_rate < 0.15) {
            fprintf(f_flag, "yes\n");
            // printf("Decision: [YES] Decoding possible.\n");
        } else {
            fprintf(f_flag, "no\n");
            // printf("Decision: [NO] Need more reads.\n");
        }
        fclose(f_flag);
    } else {
        perror("Failed to write flag file");
    }
    // =========================================================
    // 5. Output erasure rate file
    // =========================================================
    FILE *f_rate = fopen(rate_file, "w");
    if (f_rate) {
        // Save as a decimal number (e.g., 0.123456)
        fprintf(f_rate, "%.6f\n", erasure_rate);
        // printf("Erasure rate saved to %s\n", rate_file);
        fclose(f_rate);
    } else {
        perror("Failed to write rate file");
    }
    return 0;
}

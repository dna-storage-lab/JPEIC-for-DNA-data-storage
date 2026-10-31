#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>

// =====================================================================
// Configuration
// =====================================================================
#define MAX_ADDR 1000000
#define RS_CODE_NUM 11
#define BITS_PER_SYMBOL 12
#define ROWS_PER_CODEWORD 4000
#define NUM_GROUPS 250
#define ERASURE_THRESHOLD 600

#define TOTAL_ELEMENTS (NUM_GROUPS * RS_CODE_NUM * ROWS_PER_CODEWORD)

// =====================================================================
// High-performance I/O configuration
// =====================================================================
#define OUT_BUF_SIZE (1024 * 1024) // 1 MB output buffer
#define MAX_VAL 4096               // maximum 12-bit symbol value + 1

// Lookup tables for fast integer-to-string conversion:
//   "0 " ... "4095 " and "-1 "
char str_lut[MAX_VAL][8];
int  str_len[MAX_VAL];

const char *STR_ERASURE = "-1 ";
const int   LEN_ERASURE = 3;

// =====================================================================
// Global buffers
// =====================================================================
int *raw_matrix = NULL;
int *deinterleaved_array = NULL;
int *interleave_map = NULL;

// =====================================================================
// Initialize integer-to-string lookup table
// =====================================================================
void init_lookup_table() {
    for (int i = 0; i < MAX_VAL; i++) {
        sprintf(str_lut[i], "%d ", i);
        str_len[i] = strlen(str_lut[i]);
    }
}

// =====================================================================
// Utility: create directory if it does not exist
// =====================================================================
void create_directory_if_not_exists(const char *dir) {
    struct stat st = {0};
    if (stat(dir, &st) == -1) {
        mkdir(dir, 0777);
    }
}

// =====================================================================
// Load BCH-decoded bitstream data (addressed format)
// =====================================================================
void load_bch_data(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (!fp) { perror("Open file"); exit(1); }

    // Large buffered I/O to accelerate sequential reads
    char *io_buf = (char*)malloc(4 * 1024 * 1024);
    if (io_buf) setvbuf(fp, io_buf, _IOFBF, 4 * 1024 * 1024);

    char line[1024];
    int addr = -1;
    int loaded_count = 0;

    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '>') {
            // Address header line
            addr = atoi(line + 1);
        } else {
            // Data line: decode 12-bit symbols
            if (addr >= 0 && addr < MAX_ADDR) {
                int *row_ptr = raw_matrix + addr * RS_CODE_NUM;
                const char *ptr = line;

                for (int i = 0; i < RS_CODE_NUM; i++) {
                    int sym_val = 0;
                    for (int bit = 0; bit < BITS_PER_SYMBOL; bit++) {
                        if (*ptr == '1') {
                            sym_val |= (1 << (BITS_PER_SYMBOL - 1 - bit));
                        }
                        ptr++;
                    }
                    row_ptr[i] = sym_val;
                }
                loaded_count++;
            }
            addr = -1;
        }
    }

    if (io_buf) free(io_buf);
    fclose(fp);

    printf("  Loaded %d sequences.\n", loaded_count);
}

// =====================================================================
// Main routine
// =====================================================================
int main(int argc, char *argv[]) {
    if (argc != 6) {
        printf("Usage: %s <bch_file> <interleave_order_file> <output_dir> <flag_file> <erasure_report_file>\n", argv[0]);
        return 1;
    }

    // --------------------------------------------------------------
    // Initialize lookup tables (critical for fast output)
    // --------------------------------------------------------------
    init_lookup_table();

    char *bch_file = argv[1];
    char *order_file = argv[2];
    char *output_dir = argv[3];
    char *flag_file = argv[4];
    char *erasure_report_file = argv[5];

    // --------------------------------------------------------------
    // Allocate global buffers
    // --------------------------------------------------------------
    raw_matrix = (int *)malloc(TOTAL_ELEMENTS * sizeof(int));
    deinterleaved_array = (int *)malloc(TOTAL_ELEMENTS * sizeof(int));
    interleave_map = (int *)malloc(TOTAL_ELEMENTS * sizeof(int));

    if (!raw_matrix || !deinterleaved_array || !interleave_map) {
        fprintf(stderr, "Memory allocation failed!\n");
        return 1;
    }

    // Initialize raw matrix with erasure markers (-1)
    memset(raw_matrix, 0xFF, TOTAL_ELEMENTS * sizeof(int));

    // --------------------------------------------------------------
    // Load interleaving map
    // --------------------------------------------------------------
    FILE *f_order = fopen(order_file, "rb");
    if (!f_order) { perror("Order file"); return 1; }

    size_t read_count = fread(interleave_map, sizeof(int), TOTAL_ELEMENTS, f_order);
    if (read_count != (size_t)TOTAL_ELEMENTS) {
        fprintf(stderr,
                "Error: Interleave map size mismatch! Expected %d ints, got %zu\n",
                TOTAL_ELEMENTS, read_count);
        fclose(f_order);
        return 1;
    }
    fclose(f_order);

    // --------------------------------------------------------------
    // Load BCH-decoded data
    // --------------------------------------------------------------
    load_bch_data(bch_file);

    // --------------------------------------------------------------
    // Deinterleaving
    // --------------------------------------------------------------
    for (int i = 0; i < TOTAL_ELEMENTS; i++) {
        int target_idx = interleave_map[i];
        if (target_idx >= 0 && target_idx < TOTAL_ELEMENTS) {
            deinterleaved_array[target_idx] = raw_matrix[i];
        }
    }

    // --------------------------------------------------------------
    // Split into groups and write output (optimized I/O path)
    // --------------------------------------------------------------
    create_directory_if_not_exists(output_dir);

    FILE *fp_report = fopen(erasure_report_file, "w");
    if (!fp_report) { perror("Report file"); exit(1); }

    int all_codewords_pass = 1;
    char out_path[512];

    char *custom_buffer = (char*)malloc(OUT_BUF_SIZE);
    int failed_group_global = 0;

    for (int g = 0; g < NUM_GROUPS; g++)
    {
        sprintf(out_path, "%s/group_%d.txt", output_dir, g);
        FILE *fp_out = fopen(out_path, "w");
        if (!fp_out) { perror("Output file"); exit(1); }

        int buf_offset = 0;

        for (int c = 0; c < RS_CODE_NUM; c++)
        {
            int erasure_count = 0;

            for (int r = 0; r < ROWS_PER_CODEWORD; r++)
            {
                long index = (long)g * (RS_CODE_NUM * ROWS_PER_CODEWORD)
                           + (long)c * ROWS_PER_CODEWORD + r;

                int val = deinterleaved_array[index];

                if (val == -1) {
                    erasure_count++;

                    if (buf_offset + LEN_ERASURE >= OUT_BUF_SIZE) {
                        fwrite(custom_buffer, 1, buf_offset, fp_out);
                        buf_offset = 0;
                    }
                    memcpy(custom_buffer + buf_offset, STR_ERASURE, LEN_ERASURE);
                    buf_offset += LEN_ERASURE;
                }
                else {
                    int len = str_len[val];
                    if (buf_offset + len >= OUT_BUF_SIZE) {
                        fwrite(custom_buffer, 1, buf_offset, fp_out);
                        buf_offset = 0;
                    }
                    memcpy(custom_buffer + buf_offset, str_lut[val], len);
                    buf_offset += len;
                }
            }

            // End of one RS codeword
            if (buf_offset + 1 >= OUT_BUF_SIZE) {
                fwrite(custom_buffer, 1, buf_offset, fp_out);
                buf_offset = 0;
            }
            custom_buffer[buf_offset++] = '\n';

            // Record erasure statistics
            fprintf(fp_report, "%d ", erasure_count);
            if (erasure_count >= ERASURE_THRESHOLD) {
                all_codewords_pass = 0;
                failed_group_global++;
            }
        }

        if (buf_offset > 0) {
            fwrite(custom_buffer, 1, buf_offset, fp_out);
        }
        fclose(fp_out);

        fprintf(fp_report, "\n");
    }

    printf("  Number of Reed-Solomon codewords with erasure < 15%%: %d (%.6f)\n",
           NUM_GROUPS * RS_CODE_NUM - failed_group_global,
           (double)(NUM_GROUPS * RS_CODE_NUM - failed_group_global)
           / (NUM_GROUPS * RS_CODE_NUM));

    fclose(fp_report);
    if (custom_buffer) free(custom_buffer);

    // --------------------------------------------------------------
    // Write global pass/fail flag
    // --------------------------------------------------------------
    FILE *fp_flag = fopen(flag_file, "w");
    if (fp_flag) {
        if (all_codewords_pass) fprintf(fp_flag, "yes\n");
        else                   fprintf(fp_flag, "no\n");
        fclose(fp_flag);
    }

    // --------------------------------------------------------------
    // Cleanup
    // --------------------------------------------------------------
    free(raw_matrix);
    free(deinterleaved_array);
    free(interleave_map);

    return 0;
}

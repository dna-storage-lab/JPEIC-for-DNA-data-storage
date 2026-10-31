#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// =====================================================================
// Configuration
// =====================================================================
#define DESIGNED_CODEGROUP 250
#define INFORM_LEN 448800           // bits per group
#define TARGET_TOTAL_BYTES 14024071 // target output size (bytes)

// Scrambling (watermark) sequence
char buff_scramble[INFORM_LEN];

// Input bitstream buffer
char bit_buffer[INFORM_LEN];

// Output buffer (~14 MB), allocated statically
static unsigned char output_bytes[TARGET_TOTAL_BYTES];

// XOR on ASCII bits ('0'/'1')
inline char dexor(char a, char b) {
    return (a == b) ? '0' : '1';
}

int main(int argc, char *argv[])
{
    // --------------------------------------------------------------
    // Argument parsing
    // --------------------------------------------------------------
    if (argc != 4) {
        printf("Usage: %s <watermark_file> <bitstream_file> <output_video_file>\n", argv[0]);
        return 1;
    }

    char *watermark_file = argv[1];
    char *bitstream_file = argv[2];
    char *output_file = argv[3];

    // --------------------------------------------------------------
    // Load scrambling code
    // --------------------------------------------------------------
    FILE *fp_scramble = fopen(watermark_file, "r");
    if (fp_scramble == NULL) {
        printf("Error: Cannot open watermark file: %s\n", watermark_file);
        return 1;
    }

    size_t scramble_len = fread(buff_scramble, sizeof(char), INFORM_LEN, fp_scramble);
    fclose(fp_scramble);

    if (scramble_len != INFORM_LEN) {
        printf("Warning: Watermark length mismatch! Expected %d, got %zu\n",
               INFORM_LEN, scramble_len);
    }

    // --------------------------------------------------------------
    // Open input bitstream
    // --------------------------------------------------------------
    FILE *fp_in = fopen(bitstream_file, "r");
    if (fp_in == NULL) {
        printf("Error: Cannot open bitstream file: %s\n", bitstream_file);
        return 1;
    }

    // --------------------------------------------------------------
    // Bitstream recovery (descramble + bit packing)
    // --------------------------------------------------------------
    long long global_bit_counter = 0;
    int code_group = 0;

    printf("Starting recovery process (Target: %d bytes)...\n", TARGET_TOTAL_BYTES);

    while (code_group < DESIGNED_CODEGROUP)
    {
        // Read one group
        size_t read_len = fread(bit_buffer, 1, INFORM_LEN, fp_in);
        if (read_len == 0) break;

        // Descramble and pack bits (MSB-first)
        for (int i = 0; i < read_len; i++)
        {
            char descrambled_bit = dexor(bit_buffer[i], buff_scramble[i]);

            if (descrambled_bit == '1')
            {
                long long byte_idx = global_bit_counter / 8;
                int bit_offset = global_bit_counter % 8;

                if (byte_idx < TARGET_TOTAL_BYTES) {
                    output_bytes[byte_idx] |= (1 << (7 - bit_offset));
                }
            }
            global_bit_counter++;
        }

        // Skip line delimiters between groups
        int ch = fgetc(fp_in);
        while (ch == '\n' || ch == '\r') {
            ch = fgetc(fp_in);
        }
        if (ch != EOF) ungetc(ch, fp_in);

        code_group++;
    }

    fclose(fp_in);

    // --------------------------------------------------------------
    // Write recovered output
    // --------------------------------------------------------------
    FILE *fp_out = fopen(output_file, "wb");
    if (fp_out == NULL) {
        printf("Error: Cannot create output video file: %s\n", output_file);
        return 1;
    }

    size_t written = fwrite(output_bytes, 1, TARGET_TOTAL_BYTES, fp_out);
    printf("Successfully wrote %zu bytes to %s\n", written, output_file);

    fclose(fp_out);
    return 0;
}

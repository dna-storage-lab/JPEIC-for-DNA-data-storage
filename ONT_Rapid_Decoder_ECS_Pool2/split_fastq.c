#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 4000000 

void print_usage(const char *prog_name) {
    printf("Usage: %s <input_fastq> <output_dir> <prefix> <reads_per_file> <num_files>\n", prog_name);
    printf("Example: %s input.fastq ./output subfastq 200000 900\n", prog_name);
}

int main(int argc, char *argv[]) {
    if (argc != 6) {
        print_usage(argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    const char *output_dir = argv[2];
    const char *prefix = argv[3];
    long reads_per_file = atol(argv[4]);
    int max_files = atoi(argv[5]);

    if (reads_per_file <= 0 || max_files <= 0) {
        fprintf(stderr, "Error: reads_per_file and num_files must be positive integers.\n");
        return 1;
    }

    FILE *f_in = fopen(input_path, "r");
    if (!f_in) {
        perror("Error opening input file");
        return 1;
    }

    char line[MAX_LINE_LENGTH];
    char output_path[1024];
    int file_index = 0;
    
    for (file_index = 0; file_index < max_files; file_index++) {
        size_t dir_len = strlen(output_dir);
        if (dir_len > 0 && (output_dir[dir_len - 1] == '/' || output_dir[dir_len - 1] == '\\')) {
            snprintf(output_path, sizeof(output_path), "%s%s_%d.fastq", output_dir, prefix, file_index);
        } else {
            snprintf(output_path, sizeof(output_path), "%s/%s_%d.fastq", output_dir, prefix, file_index);
        }

        printf("Writing to: %s ... ", output_path);
        FILE *f_out = fopen(output_path, "w");
        if (!f_out) {
            perror("Error opening output file");
            fclose(f_in);
            return 1;
        }

        long current_reads = 0;
        int end_of_input = 0;

        while (current_reads < reads_per_file) {
            for (int k = 0; k < 4; k++) {
                if (fgets(line, MAX_LINE_LENGTH, f_in) == NULL) {
                    end_of_input = 1;
                    break;
                }
                fputs(line, f_out);
            }

            if (end_of_input) {
                break;
            }
            current_reads++;
        }

        fclose(f_out);
        printf("Done (%ld reads)\n", current_reads);

        if (end_of_input) {
            printf("End of input file reached at file index %d.\n", file_index);
            break;
        }
    }

    fclose(f_in);
    printf("Splitting complete.\n");

    return 0;
}
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#define MAX_LINE_LENGTH (1 * 1024 * 1024) 

#define FILE_BUFFER_SIZE (1 * 1024 * 1024) 

int main(int argc, char *argv[]) 
{
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <original_sam> <reads_output> <log_output>\n", argv[0]);
        return 1;
    }

    const char *sam_filename = argv[1];
    const char *reads_filename = argv[2];
    const char *log_filename = argv[3];

    FILE *sam_file = fopen(sam_filename, "r");
    if (!sam_file) { perror("Error opening input file"); return 1; }

    FILE *reads_file = fopen(reads_filename, "w");
    if (!reads_file) { perror("Error opening output file"); fclose(sam_file); return 1; }
    
    FILE *log_file = fopen(log_filename, "w");
    if (!log_file) { perror("Error opening log file"); fclose(sam_file); fclose(reads_file); return 1; }

    char *in_buf = (char*)malloc(FILE_BUFFER_SIZE);
    char *out_buf = (char*)malloc(FILE_BUFFER_SIZE);
    char *sam_line = (char*)malloc(MAX_LINE_LENGTH);

    if (!in_buf || !out_buf || !sam_line) {
        fprintf(stderr, "Critical Error: Memory allocation failed! Your system might be OOM.\n");
        return 1;
    }

    setvbuf(sam_file, in_buf, _IOFBF, FILE_BUFFER_SIZE);
    setvbuf(reads_file, out_buf, _IOFBF, FILE_BUFFER_SIZE);

    int read_count = 0;
    int total_read_count = 0;
    int truncated_lines = 0;

    // Loop through the SAM file line by line
    while (fgets(sam_line, MAX_LINE_LENGTH, sam_file) != NULL) 
    {
        // 1. [Safety check] Handle empty or very short lines
        size_t len = strlen(sam_line);
        if (len < 3) continue; // Skip empty or extremely short lines to prevent pointer overflow

        // 2. [Safety check] Check for truncated lines
        if (sam_line[len - 1] != '\n' && !feof(sam_file)) {
            truncated_lines++;
            // Discard remaining characters to avoid contaminating the next line
            int ch;
            while ((ch = fgetc(sam_file)) != '\n' && ch != EOF);
            continue; 
        }

        if (sam_line[0] == '@') continue;

        total_read_count++;

        // --- Parsing logic (reinforced version) ---
        
        char *ptr = sam_line;
        char *qname_start = ptr;
        
        // Find the end of QNAME
        ptr = strchr(ptr, '\t');
        if (!ptr) continue; // No Tab, serious format error, skip this line
        
        // 
        // Ensure that ptr is valid when calculating length
        int qname_len = ptr - qname_start;
        if (qname_len <= 0) continue; 

        ptr++; // Skip Tab

        // Parse FLAG
        int flag = atoi(ptr);
        if (flag != 4) continue;

        // Jump to SEQ (Col 10)
        // We start from Col 2 (FLAG), need to skip 8 Tabs
        int cols_to_skip = 8;
        int format_error = 0;
        
        while (cols_to_skip > 0) {
            ptr = strchr(ptr, '\t'); 
            if (!ptr) {
                format_error = 1;
                break; // Missing columns, skip this line
            }
            ptr++;
            cols_to_skip--;
        }

        if (format_error || !*ptr || *ptr == '\n') continue; // No SEQ column found

        char *seq_start = ptr;
        
        // Find the end of SEQ (handle Tab, newline, or Windows carriage return)
        char *seq_end = strpbrk(seq_start, "\t\n\r");
        
        size_t seq_len;
        if (seq_end) {
            seq_len = seq_end - seq_start;
        } else {
            seq_len = strlen(seq_start);
        }

        if (seq_len == 0) continue; // Skip empty sequences

        // --- Output ---
        fputc('>', reads_file);
        fwrite(qname_start, 1, qname_len, reads_file);
        fputc('\n', reads_file);
        fwrite(seq_start, 1, seq_len, reads_file);
        fputc('\n', reads_file);
        
        read_count++;
    
    } 

    if (truncated_lines > 0) {
        fprintf(stderr, "\n[WARNING] %d lines were too long and skipped!\n", truncated_lines);
        fprintf(log_file, "WARNING: %d lines were too long and skipped.\n", truncated_lines);
    }

    double retention_rate = (total_read_count > 0) ? (double)read_count / total_read_count * 100.0 : 0.0;
    
    printf("  Total processed: %d, Extracted: %d, Retention: %.2f%%\n", total_read_count, read_count, retention_rate);

    fprintf(log_file, "Summary:\nTotal: %d\nExtracted: %d\nRetention: %.2f%%\n", 
            total_read_count, read_count, retention_rate);
    
    fclose(sam_file);
    fclose(reads_file);
    fclose(log_file);
    // Clean up
    if (in_buf) free(in_buf);
    if (out_buf) free(out_buf);
    if (sam_line) free(sam_line);

    return 0;
}

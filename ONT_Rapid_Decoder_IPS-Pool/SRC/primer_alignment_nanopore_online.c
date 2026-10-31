#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>
#include "edlib.h"

#define MAXLEN 400000
#define min_edit 2
int MAX_THREADS=32; // Number of threads
int payload_size=208; // Number of payloads
int valid_base_size=10;
void update_progress(long processed_count);
long count_total_reads(const char *filename);
// Add progress tracking variables
typedef struct {
    long total_reads;
    long processed_reads;
    pthread_mutex_t *progress_mutex;
    int show_progress;
} ProgressData;

ProgressData *global_progress = NULL;

typedef struct
{
    int thread_id;
    char **primers;
    int primer_num;
    int primer_len;
    int min_length;
    FILE **f_outputs;
    long *read_counts;
    int *match_counts;
    int *primer_edit_record;
    char *output_prefix;
    char *result_output_prefix;
    pthread_mutex_t *mutex;
} ThreadData;

void inv(char *x, int n)
{
    int i, j, m = (n - 1) / 2;
    char temp;
    for (i = 0; i <= m; i++) 
    {
        j = n - 1 - i;
        temp = x[i];
        x[i] = x[j];
        x[j] = temp;
    }
}

void complement(char *input, char *output, int len) {
    int i;
    for (i = 0; i < len; i++) {
        if (input[i] == 'A') {
            output[i] = 'T';
        } else if (input[i] == 'T') {
            output[i] = 'A';
        } else if (input[i] == 'G') {
            output[i] = 'C';
        } else if (input[i] == 'C') {
            output[i] = 'G';
        } else {
            output[i] = 'N';
        }
    }
}

int fileSplitting(const char * filename,const char *output_prefix) // File splitting
{
    int a = 0; // Base number of groups per thread (4 lines per group)
    int b = 0; // Number of first b threads that need additional groups
    int line_groups = 0; // Total groups in file (4 lines per group)
    char line1[MAXLEN] = {0};
    char line2[MAXLEN] = {0};
    char line3[MAXLEN] = {0};
    char line4[MAXLEN] = {0};
    char c;

    // Step 1: Calculate total groups in file (4 lines per group)
    FILE *testfile = fopen(filename, "rt");
    if (testfile == NULL) {
        printf("Can't find testfile %s \n", filename);
        return 0;
    }

    // Count groups (4 lines per group)
    while (fgets(line1, MAXLEN, testfile) != NULL) 
    {
        if (fgets(line2, MAXLEN, testfile) != NULL &&
            fgets(line3, MAXLEN, testfile) != NULL &&
            fgets(line4, MAXLEN, testfile) != NULL) 
            {
            line_groups++;
        }
    }
    fclose(testfile);

    // printf("Total number of line groups (4 lines each): %d\n", line_groups);
    // Adjust thread count if file groups are less than threads
    if (line_groups < MAX_THREADS) {
        printf("Adjusting number of threads from %d to %d to match number of line groups\n", MAX_THREADS, line_groups);
        MAX_THREADS = line_groups;
    }
    // Calculate groups per thread a and additional groups for first b threads
    a = line_groups / MAX_THREADS;
    b = line_groups % MAX_THREADS;

    // printf("Threads: %d, Each thread processes at least %d groups.\n", MAX_THREADS, a);
    // printf("First %d threads will process %d groups each.\n", b, a + 1);
    // printf("Remaining %d threads will process %d groups each.\n", MAX_THREADS - b, a);

    // Step 2: Open all sub files
    FILE **subtestfile = (FILE **)malloc(MAX_THREADS * sizeof(FILE *));
    if (subtestfile == NULL) {
        printf("Memory allocation failed for subtestfile.\n");
        return 0;
    }

    char subtestfilename[200] = {0};
    for (int i = 0; i < MAX_THREADS; i++) 
    {
        sprintf(subtestfilename, "%s/Temporary/sub_%d.fastq", output_prefix,i);
        subtestfile[i] = fopen(subtestfilename, "wt");
        if (subtestfile[i] == NULL) {
            printf("Failed to open subtestfile %s\n", subtestfilename);
            // Close opened files and free memory
            for (int j = 0; j < i; j++) {
                fclose(subtestfile[j]);
            }
            free(subtestfile);
            return 0;
        }
    }

    // Step 3: Reopen main file for allocation
    testfile = fopen(filename, "rt");
    if (testfile == NULL) {
        printf("Can't find testfile %s \n", filename);
        // Close all sub files and free memory
        for (int i = 0; i < MAX_THREADS; i++) {
            fclose(subtestfile[i]);
        }
        free(subtestfile);
        return 0;
    }

    int current_thread = 0; // Current thread index for allocation
    int groups_assigned = 0; // Groups assigned to current thread
    int groups_to_assign = (b > 0) ? (a + 1) : a; // Groups to assign to current thread

    // Read 4 lines (one group) each time
    while (fgets(line1, MAXLEN, testfile) != NULL) 
    {
        if (fgets(line2, MAXLEN, testfile) != NULL &&
            fgets(line3, MAXLEN, testfile) != NULL &&
            fgets(line4, MAXLEN, testfile) != NULL) {
            // Write current group to current thread's sub file
            fprintf(subtestfile[current_thread], "%s", line1);
            fprintf(subtestfile[current_thread], "%s", line2);
            fprintf(subtestfile[current_thread], "%s", line3);
            fprintf(subtestfile[current_thread], "%s", line4);
            groups_assigned++;
            
            // Check if current thread has been assigned required groups
            if (groups_assigned >= groups_to_assign) 
            {
                current_thread++;
                groups_assigned = 0;
                if (current_thread < b) {
                    groups_to_assign = a + 1;
                } else if (current_thread < MAX_THREADS) 
                {
                    groups_to_assign = a;
                }
            }
        }
    }

    fclose(testfile);

    // Close all sub files
    for (int i = 0; i < MAX_THREADS; i++) {
        fclose(subtestfile[i]);
    }

    free(subtestfile);

    // printf("Data has been successfully split into %d subfiles.\n", MAX_THREADS);
	return a;
}

void recursive_primer_matching_forward(ThreadData *data, char *line, int seqlen, int read_num, long *valid_base_count, long *high_error_base_count, int *match_count, int primer_index)
{
    char primerL[data->primer_len + 1];
    char primerR[data->primer_len + 1];
    EdlibAlignResult result1;
    EdlibAlignResult result2;

    memcpy(primerL, data->primers[primer_index * 2], data->primer_len);
    memcpy(primerR, data->primers[primer_index * 2 + 1], data->primer_len);
    primerL[data->primer_len] = '\0';
    primerR[data->primer_len] = '\0';

    // Align primers (forward)
    result1 = edlibAlign(primerL, data->primer_len, line, seqlen, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    result2 = edlibAlign(primerR, data->primer_len, line, seqlen, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

    int start = *result1.startLocations;
    int end_forword = *result1.endLocations;
    int end = *result2.endLocations;
    int start_backword = *result2.startLocations;

    if (result1.editDistance <= min_edit && result2.editDistance <= min_edit)
    {
        
        if(start_backword>end_forword && (start_backword-end_forword-1)>=payload_size-valid_base_size && (start_backword-end_forword-1)<=payload_size+valid_base_size)
        {
            (*match_count)++;
            (*valid_base_count)+=(start_backword-end_forword-1);
            pthread_mutex_lock(data->mutex);
            fprintf(data->f_outputs[primer_index], ">%d_Read_%d_Match_%d_forward\n", data->thread_id, read_num, *match_count);
            for (int j = end_forword+1; j < start_backword; j++)
            {
                fputc(line[j], data->f_outputs[primer_index]);
            }
            fputc('\n', data->f_outputs[primer_index]);
            fflush(data->f_outputs[primer_index]);
            data->primer_edit_record[primer_index]++; 
            pthread_mutex_unlock(data->mutex);
        }
        // If middle region is long, remove both primers and iterate again
        if((start_backword-end_forword-1)>payload_size+(2*valid_base_size))
        {
            char middle_subseq[start_backword - end_forword];
            strncpy(middle_subseq, line + end_forword + 1, start_backword - end_forword - 1);
            middle_subseq[start_backword - end_forword - 1] = '\0';
            
            recursive_primer_matching_forward(data, middle_subseq, start_backword - end_forword - 1, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
        }
        if(start_backword>end_forword)
        {
            // Middle region has been read, continue recursion to process left and right remaining sequences
            if(start!=0)
            {
                char left_subseq[start + 1];
                strncpy(left_subseq, line, start);
                left_subseq[start] = '\0';
                recursive_primer_matching_forward(data, left_subseq, start, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
            }
            if(end!=seqlen-1)
            {
                char right_subseq[seqlen - end];
                strncpy(right_subseq, line + end + 1, seqlen - end-1);
                right_subseq[seqlen - end-1] = '\0';
                recursive_primer_matching_forward(data, right_subseq, seqlen - end - 1, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
            }
        }
        else
        {
            // Continue recursion to process left and right remaining sequences
            if(start_backword!=0)
            {
                char left_subseq[end + 2];
                strncpy(left_subseq, line, end + 1);
                left_subseq[end + 1] = '\0';
                recursive_primer_matching_forward(data, left_subseq, end + 1, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
            }
            if(end_forword!=seqlen-1)
            {
                char right_subseq[seqlen - start+1];
                strncpy(right_subseq, line + start, seqlen - start);
                right_subseq[seqlen - start] = '\0';
                recursive_primer_matching_forward(data, right_subseq, seqlen - start, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
            }
        }
    }

    edlibFreeAlignResult(result1);
    edlibFreeAlignResult(result2);
}

void recursive_primer_matching_reverse(ThreadData *data, char *line, int seqlen, int read_num, long *valid_base_count, long *high_error_base_count, int *match_count, int primer_index)
{
    char primerL[data->primer_len + 1];
    char primerR[data->primer_len + 1];
    char primerL_rv[data->primer_len + 1];
    char primerR_rv[data->primer_len + 1];
    char line_temp[MAXLEN]={0};
    EdlibAlignResult result3;
    EdlibAlignResult result4;

    memcpy(primerL, data->primers[primer_index * 2], data->primer_len);
    memcpy(primerR, data->primers[primer_index * 2 + 1], data->primer_len);
    primerL[data->primer_len] = '\0';
    primerR[data->primer_len] = '\0';

    // Handle reverse complement of primers
    inv(primerL, data->primer_len);
    complement(primerL, primerL_rv, data->primer_len);
    inv(primerR, data->primer_len);
    complement(primerR, primerR_rv, data->primer_len);
    primerL_rv[data->primer_len] = '\0';
    primerR_rv[data->primer_len] = '\0';

    // Align reverse complement primers
    result4 = edlibAlign(primerL_rv, data->primer_len, line, seqlen, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    result3 = edlibAlign(primerR_rv, data->primer_len, line, seqlen, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    
    int start = *result3.startLocations;
    int end_forword = *result3.endLocations;
    int start_backword = *result4.startLocations;
    int end = *result4.endLocations;

    if (result3.editDistance <= min_edit && result4.editDistance <= min_edit) 
    {
        if(start_backword>end_forword && (start_backword-end_forword-1)>=payload_size-valid_base_size && (start_backword-end_forword-1)<=payload_size+valid_base_size)
        {
            (*match_count)++;
            (*valid_base_count)+=(start_backword-end_forword-1);
            for (int j = end_forword+1; j < start_backword; j++)
            {
                line_temp[j - end_forword-1] = line[j];
            }
            inv(line_temp, start_backword-end_forword-1);
            complement(line_temp, line_temp, start_backword-end_forword-1);
            
            pthread_mutex_lock(data->mutex);
            fprintf(data->f_outputs[primer_index], ">%d_Read_%d_Match_%d_backward\n", data->thread_id, read_num, *match_count);
               
            for (int j = 0; j < start_backword-end_forword-1; j++) 
            {
                fputc(line_temp[j], data->f_outputs[primer_index]);
            }

            fputc('\n', data->f_outputs[primer_index]);
            fflush(data->f_outputs[primer_index]);
            data->primer_edit_record[primer_index]++;   
            pthread_mutex_unlock(data->mutex);
        }
        // If middle region is long, remove both primers and iterate again
        if((start_backword-end_forword-1)>payload_size+(2*valid_base_size))
        {
            char middle_subseq[start_backword - end_forword];
            strncpy(middle_subseq, line + end_forword + 1, start_backword - end_forword - 1);
            middle_subseq[start_backword - end_forword - 1] = '\0';
            recursive_primer_matching_reverse(data, middle_subseq, start_backword - end_forword - 1, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
        }

        if(start_backword>end_forword)
        {
           // Middle region has been read, continue recursion to process left and right remaining sequences
            if(start!=0)
            {
                char left_subseq[start + 1];
                strncpy(left_subseq, line, start);
                left_subseq[start] = '\0';
                recursive_primer_matching_reverse(data, left_subseq, start, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
            }
            if(end!=seqlen-1)
            {
                char right_subseq[seqlen - end];
                strncpy(right_subseq, line + end + 1, seqlen - end-1);
                right_subseq[seqlen - end-1] = '\0';
                recursive_primer_matching_reverse(data, right_subseq, seqlen - end - 1, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
            }
        }
        else
        {
            // Continue recursion to process left and right remaining sequences
            if(start_backword!=0)
            {
                char left_subseq[end + 2];
                strncpy(left_subseq, line, end + 1);
                left_subseq[end + 1] = '\0';
                recursive_primer_matching_reverse(data, left_subseq, end + 1, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
            }
            if(end_forword!=seqlen-1)
            {
                char right_subseq[seqlen - start+1];
                strncpy(right_subseq, line + start, seqlen - start);
                right_subseq[seqlen - start] = '\0';
                recursive_primer_matching_reverse(data, right_subseq, seqlen - start, read_num, valid_base_count, high_error_base_count, match_count, primer_index);
            }
        } 
    }

    edlibFreeAlignResult(result3);
    edlibFreeAlignResult(result4);
}

void* process_reads(void *arg) 
{
    ThreadData *data = (ThreadData*) arg;
    char line[MAXLEN];
    int line_number = 0;
    int read_num = 0;
    int progress_update_interval = 500; // Update progress every 500 reads

    // Open FASTQ file
    char subfastq_file[100] = {0};
    sprintf(subfastq_file, "%s/Temporary/sub_%d.fastq", data->output_prefix, data->thread_id);
    FILE *f_fastq = fopen(subfastq_file, "r");
    if (!f_fastq) {
        fprintf(stderr, "Error opening FASTQ file: %s\n", subfastq_file);
        exit(1);
    }
    long sequencing_base_count = 0;
    long valid_base_count = 0;
    long high_error_base_count = 0;
    while (fgets(line, MAXLEN, f_fastq) != NULL) 
    {
        
        line_number++;
        if (line_number % 4 == 2) 
        { // Process sequence line
            read_num++;
            
            int seqlen = strlen(line) - 1; // Remove newline character
            line[seqlen] = '\0';
            sequencing_base_count+=seqlen; // Count total sequencing fragment length
            data->read_counts[0]=sequencing_base_count;
            if (seqlen < data->min_length-10) // Filter out too short fragments <160
            {
                data->match_counts[0]++;
                continue;
            }
            

            // Recursively align each primer pair
            int match_count = 0;
            valid_base_count=0;
            high_error_base_count=0;
            for (int i = 0; i < data->primer_num / 2; i++)
            {
                // Call forward and reverse alignment functions respectively
                recursive_primer_matching_forward(data, line, seqlen, read_num, &valid_base_count, &high_error_base_count, &match_count, i);
                recursive_primer_matching_reverse(data, line, seqlen, read_num, &valid_base_count, &high_error_base_count, &match_count, i);
            }
            data->read_counts[1]+=valid_base_count;
            data->read_counts[2]+=high_error_base_count;
            if (match_count >= 0 && match_count < 100) 
            {
                data->match_counts[match_count]++;
            }  
            else {
                data->match_counts[99]++;
            }
        }
    }
    
    // Free resources
    fclose(f_fastq);
    return NULL;
}

// Function to count total reads
long count_total_reads(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file to count reads: %s\n", filename);
        return 0;
    }
    
    long line_count = 0;
    char buffer[MAXLEN];
    
    while (fgets(buffer, MAXLEN, file) != NULL) {
        line_count++;
    }
    
    fclose(file);
    return line_count / 4; // FASTQ file has 4 lines per read
}

void update_progress(long processed_count) 
{
    if (global_progress && global_progress->show_progress) 
    {
        pthread_mutex_lock(global_progress->progress_mutex);
        global_progress->processed_reads += processed_count;
        
        double percentage = (double)global_progress->processed_reads / global_progress->total_reads * 100.0;
        if (percentage > 100.0) percentage = 100.0;
        
        // Draw progress bar (single line update)
        int bar_width = 50;
        int filled = (int)(percentage / 100.0 * bar_width);
        
        printf("\r");  // Return to the beginning of the line
        printf("Progress: [");
        for (int i = 0; i < bar_width; i++) {
            if (i < filled) {
            printf("=");
            } else if (i == filled) {
            printf(">");
            } else {
            printf(" ");
            }
        }
        printf("] %.1f%% (%ld/%ld reads)", percentage, global_progress->processed_reads, global_progress->total_reads);
        fflush(stdout);
        
        if (global_progress->processed_reads >= global_progress->total_reads) {
            printf("\n");  // New line only when complete
        }
        
        pthread_mutex_unlock(global_progress->progress_mutex);
    }
}

int main(int argc, char *argv[])
{
    
    if (argc != 14)
    {
        fprintf(stderr, "Usage: %s <fastq_file> <primer_list> <primer_num> <primer_len> <output_file> <min_length> <match_count_file>\n", argv[0]);
        return 1;
    }

    char *fastq_file = argv[1];
    char *primer_list = argv[2];
    int primer_num = atoi(argv[3]);
    int primer_len = atoi(argv[4]);
    char *output_prefix = argv[5];
    int min_length = atoi(argv[6]);
    char *reads_primer_match_counts = argv[7];
    char *result_output_prefix = argv[8];
    MAX_THREADS = atoi(argv[9]);
    char *valid_base_percentage_file = argv[10];
    payload_size=atoi(argv[11]);
    valid_base_size=atoi(argv[12]);
    int file_index_number=atoi(argv[13]);

    // Initialize progress tracking
    ProgressData progress_data;
    pthread_mutex_t progress_mutex = PTHREAD_MUTEX_INITIALIZER;
    
    progress_data.total_reads = count_total_reads(fastq_file);
    progress_data.processed_reads = 0;
    progress_data.progress_mutex = &progress_mutex;
    progress_data.show_progress = 1;
    global_progress = &progress_data;
    
    printf("Number of reads to process: %ld\n", progress_data.total_reads);

    // Read primer list
    char **primers = (char **)malloc(primer_num * sizeof(char *));
    for (int i = 0; i < primer_num; i++) {
        primers[i] = (char *)malloc((primer_len + 1) * sizeof(char));
    }

    FILE *f_primer_list = fopen(primer_list, "r");
    if (!f_primer_list) {
        fprintf(stderr, "Error opening primer list file: %s\n", primer_list);
        return 1;
    }

    for (int i = 0; i < primer_num; i++) 
    {
        if (fscanf(f_primer_list, "%s", primers[i]) != 1) 
        {
            fprintf(stderr, "Error reading primer %d\n", i);
            return 1;
        }
    }
    fclose(f_primer_list);

    FILE *f_reads_primer_match_counts = fopen(reads_primer_match_counts, "wt");
    if (!f_reads_primer_match_counts) {
        fprintf(stderr, "Error opening match_counts file: %s\n", reads_primer_match_counts);
        return 1;
    }
    FILE *f_valid_base_percentage = fopen(valid_base_percentage_file, "wt");
    if (!f_valid_base_percentage) {
            fprintf(stderr, "Error opening valid_base_percentage file: %s\n", valid_base_percentage_file);
            return 1;
    }

    // Create output files
    FILE **f_outputs = (FILE **)malloc((primer_num / 2) * sizeof(FILE *));
    char output_filename[256];
    for (int i = 0; i < primer_num / 2; i++)
    {
        snprintf(output_filename, sizeof(output_filename), "%s/primer_pair_%d_%d.fasta", result_output_prefix, i + 1,file_index_number);
        f_outputs[i] = fopen(output_filename, "w");
        if (!f_outputs[i]) 
        {
            fprintf(stderr, "Error opening output file: %s\n", output_filename);
            return 1;
        }
    }

    fileSplitting(fastq_file,output_prefix);
    
    // Create and initialize thread data
    pthread_t threads[MAX_THREADS];
    ThreadData thread_data[MAX_THREADS];
    
    int **match_counts = (int **)malloc(MAX_THREADS * sizeof(int *));
    for (int i = 0; i < MAX_THREADS; i++)
    {
        match_counts[i] = (int *)calloc(100, sizeof(int)); // Preset maximum reads per thread, e.g., 100000
    }
    // Initialize match_counts array for each thread
    for (int i = 0; i < MAX_THREADS; i++)
    {
        for (int j = 0; j < 100; j++) 
        {
            match_counts[i][j] = 0;
        }
    }

    long **read_counts = (long **)malloc(MAX_THREADS * sizeof(long *));
    for (int i = 0; i < MAX_THREADS; i++) 
    {
        read_counts[i] = (long *)calloc(3, sizeof(long)); // Preset maximum reads per thread, e.g., 100000
    }
    // Initialize read_counts array for each thread
    for (int i = 0; i < MAX_THREADS; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            read_counts[i][j] = 0;
        }
    }  

    int primer_edit_record[10] = {0};
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    for (int i = 0; i < MAX_THREADS; i++) 
	{
        thread_data[i].thread_id = i;
        thread_data[i].primers = primers;
        thread_data[i].primer_num = primer_num;
        thread_data[i].primer_len = primer_len;
        thread_data[i].min_length = min_length;
        thread_data[i].f_outputs = f_outputs;  // valid reads
        thread_data[i].match_counts = match_counts[i]; // Each thread only sees its own array
        thread_data[i].read_counts = read_counts[i]; // Each thread only sees its own array
        thread_data[i].primer_edit_record = primer_edit_record;
        thread_data[i].output_prefix = output_prefix;
        thread_data[i].result_output_prefix = result_output_prefix;
        thread_data[i].mutex = &mutex;

        pthread_create(&threads[i], NULL, process_reads, (void*)&thread_data[i]);
    }

    // Join threads
    for (int i = 0; i < MAX_THREADS; i++) 
	{
        pthread_join(threads[i], NULL);
    }

    for (int i = 0; i < primer_num / 2; i++)
	{
        printf("The number of reads identified by primer %d is %d\n", i+1, primer_edit_record[i]);
    }
    int summary_match_counts[100] = {0};
    long summary_base_counts[3] = {0};
    for (int i = 0; i < MAX_THREADS; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            summary_base_counts[j] += read_counts[i][j];
        }
    }
    printf("Total_base_counts\tvalid_counts\tinvalid_counts\n");
    printf("%ld\t%.6f\t%.6f\n", summary_base_counts[0], (double)summary_base_counts[1]/summary_base_counts[0], 1-((double)summary_base_counts[1]/summary_base_counts[0]));

    fprintf(f_valid_base_percentage, "Total base counts: %ld\n", summary_base_counts[0]);
    fprintf(f_valid_base_percentage, "Valid base counts: %lf\n", (double)summary_base_counts[1]/summary_base_counts[0]);
    for (int i = 0; i < primer_num / 2; i++)
    {
        fprintf(f_valid_base_percentage, "primer pair %d matched reads: %d\n", i + 1, primer_edit_record[i]);
    }
    for (int i = 0; i < MAX_THREADS; i++) 
    {
        for (int j = 0; j < 100; j++) 
        {
            summary_match_counts[j] += match_counts[i][j];
        }
    }
    for (int i = 0; i < 100; i++)
	{
        fprintf(f_reads_primer_match_counts, "%d\n", summary_match_counts[i]);
    }

    for (int i = 0; i < primer_num / 2; i++) 
	{
        fclose(f_outputs[i]);
    }
    free(f_outputs);

    for (int i = 0; i < primer_num; i++) 
	{
        free(primers[i]);
    }
    free(primers);
    fclose(f_reads_primer_match_counts);

    for (int i = 0; i < MAX_THREADS; i++) 
    {
        free(match_counts[i]);
    }
    free(match_counts);
    fclose(f_valid_base_percentage);
    return 0;
}

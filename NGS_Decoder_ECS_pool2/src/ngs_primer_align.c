#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>
#include "edlib.h"

// ==== Configuration Parameters ====
#define MAXLEN 4000000
#define MIN_EDIT 2
#define PRIMER_LEN 10
#define BATCH_SIZE 100  // Batch size

// ==== Global Variables ====
int min_length = 130; 
int oligo_len = 170;
int payload_size = 130;
int valid_bias_size = 5; //valid bias range

// ==== Data Structures ====

// Task batch structure
typedef struct TaskBatch {
    char **sequences;       // Store multiple lines of data
    int *read_ids;          // Store read_num for each line
    int count;              // Actual count in the current batch
    struct TaskBatch *next;
} TaskBatch;

// Thread pool queue
typedef struct {
    TaskBatch *head;
    TaskBatch *tail;
    int count;          // Number of batches currently in the queue
    int max_size;       // Maximum queue capacity (e.g., 1000 batches)
    pthread_mutex_t mutex;
    pthread_cond_t cond_empty; // Consumer waits on this condition
    pthread_cond_t cond_full;  // Producer waits on this condition
    int shutdown;
} TaskQueue;

// Initialization
void queue_init(TaskQueue *q, int max_size) {
    q->head = q->tail = NULL;
    q->count = 0;
    q->max_size = max_size;
    pthread_mutex_init(&q->mutex, NULL);
    pthread_cond_init(&q->cond_empty, NULL);
    pthread_cond_init(&q->cond_full, NULL);
    q->shutdown = 0;
}

// Thread-private data
typedef struct {
    int thread_id;
    TaskQueue *queue;
    
    char p_fwd_L[PRIMER_LEN + 1];
    char p_fwd_R[PRIMER_LEN + 1];
    char p_rev_L[PRIMER_LEN + 1];
    char p_rev_R[PRIMER_LEN + 1];
    
    int processed_count;
    int valid_match_count;
    int *local_match_counts;
    
    FILE *f_out_result;
    FILE *f_out_len;
    char *temp_buf; 
} ThreadWorkerData;

// ==== Helper Functions (time) ====

double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time, NULL)) return 0;
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

void inv(char *x, int n) {
    int i, j;
    char temp;
    for (i = 0; i < n / 2; i++) {
        j = n - 1 - i;
        temp = x[i];
        x[i] = x[j];
        x[j] = temp;
    }
}

void complement(char *input, char *output, int len) {
    for (int i = 0; i < len; i++) {
        switch(input[i]) {
            case 'A': output[i] = 'T'; break;
            case 'T': output[i] = 'A'; break;
            case 'G': output[i] = 'C'; break;
            case 'C': output[i] = 'G'; break;
            default:  output[i] = 'N'; break;
        }
    }
}

// Check whether the length meets the requirements
inline int is_valid_length(int len) {
    return abs(len - payload_size) <= valid_bias_size;
}
// Define a structure to temporarily store match results for easy comparison
typedef struct {
    int valid;      // 0:Invalid, 1:Valid
    int is_reverse; // 0:Forward, 1:Reverse
    int is_dual;    // 0:Single, 1:Dual
    int start_idx;  // Start position in the current seq (inclusive)
    int end_idx;    // End position in the current seq (inclusive)
    int payload_start; // Start of the payload segment
    int payload_len;   // Length of the payload
    const char *type_str; // For logging: "forward_dual", etc.
} MatchResult;


// Helper: Check match quality and return result (Final Version)
// Helper: Check match quality and return result (Final Version)
MatchResult check_alignment(ThreadWorkerData *data, const char *line, int seqlen, int is_reverse)
{
    MatchResult res = {0};
    res.is_reverse = is_reverse;

    // 1. Get primer sequences
    const char *primerL = is_reverse ? data->p_rev_L : data->p_fwd_L;
    const char *primerR = is_reverse ? data->p_rev_R : data->p_fwd_R;

    // 2. Edlib global alignment(HW)
    static const EdlibAlignConfig config = { .k = -1, .mode = EDLIB_MODE_HW, .task = EDLIB_TASK_PATH };
    EdlibAlignResult resL = edlibAlign(primerL, PRIMER_LEN, line, seqlen, config);
    EdlibAlignResult resR = edlibAlign(primerR, PRIMER_LEN, line, seqlen, config);

    int found_L = (resL.numLocations > 0 && resL.editDistance != -1 && resL.editDistance <= MIN_EDIT);
    int found_R = (resR.numLocations > 0 && resR.editDistance != -1 && resR.editDistance <= MIN_EDIT);

    // Default to the first location
    int start_L = found_L ? *resL.startLocations : -1;
    int end_L   = found_L ? *resL.endLocations : -1;
    int start_R = found_R ? *resR.startLocations : -1;
    int end_R   = found_R ? *resR.endLocations : -1;

    int solved_perfectly = 0; // Mark whether a valid solution was found in the first step

    // =========================================================
    // Strategy A: Global brute-force combination (Find best length match)
    // Solves "ambiguous positioning" issues (e.g., both 130bp and 132bp candidates exist)
    // =========================================================
    if (found_L && found_R) {
        
        int best_diff = 0x7FFFFFFF; // Record minimum length deviation
        int found_valid_in_A = 0;   // Flag if a solution with valid length was found
        
        // Temporarily store best coordinates
        int best_sL = -1, best_eL = -1;
        int best_sR = -1, best_eR = -1;

        for (int i = 0; i < resL.numLocations; i++) {
            for (int j = 0; j < resR.numLocations; j++) {
                
                int tmp_sL = resL.startLocations[i];
                int tmp_eL = resL.endLocations[i];
                int tmp_sR = resR.startLocations[j];
                int tmp_eR = resR.endLocations[j];

                // Structural check: L must be before R
                if (tmp_sR <= tmp_eL) continue;

                int dist = tmp_sR - tmp_eL - 1;

                // Check if within allowed range (is_valid_length)
                if (is_valid_length(dist)) {
                    int current_diff = abs(dist - payload_size);
                    
                    // Greedy logic: If it's the first valid solution, or closer to payload_size than the previous one
                    if (!found_valid_in_A || current_diff < best_diff) {
                        best_diff = current_diff;
                        
                        best_sL = tmp_sL; best_eL = tmp_eL;
                        best_sR = tmp_sR; best_eR = tmp_eR;
                        
                        found_valid_in_A = 1;

                        // [Fast Optimization] If deviation is 0 (perfect), lock it immediately and break loop
                        if (best_diff == 0) {
                            start_L = best_sL; end_L = best_eL;
                            start_R = best_sR; end_R = best_eR;
                            solved_perfectly = 1;
                            goto FINISH_CHECK; // The use of goto here is professional and efficient
                        }
                    }
                }
            }
        }

        // If Strategy A found a valid solution (even if not 0 deviation), apply it
        if (found_valid_in_A) {
            start_L = best_sL; end_L = best_eL;
            start_R = best_sR; end_R = best_eR;
            solved_perfectly = 1;
            goto FINISH_CHECK; 
        }
    }

    // =========================================================
    // Strategy B: Local Re-search (Fallback)
    // Resolves false positive R caused by "crossing" or "too short" distances
    // Only execute if Strategy A did not solve the problem
    // =========================================================
    if (!solved_perfectly && found_L && found_R) 
    {
        int current_dist = start_R - end_L - 1;
        // Minimum reasonable distance floor (payload_size / 2)
        int min_expected_dist = payload_size / 2;

        // Trigger condition: R is before L, or R is too close to L
        if (start_R <= end_L || current_dist < min_expected_dist) 
        {
            // 1. Calculate deviation of the original result first (as a baseline for comparison)
            int old_diff = 0x7FFFFFFF;
            if (start_R > end_L) {
                old_diff = current_dist;
            }
            // If start_R <= end_L (crossing), old_diff remains max value, new result wins

            // Determine search start point
            int search_offset = 0;
            if (start_R > end_L) search_offset = end_R + 1; // Avoid R that is too close
            else search_offset = end_L + 1; // Avoid L

            if (search_offset < seqlen) {
                int search_len = seqlen - search_offset;

                if (search_len >= PRIMER_LEN) {
                    // Local re-search for R
                    EdlibAlignResult resR_2 = edlibAlign(primerR, PRIMER_LEN, line + search_offset, search_len, config);
                    
                    if (resR_2.numLocations > 0 && resR_2.editDistance != -1 && resR_2.editDistance <= MIN_EDIT) {
                        
                        int best_idx = -1;
                        int min_len_diff = 0x7FFFFFFF; // Minimum deviation in new search results

                        // Iterate through re-search results to find best match
                        for (int k = 0; k < resR_2.numLocations; k++) {
                            int tmp_start = resR_2.startLocations[k] + search_offset;
                            int tmp_dist = tmp_start - end_L - 1;

                            // Re-confirm positional relationship
                            if (tmp_start <= end_L) continue;

                            // Calculate deviation
                            int current_diff = tmp_dist;

                            // Find the best among current re-search results
                            // Prioritize is_valid_length; if none, find the one with minimum deviation
                            int is_perfect = is_valid_length(tmp_dist);
                            
                            // If a perfect one is found, or one better than current candidate
                            if (best_idx == -1 || is_perfect || current_diff < min_len_diff) {
                                // If a perfect one was already found and current one is not, skip
                                if (best_idx != -1 && is_valid_length(resR_2.startLocations[best_idx] + search_offset - end_L - 1) && !is_perfect) {
                                    continue;
                                }

                                best_idx = k;
                                min_len_diff = current_diff;
                            }
                        }

                        // [Critical Fix] Compare New vs Old results
                        // Only update if a valid new result is found AND it is closer to payload_size than the old result
                        if (best_idx != -1) {
                            if (min_len_diff < old_diff) {
                                start_R = resR_2.startLocations[best_idx] + search_offset;
                                end_R   = resR_2.endLocations[best_idx]   + search_offset;
                            }
                        }
                    }
                    edlibFreeAlignResult(resR_2);
                }
            }
        }
    }
FINISH_CHECK:

    // 3. Package result (Dual)
    if (found_L && found_R && start_R > end_L) {
        int len = start_R - end_L - 1;
        // As long as structure is Dual, mark it as Dual
        // Whether length is valid determines if it's saved in outer recursive_process
        res.valid = 1;
        res.is_dual = 1;
        res.start_idx = start_L; 
        res.end_idx = end_R;     
        res.payload_start = end_L + 1;
        res.payload_len = len;
        res.type_str = is_reverse ? "backward_dual" : "forward_dual";
    }

    // 4. Package result (Single - only if Dual is invalid)
    if (!res.valid) {
        if (found_L) {
            int target_start = end_L + 1;
            if (target_start + payload_size <= seqlen) {
                res.valid = 1; res.is_dual = 0;
                res.start_idx = start_L;
                res.end_idx = target_start + payload_size - 1; 
                res.payload_start = target_start; res.payload_len = payload_size;
                res.type_str = is_reverse ? "backward_single_L" : "forward_single_L";
            }
        } else if (found_R) {
            int target_end = start_R - 1;
            int target_start = target_end - payload_size + 1;
            if (target_start >= 0) {
                res.valid = 1; res.is_dual = 0;
                res.start_idx = target_start; res.end_idx = end_R;
                res.payload_start = target_start; res.payload_len = payload_size;
                res.type_str = is_reverse ? "backward_single_R" : "forward_single_R";
            }
        }
    }

    edlibFreeAlignResult(resL);
    edlibFreeAlignResult(resR);
    return res;
}

// Core recursive function
void recursive_process(const char *line, int seqlen, int read_num, 
                       int *match_counter, ThreadWorkerData *data)
{
    if (seqlen <= min_length) return;

    // 1. align both directions
    MatchResult fwd = check_alignment(data, line, seqlen, 0);
    MatchResult rev = check_alignment(data, line, seqlen, 1);

    MatchResult *choice = NULL;

    // 2. save priority decision logic
    
    // Check if both have perfect lengths
    int fwd_perfect = (fwd.valid && fwd.is_dual && is_valid_length(fwd.payload_len));
    int rev_perfect = (rev.valid && rev.is_dual && is_valid_length(rev.payload_len));

    // === Priority Tier 1: Perfect length dual (Preferred) ===
    if (fwd_perfect && rev_perfect) {
        // If both directions are perfect dual, choose Forward (default preference)
        choice = &fwd;
    }
    else if (fwd_perfect) {
        choice = &fwd;
    }
    else if (rev_perfect) {
        choice = &rev;
    }
    // === Priority Tier 2: Fallback Dual (Invalid length) ===
    // (Although invalid length won't be saved, we need to select it to enter recursion)
    else if (fwd.valid && fwd.is_dual) {
        choice = &fwd;
    }
    else if (rev.valid && rev.is_dual) {
        choice = &rev;
    }
    // === Priority Tier 3: Single end ===
    else if (fwd.valid) {
        choice = &fwd;
    }
    else if (rev.valid) {
        choice = &rev;
    }

    // 3. processing
    if (choice != NULL) {
        // Save logic: strictly check length
        if (is_valid_length(choice->payload_len))
        {
            (*match_counter)++;
            data->valid_match_count++;

            memcpy(data->temp_buf, line + choice->payload_start, choice->payload_len);
            
            if (choice->is_reverse) {
                inv(data->temp_buf, choice->payload_len);
                complement(data->temp_buf, data->temp_buf, choice->payload_len);
            }
            data->temp_buf[choice->payload_len] = '\0';

            fprintf(data->f_out_result, ">Read_%d_Match_%d_%s\n%s\n", 
                    read_num, *match_counter, choice->type_str, data->temp_buf);
            fprintf(data->f_out_len, "%d\n", choice->payload_len);
        }

        // Recursive logic
        // A. Recursive inside Payload
        if (choice->is_dual && choice->payload_len >= oligo_len) {
            recursive_process(line + choice->payload_start, choice->payload_len, read_num, match_counter, data);
        }

        // B. Recursive left region
        if (choice->start_idx > min_length) {
            recursive_process(line, choice->start_idx, read_num, match_counter, data);
        }

        // C. Recursive right region
        int right_start = choice->end_idx + 1;
        if (right_start < seqlen) {
            int right_len = seqlen - right_start;
            if (right_len > min_length) {
                recursive_process(line + right_start, right_len, read_num, match_counter, data);
            }
        }
    }
}

// ==== Worker thread logic ====
void* worker_routine(void *arg)
{
    ThreadWorkerData *data = (ThreadWorkerData*)arg;
    
    while (1) {
        TaskBatch *batch = NULL;

        // Lock and fetch task
        pthread_mutex_lock(&data->queue->mutex);
        // Wait for available tasks (using cond_empty)
        while (data->queue->head == NULL && !data->queue->shutdown) {
            pthread_cond_wait(&data->queue->cond_empty, &data->queue->mutex);
        }

        if (data->queue->head != NULL) {
            batch = data->queue->head;
            data->queue->head = batch->next;
            if (data->queue->head == NULL) data->queue->tail = NULL;
            data->queue->count--; // Consumed one, decrement count

            pthread_cond_signal(&data->queue->cond_full);
        } else if (data->queue->shutdown) {
            pthread_mutex_unlock(&data->queue->mutex);
            break;
        }
        pthread_mutex_unlock(&data->queue->mutex);

        // Process Batch
        if (batch) {
            for (int i = 0; i < batch->count; i++) {
                int match_count = 0;
                char *seq = batch->sequences[i];
                
                // Get original length
                int seq_len = strlen(seq);
                
                // loop to remove trailing \n and \r (compatible with Windows/Linux)
                while (seq_len > 0 && (seq[seq_len-1] == '\n' || seq[seq_len-1] == '\r')) {
                    seq[--seq_len] = '\0';
                }

                recursive_process(seq, seq_len, batch->read_ids[i], &match_count, data);
                
                // Count match_count distribution
                int idx = (match_count >= 100) ? 99 : match_count;
                data->local_match_counts[idx]++;
                
                data->processed_count++;
                
                free(batch->sequences[i]); 
            }
            free(batch->sequences);
            free(batch->read_ids);
            free(batch);
        }
    }
    
    // Close local files
    fclose(data->f_out_result);
    fclose(data->f_out_len);
    free(data->temp_buf);
    
    return NULL;
}

// ==== Main function ====

int main(int argc, char *argv[]) 
{
    if (argc != 12) {
        fprintf(stderr, "Usage: %s <fastq> <primer_list> <primer_num> <min_len> <match_cnt_file> <res_prefix> <valid_pct_file> <payload_size> <valid_base_size> <len_rec_file> <threads>\n", argv[0]);
        return 1;
    }

    char *fastq_file = argv[1];
    char *primer_list = argv[2];
    int primer_num = atoi(argv[3]);
    int min_length_param = atoi(argv[4]);
    char *reads_primer_match_counts = argv[5];
    char *result_output_prefix = argv[6];
    char *valid_base_percentage_file = argv[7];
    oligo_len = atoi(argv[8]);
    payload_size = min_length_param;
    valid_bias_size = atoi(argv[9]);
    
    char *length_record_file = argv[10];
    int num_threads = atoi(argv[11]);

    min_length = min_length_param; // Minimum length filter

    long nproc = sysconf(_SC_NPROCESSORS_ONLN);  // Number of available logical CPUs
    if (nproc < 1) {
        nproc = 1;   // Equivalent to std::max(1u, ...)
    }

    // num_threads = (unsigned int)nproc-5;
    // printf("\nnum_threads = %u\n", num_threads);

    // Read primers
    char **primers = (char **)malloc(primer_num * sizeof(char *));
    FILE *f_primer = fopen(primer_list, "r");
    if(!f_primer) { perror("Primer file"); return 1; }
    for(int i=0; i<primer_num; i++) {
        primers[i] = (char*)malloc(PRIMER_LEN + 1);
        int a = fscanf(f_primer, "%20s", primers[i]);
    }
    fclose(f_primer);


    // Precompute 4 primer sequences
    char p_fwd_L[PRIMER_LEN + 1];
    char p_fwd_R[PRIMER_LEN + 1];
    char p_rev_L[PRIMER_LEN + 1];
    char p_rev_R[PRIMER_LEN + 1];

    // 1. Forward (direct copy)
    memcpy(p_fwd_L, primers[0], PRIMER_LEN); p_fwd_L[PRIMER_LEN] = 0;
    memcpy(p_fwd_R, primers[1], PRIMER_LEN); p_fwd_R[PRIMER_LEN] = 0;

    // 2. Reverse
    // Rev L = RC of Primer 1
    char tempP[PRIMER_LEN + 1];
    memcpy(tempP, primers[1], PRIMER_LEN); 
    inv(tempP, PRIMER_LEN); 
    complement(tempP, p_rev_L, PRIMER_LEN); 
    p_rev_L[PRIMER_LEN] = 0;

    // Rev R = RC of Primer 0
    memcpy(tempP, primers[0], PRIMER_LEN); 
    inv(tempP, PRIMER_LEN); 
    complement(tempP, p_rev_R, PRIMER_LEN); 
    p_rev_R[PRIMER_LEN] = 0;

    // Initialize task queue
    TaskQueue queue = {0};
    queue_init(&queue, 500);

    // Start thread pool
    pthread_t *threads = (pthread_t*)malloc(sizeof(pthread_t) * num_threads);
    ThreadWorkerData *workers = (ThreadWorkerData*)malloc(sizeof(ThreadWorkerData) * num_threads);
    
    for(int i=0; i<num_threads; i++) 
    {
        workers[i].thread_id = i;
        workers[i].queue = &queue;
        // workers[i].primers = primers;
        // workers[i].primer_num = primer_num;
        workers[i].processed_count = 0;
        workers[i].valid_match_count = 0;
        workers[i].local_match_counts = (int*)calloc(100, sizeof(int));
        workers[i].temp_buf = (char*)malloc(MAXLEN);
        // Precompute primers into worker data
        memcpy(workers[i].p_fwd_L, p_fwd_L, PRIMER_LEN + 1);
        memcpy(workers[i].p_fwd_R, p_fwd_R, PRIMER_LEN + 1);
        memcpy(workers[i].p_rev_L, p_rev_L, PRIMER_LEN + 1);
        memcpy(workers[i].p_rev_R, p_rev_R, PRIMER_LEN + 1);

        // Each thread writes its own temporary file
        char fname[256];
        sprintf(fname, "%s_temp_%d", result_output_prefix, i);
        workers[i].f_out_result = fopen(fname, "w");
        // Set large buffer to reduce IO operations
        char *io_buf = (char*)malloc(64*1024);
        setvbuf(workers[i].f_out_result, io_buf, _IOFBF, 64*1024);
        
        sprintf(fname, "%s_len_%d", length_record_file, i);
        workers[i].f_out_len = fopen(fname, "w");
        
        pthread_create(&threads[i], NULL, worker_routine, &workers[i]);
    }

    // =========================================================
    // Auto-detect FASTA/FASTQ
    // =========================================================
    FILE *f_in = fopen(fastq_file, "r");
    if(!f_in) { perror("Input File"); return 1; }
    
    // 1. Auto-detect format
    int is_fasta = 0; // Default to FASTQ
    int first_char = fgetc(f_in); // Read the first character
    
    if (first_char == '>') {
        is_fasta = 1;
        // printf("Format detected: FASTA (Processing 2 lines per read)\n");
    } else if (first_char == EOF) {
        // printf("Warning: Input file is empty.\n");
    } else {
        // printf("Format detected: FASTQ (Processing 4 lines per read)\n");
    }

    // 2. Reset file pointer to the beginning
    rewind(f_in);

    char line[MAXLEN];
    int read_num = 0;
    
    TaskBatch *batch = (TaskBatch*)malloc(sizeof(TaskBatch));
    batch->sequences = (char**)malloc(sizeof(char*) * BATCH_SIZE);
    batch->read_ids = (int*)malloc(sizeof(int) * BATCH_SIZE);
    batch->count = 0;
    batch->next = NULL;

    // printf("Reading file...\n");
    
    // 3. reading loop
    while(fgets(line, MAXLEN, f_in)) // Line 1: ID (Header)
    { 
        // Line 2: Sequence
        if(!fgets(line, MAXLEN, f_in)) break; 
        
        read_num++;
        
        // Add to batch
        int idx = batch->count;
        batch->sequences[idx] = strdup(line);
        batch->read_ids[idx] = read_num;
        batch->count++;
        
        if(batch->count >= BATCH_SIZE) {
            pthread_mutex_lock(&queue.mutex);

            while (queue.count >= queue.max_size) {
                pthread_cond_wait(&queue.cond_full, &queue.mutex);
            }

            if(queue.tail == NULL) { queue.head = batch; queue.tail = batch; }
            else { queue.tail->next = batch; queue.tail = batch; }
            
            queue.count++;

            pthread_cond_signal(&queue.cond_empty);
            pthread_mutex_unlock(&queue.mutex);
    
            // Reset batch
            batch = (TaskBatch*)malloc(sizeof(TaskBatch));
            batch->sequences = (char**)malloc(sizeof(char*) * BATCH_SIZE);
            batch->read_ids = (int*)malloc(sizeof(int) * BATCH_SIZE);
            batch->count = 0;
            batch->next = NULL;
           
        }
        
        // Skip quality lines based on format
        if (!is_fasta) {
            // If FASTQ, consume the next two lines
            fgets(line, MAXLEN, f_in); // Line 3: +
            fgets(line, MAXLEN, f_in); // Line 4: Quality
        }

    }
    
    // Send the last Batch
    if(batch->count > 0) {
        pthread_mutex_lock(&queue.mutex);
        if(queue.tail == NULL) { queue.head = batch; queue.tail = batch; }
        else { queue.tail->next = batch; queue.tail = batch; }
        pthread_cond_signal(&queue.cond_empty);
        pthread_mutex_unlock(&queue.mutex);
    } else {
        free(batch->sequences); free(batch->read_ids); free(batch);
    }
    fclose(f_in);

    // Shutdown notification
    pthread_mutex_lock(&queue.mutex);
    queue.shutdown = 1;
    pthread_cond_broadcast(&queue.cond_empty);
    pthread_mutex_unlock(&queue.mutex);

    for(int i=0; i<num_threads; i++) pthread_join(threads[i], NULL);

    // ==== Merging results ====
    // printf("Merging results...\n");
    
    // 1. Merge result files
    char cmd[1024];
    // Use cat to merge all temp files into the final output
    sprintf(cmd, "cat %s_temp_* > %s && rm %s_temp_*", result_output_prefix, result_output_prefix, result_output_prefix);
    system(cmd);
    
    // 2. Merge length files
    sprintf(cmd, "cat %s_len_* > %s && rm %s_len_*", length_record_file, length_record_file, length_record_file);
    system(cmd);

    // 3. Aggregate statistics
    int total_valid_reads = 0;
    int global_match_counts[100] = {0};
    
    for(int i=0; i<num_threads; i++) {
        total_valid_reads += workers[i].valid_match_count;
        for(int j=0; j<100; j++) global_match_counts[j] += workers[i].local_match_counts[j];
        free(workers[i].local_match_counts);
    }

    // 4. Output statistics report
    FILE *f_stat = fopen(valid_base_percentage_file, "w");
    if(f_stat) {
        fprintf(f_stat, "Total reads: %d\n", read_num);
        fprintf(f_stat, "Valid extracted segments: %d\n", total_valid_reads);
        fprintf(f_stat, "Ratio: %.6f\n", (double)total_valid_reads/read_num);
        fclose(f_stat);
    }
    
    FILE *f_match = fopen(reads_primer_match_counts, "w");
    if(f_match) {
        for(int i=0; i<100; i++) fprintf(f_match, "%d\n", global_match_counts[i]);
        fclose(f_match);
    }
    // Free resources
    free(threads);
    free(workers);
    for(int i=0; i<primer_num; i++) free(primers[i]);
    free(primers);
    
    printf("  Total processed: %d, Valid segments: %d, percentage: %.6f\n", read_num, total_valid_reads, (double)total_valid_reads/read_num);
    return 0;
}
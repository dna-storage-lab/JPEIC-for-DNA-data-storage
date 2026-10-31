#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <sys/time.h>

// =====================================================================
// Parameter configuration
// =====================================================================
#define RS_N 4000
#define RS_K 3400
#define RS_M 12                  // GF(2^12)
#define NUM_GROUPS 250           // number of groups
#define RS_PER_GROUP 11          // codewords per group
#define BITS_PER_SYMBOL 12
#define PARITY_LEN (RS_N - RS_K) // parity symbols (600)
#define T_VAL (PARITY_LEN / 2)   // error-correction capability (300)
#define GF_SIZE 4096             // 2^12
#define init_zero 0

// =====================================================================
// Global read-only RS tables (shared across threads)
// =====================================================================
int *alpha_to, *index_of, *g;
int pp[RS_M + 1];                // primitive polynomial coefficients (replaces p[])

// =====================================================================
// Work scheduling (thread-safe group assignment)
// =====================================================================
int global_next_group_id = 0;    // next group index to process
pthread_mutex_t job_mutex;       // protects global_next_group_id

// Decoding status: NUM_GROUPS × RS_PER_GROUP
// 0 = not yet decoded, 1 = successfully decoded
int global_rs_status[NUM_GROUPS][RS_PER_GROUP];
pthread_mutex_t status_mutex;    // protects status updates

// =====================================================================
// Thread-private workspace
// =====================================================================
typedef struct
{
    int thread_id;
    char input_dir[256];
    char output_dir[256];

    // --- RS decoder working arrays ---
    int *recd;       // received codeword (full-length buffer)
    int *s;          // syndromes

    // Berlekamp–Massey / Forney scratch buffers
    int *d;          // discrepancy
    int *l;          // ELP degree
    int *u_lu;       // auxiliary for degree update
    int *forney;     // Forney-modified syndromes
    int *err;        // error magnitudes
    int *omega;      // error evaluator polynomial
    int *phi;        // errata locator polynomial
    int *phiprime;   // derivative of phi
    int *tau;        // erasure locator polynomial
    int *root;       // roots of the locator polynomial
    int *loc;        // error/erasure locations
    int *z;          // scratchpad
    int *reg;        // register for Chien search
    int *aux;        // auxiliary scratch

    // Flattened storage for the ELP matrix (row-major)
    // Backed by one contiguous array and sliced into per-row pointers.
    int *elp_storage;

    // --- I/O buffers ---
    int *data_buf;   // symbol matrix buffer
    char *out_str;   // output bitstream buffer

} ThreadContext;

// =====================================================================
// Load decoding status from disk (resume support)
// =====================================================================
void load_status(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        // Status file not found: treat as a fresh run
        memset(global_rs_status, 0, sizeof(global_rs_status));
        return;
    }

    // Read status entries (NUM_GROUPS × RS_PER_GROUP)
    for (int g = 0; g < NUM_GROUPS; g++) {
        for (int c = 0; c < RS_PER_GROUP; c++)
        {
            if (fscanf(fp, "%d", &global_rs_status[g][c]) != 1)
            {
                // Truncated/corrupted status file: keep the remaining entries as 0
                fprintf(stderr,
                        "Warning: Status file corrupted or truncated at Group %d, Code %d. Stopping load.\n",
                        g, c);
                fclose(fp);
                return;
            }
        }
    }
    fclose(fp);
}

// =====================================================================
// Save decoding status to disk
// =====================================================================
void save_status(const char *filename)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) return;

    for (int g = 0; g < NUM_GROUPS; g++) {
        for (int c = 0; c < RS_PER_GROUP; c++) {
            fprintf(fp, "%d ", global_rs_status[g][c]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

// =====================================================================
// Function declarations
// =====================================================================
void read_p();
void generate_gf();
void gen_poly();
int decode_rs_block(int *input_data, int *decoded_data, ThreadContext *ctx);

// =====================================================================
// GF initialization and generator polynomial construction
// =====================================================================
void read_p()
{
    // Primitive polynomial for GF(2^12):
    //   p(x) = 1 + x^3 + x^4 + x^7 + x^12
    memset(pp, 0, sizeof(pp));
    pp[0] = pp[RS_M] = 1;
    if (RS_M == 12) { pp[3] = pp[4] = pp[7] = 1; }
}

void generate_gf()
{
    int i, mask = 1;

    // Lookup tables:
    //   alpha_to[i] = α^i in polynomial form
    //   index_of[x] = log_α(x) in index form (index_of[0] = -1)
    alpha_to = (int*)malloc(GF_SIZE * sizeof(int));
    index_of = (int*)malloc(GF_SIZE * sizeof(int));

    alpha_to[RS_M] = 0;
    for (i = 0; i < RS_M; i++) {
        alpha_to[i] = mask;
        index_of[alpha_to[i]] = i;
        if (pp[i] != 0) alpha_to[RS_M] ^= mask;
        mask <<= 1;
    }

    index_of[alpha_to[RS_M]] = RS_M;
    mask >>= 1;

    for (i = RS_M + 1; i < GF_SIZE - 1; i++) {
        if (alpha_to[i - 1] >= mask)
            alpha_to[i] = alpha_to[RS_M] ^ ((alpha_to[i - 1] ^ mask) << 1);
        else
            alpha_to[i] = alpha_to[i - 1] << 1;

        index_of[alpha_to[i]] = i;
    }

    index_of[0] = -1;
}

void gen_poly()
{
    int i, j;

    // Generator polynomial g(x) of degree PARITY_LEN
    // Coefficients are stored in log form at the end.
    g = (int*)malloc((PARITY_LEN + 1) * sizeof(int));

    // Initialize g(x) = (x + 1)
    g[0] = 1;
    g[1] = 1;

    // Multiply by (x + α^(i-1)) for i = 2..PARITY_LEN
    for (i = 2; i <= PARITY_LEN; i++) {
        g[i] = 1;
        for (j = i - 1; j > 0; j--)
            if (g[j] != 0)
                g[j] = g[j - 1] ^ alpha_to[(index_of[g[j]] + i - 1) % (GF_SIZE - 1)];
            else
                g[j] = g[j - 1];

        g[0] = alpha_to[(index_of[g[0]] + i - 1) % (GF_SIZE - 1)];
    }

    // Convert to log form
    for (i = 0; i <= PARITY_LEN; i++)
        g[i] = index_of[g[i]];
}

// =====================================================================
// Core RS decoding routine (thread-local buffers, no malloc; shortened code
// with +95 left-padding compensation).
//
// Returns:
//   0   on success (syndrome check passed)
//  -1   on failure (syndrome check failed)
// =====================================================================
int decode_rs_block(int *input_data, int *decoded_data, ThreadContext *ctx) 
{
    // =========================================================
    // 1. Thread-local buffer mapping (no dynamic allocation)
    // =========================================================
    // Reuse preallocated arrays inside ThreadContext to avoid malloc/free.
    int *d        = ctx->d;
    int *l        = ctx->l;
    int *u_lu     = ctx->u_lu;
    int *s        = ctx->s;
    int *forney   = ctx->forney;
    int *err      = ctx->err;
    int *omega    = ctx->omega;
    int *phi      = ctx->phi;
    int *phiprime = ctx->phiprime;
    int *tau      = ctx->tau;
    int *root     = ctx->root;
    int *loc      = ctx->loc;
    int *z        = ctx->z;
    int *reg      = ctx->reg;
    int *aux      = ctx->aux;

    // Local scratch to replace original code's pointer usage
    int p_local[2];
    int *p = p_local;

    // Alias for received vector (to keep original recd[j] style)
    int *recd = ctx->recd;

    // 2D ELP matrix view backed by a flattened storage buffer
    const int STRIDE = PARITY_LEN + 20;
    int *elp[PARITY_LEN + 20];
    for(int k = 0; k < PARITY_LEN + 20; k++) {
        elp[k] = ctx->elp_storage + k * STRIDE;
    }

    // Clear full-length working buffers (RS over GF(2^12): length 4096)
    memset(recd, 0, 4096 * sizeof(int));
    memset(err,  0, 4096 * sizeof(int));

    int i, j, u, q;
    int count = 0, syn_error = 0;
    int degphi, ell, temp = 0;
    int n = GF_SIZE - 1;          // n = 4095 (RS full length without the "0" symbol)
    int k = n - PARITY_LEN;
    int padding = 95;             // 4095 - 4000

    // =========================================================
    // 2. Shortening compensation and erasure collection
    // =========================================================
    // A) Left-pad 95 zeros for shortened-code alignment.
    for (i = 0; i < padding; i++) {
        recd[i] = 0;
    }

    // B) Shift input symbols [0..RS_N-1] -> [padding..padding+RS_N-1]
    //    and collect erasure locations (input symbol = -1).
    int numera = 0;
    int erasure_locs[PARITY_LEN + 2];
    int *era = erasure_locs;

    for (i = 0; i < RS_N; i++)            // RS_N = 4000
    {
        int target_idx = i + padding;     // +95 offset

        if (input_data[i] == -1) {
            // Erasure: insert placeholder and record position
            recd[target_idx] = 0;
            era[numera++] = target_idx;
        } else {
            recd[target_idx] = input_data[i];
        }
    }

    // =========================================================
    // 3. Syndrome computation on the full-length vector (n = 4095)
    // =========================================================
    for (i = 1; i <= PARITY_LEN; i++)
    {
        s[i] = 0;

        // Full scan over 0..n-1 is kept for generality.
        // (An optional optimization is scanning padding..n-1 only.)
        for (j = 0; j < n; j++)
        {
            if (recd[j] != 0) {
                s[i] ^= alpha_to[(index_of[recd[j]] + (i - 1) * j) % n];
            }
        }

        if (s[i] != 0) syn_error = 1;
        s[i] = index_of[s[i]];
    }

    // Fast path: no syndrome and no erasures => already a valid codeword
    if (!syn_error && numera == 0) {
        memcpy(decoded_data, recd + padding, RS_N * sizeof(int));
        return 0;
    }

    // =========================================================
    // 4. Berlekamp–Massey decoding (errors + erasures)
    // =========================================================
    if (syn_error)   // attempt correction only if syndromes are non-zero
    {
        s[0] = 0;    // S(x) = 1 + s1 x + ...

        // -------------------------------
        // Erasure handling (tau / Forney syndromes)
        // -------------------------------
        if (numera)
        {
            // Compute erasure locator polynomial tau(x)
            for (i = 0; i <= PARITY_LEN; i++) {
                tau[i] = 0;
                aux[i] = 0;
            }

            aux[1] = alpha_to[era[0]];
            aux[0] = 1;   // (x + era[0])

            if (numera > 1)
                for (i = 0; i < numera; i++)
                {
                    p[1] = era[i];
                    p[0] = 0;

                    for (j = 0; j < 2; j++)
                        for (ell = 0; ell <= i; ell++)
                            if ((p[j] != -1) && (aux[ell] != 0))
                                tau[j + ell] ^=
                                    alpha_to[(p[j] + index_of[aux[ell]]) % n];

                    if (i != (numera - 1))
                        for (ell = 0; ell <= (i + 1); ell++) {
                            aux[ell] = tau[ell];
                            tau[ell] = 0;
                        }
                }
            else {
                tau[0] = aux[0];
                tau[1] = aux[1];
            }

            // Convert tau(x) to log form
            for (i = 0; i <= numera; i++)
                tau[i] = index_of[tau[i]];

            // Compute Forney-modified syndromes:
            //   forney(x) = [ s(x) * tau(x) + 1 ] mod x^{PARITY_LEN}
            for (i = 0; i <= n - k; i++) forney[i] = 0;
            for (i = 0; i <= n - k; i++)
                for (j = 0; j <= numera; j++)
                    if (i + j <= (n - k))
                        if ((s[i] != -1) && (tau[j] != -1))
                            forney[i + j] ^= alpha_to[(s[i] + tau[j]) % n];

            forney[0] ^= 1;

            for (i = 0; i <= n - k; i++)
                forney[i] = index_of[forney[i]];
        }
        else
        {
            // No erasures: Forney syndromes reduce to standard syndromes
            tau[0] = 0;
            for (i = 1; i <= n - k; i++) forney[i] = s[i];
        }

        // --------------------------------------------------------------
        // Berlekamp–Massey iterations (errors + erasures)
        // --------------------------------------------------------------

        // Initialize discrepancy and ELP table
        d[0] = 0;                   
        d[1] = forney[numera + 1];  
        elp[0][0] = 0;              
        elp[1][0] = 1;              
        for (i = 1; i < PARITY_LEN; i++) {
            elp[0][i] = -1;
            elp[1][i] = 0;
        }

        l[0] = 0;
        l[1] = 0;
        u_lu[0] = -1;
        u_lu[1] = 0;
        u = 0;

        if (numera < PARITY_LEN)
        {
            do
            {
                u++;

                if (d[u] == -1)
                {
                    // Zero discrepancy: copy previous ELP
                    l[u + 1] = l[u];
                    for (i = 0; i <= l[u]; i++) {
                        elp[u + 1][i] = elp[u][i];
                        elp[u][i] = index_of[elp[u][i]];
                    }
                }
                else
                {
                    // Non-zero discrepancy: find q with maximal u_lu[q] and d[q] != 0
                    q = u - 1;
                    while ((d[q] == -1) && (q > 0)) q--;

                    if (q > 0) {
                        j = q;
                        do {
                            j--;
                            if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
                                q = j;
                        } while (j > 0);
                    }

                    // Update ELP degree
                    if (l[u] > l[q] + u - q)  l[u + 1] = l[u];
                    else                      l[u + 1] = l[q] + u - q;

                    // Compute new ELP polynomial
                    for (i = 0; i < PARITY_LEN; i++) elp[u + 1][i] = 0;
                    for (i = 0; i <= l[q]; i++)
                        if (elp[q][i] != -1)
                            elp[u + 1][i + u - q] =
                                alpha_to[(d[u] + n - d[q] + elp[q][i]) % n];

                    for (i = 0; i <= l[u]; i++) {
                        elp[u + 1][i] ^= elp[u][i];
                        elp[u][i] = index_of[elp[u][i]];
                    }
                }

                u_lu[u + 1] = u - l[u + 1];

                // Compute next discrepancy (skip last iteration)
                if (u < (PARITY_LEN - numera))
                {
                    if (forney[numera + u + 1] != -1)
                        d[u + 1] = alpha_to[forney[numera + u + 1]];
                    else
                        d[u + 1] = 0;

                    for (i = 1; i <= l[u + 1]; i++)
                        if ((forney[numera + u + 1 - i] != -1) && (elp[u + 1][i] != 0))
                            d[u + 1] ^= alpha_to[(forney[numera + u + 1 - i]
                                + index_of[elp[u + 1][i]]) % n];

                    d[u + 1] = index_of[d[u + 1]];
                }

            } while ((u < (PARITY_LEN - numera)) &&
                     (l[u + 1] <= ((PARITY_LEN - numera) / 2)));
        }

        u++;

        if (l[u] <= T_VAL - numera / 2)
        {
            // Convert final ELP coefficients to log form
            for (i = 0; i <= l[u]; i++) elp[u][i] = index_of[elp[u][i]];

            // Find roots of error-location polynomial
            for (i = 1; i <= l[u]; i++) reg[i] = elp[u][i];
            count = 0;

            for (i = 1; i <= n; i++)
            {
                q = 1;
                for (j = 1; j <= l[u]; j++)
                    if (reg[j] != -1) {
                        reg[j] = (reg[j] + j) % n;
                        q ^= alpha_to[reg[j]];
                    }

                if (!q)
                {
                    root[count] = i;
                    loc[count] = n - i;
                    count++;
                }
            }

            if (count == l[u])
            {
                // Compute omega(x): error evaluator polynomial
                forney[0] = 0;
                for (i = 0; i <= PARITY_LEN; i++) omega[i] = 0;

                for (i = 0; i <= PARITY_LEN; i++)
                    for (j = 0; j <= l[u]; j++)\
                        if (i + j <= PARITY_LEN)
                            if ((forney[i] != -1) && (elp[u][j] != -1))
                                omega[i + j] ^= alpha_to[(forney[i] + elp[u][j]) % n];

                for (i = 0; i <= PARITY_LEN; i++)
                    omega[i] = index_of[omega[i]];

                // Compute phi(x): errata locator polynomial
                degphi = numera + l[u];
                for (i = 0; i <= degphi; i++) phi[i] = 0;

                for (i = 0; i <= numera; i++)
                    for (j = 0; j <= l[u]; j++)
                        if ((tau[i] != -1) && (elp[u][j] != -1))
                            phi[i + j] ^= alpha_to[(tau[i] + elp[u][j]) % n];

                for (i = 0; i <= degphi; i++)
                    phi[i] = index_of[phi[i]];

                // Compute derivative phi'(x): odd-power terms only
                for (i = 0; i <= degphi; i++) phiprime[i] = -1;
                for (i = 0; i <= degphi; i++)
                    if (i % 2) phiprime[i - 1] = phi[i];

                // If erasures exist, append them into (loc, root)
                if (numera)
                    for (i = 0; i < numera; i++) {
                        loc[count + i] = era[i];
                        root[count + i] = (n - era[i]) % n;
                    }

                // Forney evaluation: compute error magnitudes and correct recd[]
                for (i = 0; i < degphi; i++)
                {
                    err[loc[i]] = 0;

                    for (j = 0; j <= PARITY_LEN; j++)
                        if ((omega[j] != -1) && (root[i] != -1))
                            err[loc[i]] ^= alpha_to[(omega[j] + j * root[i]) % n];

                    // Multiply by loc[i]^{2-init_zero} term as in the original implementation
                    if ((err[loc[i]] != 0) && (loc[i] != -1))
                        err[loc[i]] = alpha_to[(index_of[err[loc[i]]]
                            + loc[i] * (2 - init_zero + n)) % n];

                    if (err[loc[i]] != 0)
                    {
                        err[loc[i]] = index_of[err[loc[i]]];

                        q = 0;
                        for (j = 0; j <= degphi; j++)
                            if ((phiprime[j] != -1) && (root[i] != -1))
                                q ^= alpha_to[(phiprime[j] + j * root[i]) % n];

                        // Divide by q and apply correction
                        err[loc[i]] = alpha_to[(err[loc[i]] - index_of[q] + n) % n];
                        recd[loc[i]] ^= err[loc[i]];
                    }
                }
            }
            else
            {
                // Root count mismatch => decoding failure (> t errors)
            }
        }
        else
        {
            // Degree too large => decoding failure (> t errors)
        }
    }
    else
    {
        // All syndromes zero: no correction needed
    }

    // =========================================================
    // 5. Final syndrome verification (post-correction)
    // =========================================================
    // Recompute syndromes on corrected recd[] to validate codeword legality.
    int check_success = 1;

    if (syn_error)
    {
        for (i = 1; i <= PARITY_LEN; i++)
        {
            int final_syn = 0;

            for (j = 0; j < n; j++)
            {
                if (recd[j] != 0)
                {
                    int idx = (index_of[recd[j]] + (i - 1) * j) % n;
                    final_syn ^= alpha_to[idx];
                }
            }

            if (final_syn != 0)
            {
                check_success = 0;
                break;
            }
        }
    }

    // =========================================================
    // 6. Output extraction (remove left-padding)
    // =========================================================
    memcpy(decoded_data, recd + padding, RS_N * sizeof(int));

    // Return status:
    //   0  success
    //  -1  failure
    if (check_success) {
        return 0;
    } else {
        // Optionally, one may revert decoded_data back to input_data here.
        return -1;
    }
}

// =====================================================================
// Core worker-thread routine with cache-aware decoding
// =====================================================================
void* worker_thread(void *arg)
{
    ThreadContext *ctx = (ThreadContext*)arg;

    // Reuse preallocated buffer as a large RS symbol matrix:
    //   RS_PER_GROUP rows (codewords) × RS_N columns (symbols)
    // Assumes ctx->data_buf is sufficiently large.
    int *codeword_matrix = ctx->data_buf;

    while (1)
    {
        int my_group_id = -1;

        // -----------------------------------------------------------------
        // Fetch next unprocessed group ID (thread-safe)
        // -----------------------------------------------------------------
        pthread_mutex_lock(&job_mutex);
        if (global_next_group_id < NUM_GROUPS) {
            my_group_id = global_next_group_id;
            global_next_group_id++;
        }
        pthread_mutex_unlock(&job_mutex);

        // No remaining work
        if (my_group_id == -1) {
            break;
        }

        // -----------------------------------------------------------------
        // Construct file paths
        // -----------------------------------------------------------------
        char in_path[512], cache_path[512], out_path[512];

        // Raw noisy input
        sprintf(in_path, "%s/group_%d.txt", ctx->input_dir, my_group_id);
        // Cache file storing corrected codewords from previous rounds
        sprintf(cache_path, "%s/cache_group_%d.txt", ctx->output_dir, my_group_id);
        // Final decoded bitstream output
        sprintf(out_path, "%s/decoded_%d.txt", ctx->output_dir, my_group_id);

        // -----------------------------------------------------------------
        // Step A: Load noisy input codewords
        // -----------------------------------------------------------------
        FILE *fp_in = fopen(in_path, "r");
        if (!fp_in) {
            fprintf(stderr, "Error opening %s\n", in_path);
            continue;
        }

        for (int i = 0; i < RS_PER_GROUP; i++) {
            for (int j = 0; j < RS_N; j++) {
                // Read symbols into row i of the matrix
                if (fscanf(fp_in, "%d", &codeword_matrix[i * RS_N + j]) != 1)
                {
                    // On read failure (EOF or format error),
                    // explicitly mark as erasure (-1) to avoid stale data
                    codeword_matrix[i * RS_N + j] = -1;
                }
            }
        }
        fclose(fp_in);

        // -----------------------------------------------------------------
        // Step B: Load cached corrected codewords (if available)
        // -----------------------------------------------------------------
        FILE *fp_cache = fopen(cache_path, "r");
        if (fp_cache)
        {
            int temp_val;
            for (int i = 0; i < RS_PER_GROUP; i++)
            {
                // Check whether this codeword was already successfully decoded
                int is_already_success = global_rs_status[my_group_id][i];

                for (int j = 0; j < RS_N; j++)
                {
                    if (fscanf(fp_cache, "%d", &temp_val) == 1)
                    {
                        // Only overwrite input with cached data
                        // if the codeword was previously decoded successfully
                        if (is_already_success)
                        {
                            codeword_matrix[i * RS_N + j] = temp_val;
                        }
                    }
                }
            }
            fclose(fp_cache);
        }

        // -----------------------------------------------------------------
        // Step C: Decode RS codewords within this group
        // -----------------------------------------------------------------
        for (int i = 0; i < RS_PER_GROUP; i++)
        {
            // Skip decoding if already marked as successful
            if (global_rs_status[my_group_id][i] == 1)
            {
                // Correct data already loaded from cache
                continue;
            }

            // Pointer to the current codeword (row i)
            int *current_codeword = codeword_matrix + i * RS_N;

            // Perform RS decoding in-place
            int result = decode_rs_block(current_codeword, current_codeword, ctx);

            if (result == 0) {
                // Decoding succeeded: update global status
                pthread_mutex_lock(&status_mutex);
                global_rs_status[my_group_id][i] = 1;
                pthread_mutex_unlock(&status_mutex);
            }
        }

        // -----------------------------------------------------------------
        // Step D: Save updated cache and decoded output
        // -----------------------------------------------------------------

        // (1) Update cache file
        // Always write back the latest matrix state:
        //   - successfully decoded codewords are preserved
        //   - failed ones remain as-is for the next round
        FILE *fp_cache_write = fopen(cache_path, "w");
        if (fp_cache_write) {
            char write_buf[128 * 1024];
            setvbuf(fp_cache_write, write_buf, _IOFBF, sizeof(write_buf));

            for (int i = 0; i < RS_PER_GROUP; i++)
            {
                for (int j = 0; j < RS_N; j++)
                {
                    fprintf(fp_cache_write, "%d ",
                            codeword_matrix[i * RS_N + j]);
                }
                fprintf(fp_cache_write, "\n");
            }
            fclose(fp_cache_write);
        }

        // (2) Write interleaved decoded bitstream
        FILE *fp_out = fopen(out_path, "w");
        if (fp_out) {
            char file_buf[256 * 1024];
            setvbuf(fp_out, file_buf, _IOFBF, sizeof(file_buf));

            // out_str must be sufficiently large
            char *ptr = ctx->out_str;

            // Interleaving order:
            //   symbol 0 of all codewords,
            //   symbol 1 of all codewords, ...
            for (int k = 0; k < RS_K; k++)
            {
                for (int i = 0; i < RS_PER_GROUP; i++)
                {
                    // Use data symbols only (exclude parity region if present)
                    int val = codeword_matrix[i * RS_N + k];

                    // Safety fallback: treat erasure as zero
                    if (val == -1) val = 0;

                    // Output as 12-bit binary representation (MSB first)
                    for (int b = 11; b >= 0; b--) {
                        *ptr++ = ((val >> b) & 1) ? '1' : '0';
                    }
                }
            }
            *ptr++ = '\n';
            *ptr = '\0';

            fputs(ctx->out_str, fp_out);
            fclose(fp_out);
        }
    }

    return NULL;
}


int main(int argc, char *argv[])
{
    // ---------------------------------------------------------------------
    // Command-line arguments:
    //   argv[1]: input directory containing group-wise input files
    //   argv[2]: output directory for per-group decoded results
    //   argv[3]: number of worker threads
    //   argv[4]: status file for checkpointing decoding progress
    //   argv[5]: flag file indicating overall decoding success
    //   argv[6]: final merged output file
    // ---------------------------------------------------------------------
    if (argc != 7) {
        printf("Usage: %s <input_dir> <output_dir> <num_threads> <status_file> <flag_file> <final_output>\n", argv[0]);
        return 1;
    }
    
    char *input_dir = argv[1];
    char *output_dir = argv[2];
    int num_threads = atoi(argv[3]);
    
    // Path to persistent status file (for resume / checkpointing)
    char *status_file_path = argv[4];
    
    // Output files
    char *flag_file_path = argv[5];
    char *final_output_path = argv[6];
    
    // ---------------------------------------------------------------------
    // Reed–Solomon parameter-dependent buffer sizes
    // ---------------------------------------------------------------------

    // Full codeword length for RS over GF(2^12): indices 0–4095
    const int SIZE_N = 4096;

    // Maximum polynomial degree involved in syndrome / error locator
    // (parity length = 600, with additional margin to avoid overflow)
    const int SIZE_2T = PARITY_LEN + 20;
    
    // Number of iterations (rows) for Berlekamp–Massey ELP matrix
    const int ELP_ROWS = SIZE_2T;

    // ---------------------------------------------------------------------
    // Initialize global RS lookup tables
    // ---------------------------------------------------------------------
    read_p();
    generate_gf();
    gen_poly();
    
    // ---------------------------------------------------------------------
    // Initialize global synchronization primitives and decoding status
    // ---------------------------------------------------------------------
    pthread_mutex_init(&job_mutex, NULL);
    global_next_group_id = 0;

    pthread_mutex_init(&status_mutex, NULL);
    load_status(status_file_path);

    // ---------------------------------------------------------------------
    // Launch worker thread pool
    // ---------------------------------------------------------------------
    pthread_t *threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    ThreadContext *contexts = (ThreadContext*)malloc(num_threads * sizeof(ThreadContext));
    
    for (int i = 0; i < num_threads; i++)
    {
        contexts[i].thread_id = i;
        strcpy(contexts[i].input_dir, input_dir);
        strcpy(contexts[i].output_dir, output_dir);

        // Buffer for decoded RS symbols (RS_PER_GROUP × RS_N)
        contexts[i].data_buf = (int*)malloc((RS_PER_GROUP * RS_N + 10) * sizeof(int));
        
        // Output bitstream buffer (RS_PER_GROUP × RS_K × 12 bits)
        contexts[i].out_str = (char*)malloc(RS_PER_GROUP * RS_K * 12 + 10);

        // -----------------------------------------------------------------
        // A. Full-length codeword buffers (size = 4096)
        // -----------------------------------------------------------------
        contexts[i].recd      = (int*)malloc(SIZE_N * sizeof(int)); // received symbols
        contexts[i].err       = (int*)malloc(SIZE_N * sizeof(int)); // error values
        contexts[i].omega     = (int*)malloc(SIZE_N * sizeof(int)); // error evaluator
        contexts[i].phi       = (int*)malloc(SIZE_N * sizeof(int)); // error locator
        contexts[i].phiprime  = (int*)malloc(SIZE_N * sizeof(int)); // derivative of locator

        // -----------------------------------------------------------------
        // B. Polynomial buffers (degree ≤ 2t)
        // -----------------------------------------------------------------
        contexts[i].s      = (int*)malloc(SIZE_2T * sizeof(int)); // syndromes
        contexts[i].d      = (int*)malloc(SIZE_2T * sizeof(int)); // discrepancy
        contexts[i].l      = (int*)malloc(SIZE_2T * sizeof(int)); // polynomial degree
        contexts[i].u_lu   = (int*)malloc(SIZE_2T * sizeof(int)); // auxiliary variables
        contexts[i].forney = (int*)malloc(SIZE_2T * sizeof(int)); // Forney algorithm
        contexts[i].tau    = (int*)malloc(SIZE_2T * sizeof(int));
        contexts[i].root   = (int*)malloc(SIZE_2T * sizeof(int)); // error roots
        contexts[i].loc    = (int*)malloc(SIZE_2T * sizeof(int)); // error locations
        contexts[i].z      = (int*)malloc(SIZE_2T * sizeof(int));
        contexts[i].reg    = (int*)malloc(SIZE_2T * sizeof(int));
        contexts[i].aux    = (int*)malloc(SIZE_2T * sizeof(int));

        // -----------------------------------------------------------------
        // C. Error locator polynomial (ELP) storage
        // Each row corresponds to one BM iteration
        // -----------------------------------------------------------------
        contexts[i].elp_storage = (int*)malloc(ELP_ROWS * SIZE_2T * sizeof(int));

        // Sanity check
        if (!contexts[i].recd || !contexts[i].s || !contexts[i].elp_storage) {
            fprintf(stderr, "Critical: Memory allocation failed for thread %d\n", i);
            exit(1);
        }

        pthread_create(&threads[i], NULL, worker_thread, (void*)&contexts[i]);
    }

    // ---------------------------------------------------------------------
    // Wait for all worker threads and release per-thread resources
    // ---------------------------------------------------------------------
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
        
        // I/O buffers
        free(contexts[i].data_buf);
        free(contexts[i].out_str);

        // Group A buffers
        free(contexts[i].recd);
        free(contexts[i].err);
        free(contexts[i].omega);
        free(contexts[i].phi);
        free(contexts[i].phiprime);

        // Group B buffers
        free(contexts[i].s);
        free(contexts[i].d);
        free(contexts[i].l);
        free(contexts[i].u_lu);
        free(contexts[i].forney);
        free(contexts[i].tau);
        free(contexts[i].root);
        free(contexts[i].loc);
        free(contexts[i].z);
        free(contexts[i].reg);
        free(contexts[i].aux);

        // Group C buffer
        free(contexts[i].elp_storage);
    }

    free(threads);
    free(contexts);
    pthread_mutex_destroy(&job_mutex);
    pthread_mutex_destroy(&status_mutex);
    
    // ---------------------------------------------------------------------
    // Save decoding status and summarize global decoding results
    // ---------------------------------------------------------------------
    save_status(status_file_path);

    int total_success = 0;
    int total_codewords = NUM_GROUPS * RS_PER_GROUP;
    
    for (int g = 0; g < NUM_GROUPS; g++) {
        for (int c = 0; c < RS_PER_GROUP; c++) {
            if (global_rs_status[g][c] == 1) total_success++;
        }
    }
    
    printf("  Successfully decoded Reed-Solomon codewords: %d / %d\n",
           total_success, total_codewords);
    
    // Write success flag
    FILE *fp_flag = fopen(flag_file_path, "w");
    if (fp_flag) {
        fprintf(fp_flag, (total_success == total_codewords) ? "yes\n" : "no\n");
        fclose(fp_flag);
    } else {
        fprintf(stderr, "Error: Cannot open flag file %s\n", flag_file_path);
    }
    
    // ---------------------------------------------------------------------
    // Merge per-group decoded outputs into a single file
    // ---------------------------------------------------------------------
    FILE *fp_final = fopen(final_output_path, "w");
    if (fp_final) 
    {
        // Large I/O buffer to minimize system calls during merging
        const size_t BUF_SIZE = 512 * 1024;
        char *copy_buf = (char*)malloc(BUF_SIZE);
        if (!copy_buf) {
            perror("Malloc merge buffer");
            fclose(fp_final);
            return 1;
        }

        setvbuf(fp_final, NULL, _IOFBF, 1024 * 1024);

        char part_filename[1024];

        for (int g = 0; g < NUM_GROUPS; g++) 
        {
            sprintf(part_filename, "%s/decoded_%d.txt", output_dir, g);
            
            FILE *fp_part = fopen(part_filename, "r");
            if (fp_part) 
            {
                size_t n;
                while ((n = fread(copy_buf, 1, BUF_SIZE, fp_part)) > 0) {
                    fwrite(copy_buf, 1, n, fp_final);
                }
                fclose(fp_part);
            } 
            else 
            {
                // Missing group output: insert empty line to preserve alignment
                fprintf(stderr, "Warning: Missing partial file %s\n", part_filename);
                fputc('\n', fp_final);
            }
        }
        
        free(copy_buf); 
        fclose(fp_final);
    } 
    else 
    {
        fprintf(stderr, "Error: Cannot open final output file %s\n", final_output_path);
    }

    return 0;
}
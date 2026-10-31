#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>

#define MAX_ADDR 4000
#define BIT_LENGTH 156
#define TARGET_START 0
#define TARGET_END 11
#define TARGET_LENGTH (TARGET_END - TARGET_START + 1)
#define RESULT_LENGTH 144

// Voting structure to store the count of 0 and 1 votes for each bit position
typedef struct {
    int vote_0[RESULT_LENGTH];     // Count of '0' votes for each bit position
    int vote_1[RESULT_LENGTH];     // Count of '1' votes for each bit position
    int total_votes;               // Total votes for this address
} VoteResult;

// Convert binary string to decimal
uint32_t bin_to_dec(const char *bin_str, int start, int end)
{
    uint32_t result = 0;
    for (int i = start; i <= end; i++) 
    {
        result = (result << 1) | (bin_str[i] - '0');
    }
    return result;
}

void create_directory(const char *dirname) 
{
    // 0777 gives read/write/execute permission to all users
    if (mkdir(dirname, 0777) != 0) {
        if (errno != EEXIST) {  // If error is not "directory already exists"
            perror(dirname);
            exit(EXIT_FAILURE);
        }
    }
}

int file_exists(const char *filename) 
{
    FILE *file = fopen(filename, "r");
    if (file) {
        fclose(file);
        return 1;
    }
    return 0;
}


char dexor(char a, char b);
char mapping(char a, char b);

int main(int argc, char *argv[]) 
{
    if (argc != 8) 
    {
        printf("Usage: %s <decoded_data_file> <ifdecoding_flag> <old_found_codeword_file> <address_info> <total_erasure_record_file> <new_found_codeword_file> <result_for_rs_decoding>\n", argv[0]);
        return 1;
    }
    
    char *decoded_filename = argv[1]; // Bit sequence to extract address information (156=12+144)
    char *ifdecoding_flag = argv[2]; // Flag indicating whether the codeword can be decoded (erasure error <= 100)
    char *old_found_codeword_file = argv[3]; // New parameter: found bit sequence file
    char *address_info = argv[4]; // New parameter: address information file
    char *total_erasure_record_file = argv[5]; // New parameter: total erasure record file
    char *new_found_codeword_file = argv[6]; // New parameter: found codeword sequence file
    char *result_for_rs_decoding = argv[7]; // New parameter: result file if erasure < 100
    
    // Create a new array to store voting data read from old_found_codeword_file
    VoteResult **found_votes = (VoteResult **)calloc(MAX_ADDR, sizeof(VoteResult *));
    int *already_found = (int *)calloc(MAX_ADDR, sizeof(int)); // Mark whether sequence is already found
    
    if (found_votes == NULL || already_found == NULL) {
        printf("Memory allocation failed\n");
        return 1;
    }
    
    // Initialize voting array
    for (int i = 0; i < MAX_ADDR; i++) {
        found_votes[i] = NULL;
    }
    
    // Check if old_found_codeword_file exists and read voting data
    if (file_exists(old_found_codeword_file)) 
    {
        // printf("Found existing codeword file %s, loading voting data...\n", old_found_codeword_file);
        FILE *found_file = fopen(old_found_codeword_file, "r");
        if (found_file != NULL) 
        {
            char line[1000];
            int addr;
            
            while (fgets(line, sizeof(line), found_file) != NULL) 
            {
                // Read address line "> position"
                if (line[0] == '>') 
                {
                    sscanf(line + 1, "%d", &addr);
                    
                    if (addr >= 0 && addr < MAX_ADDR) 
                    {
                        // Allocate voting structure for this address
                        if (found_votes[addr] == NULL) 
                        {
                            found_votes[addr] = (VoteResult *)calloc(1, sizeof(VoteResult));
                            if (found_votes[addr] == NULL) 
                            {
                                printf("Memory allocation failed for address %d\n", addr);
                                continue;
                            }
                        }
                        already_found[addr] = 1;
                        
                        // Read voting data lines
                        for (int pos = 0; pos < RESULT_LENGTH; pos++) 
                        {
                            if (fgets(line, sizeof(line), found_file) != NULL) 
                            {
                                int vote0, vote1;
                                sscanf(line, "%d %d", &vote0, &vote1);
                                found_votes[addr]->vote_0[pos] = vote0;
                                found_votes[addr]->vote_1[pos] = vote1;
                                found_votes[addr]->total_votes = vote0 + vote1;
                            }
                        }
                    }
                }
            }
            
            fclose(found_file);
            // printf("Loaded existing voting data\n");
        }
    }
    char a[200]={0};
    char bit_stream[BIT_LENGTH + 2];
    int num1, num2;
    double error_rate;
    int processed = 0;
    int ture_address_num = 0;
    int error_address_num = 0;
    // Read file and process data
    int vaild_read_num = 0;
    int error_read_num = 0;

    // Open decoded data file
    FILE *file = fopen(decoded_filename, "r");
    if (file == NULL) 
    {
        printf("Cannot open file %s\n", decoded_filename);
        return 1;
    }
    while (1) 
    {
        char meta_line[100];
        fgets(meta_line, sizeof(meta_line), file);
        // sscanf(meta_line, "%s %d %lf", a, &num2, &error_rate);
        fgets(bit_stream, BIT_LENGTH + 2, file);
        if(feof(file)) 
            break;

        // Extract address
        uint32_t address = bin_to_dec(bit_stream, TARGET_START, TARGET_END);
        processed++;
        
        // Check if the address is out of range
        if (address >= MAX_ADDR)
        {
            continue;
        }
        
        // Allocate voting structure for this address (if not already allocated)
        if (found_votes[address] == NULL) 
        {
            found_votes[address] = (VoteResult *)calloc(1, sizeof(VoteResult));
            if (found_votes[address] == NULL) 
            {
                printf("Failed to allocate memory for address %u\n", address);
                continue;
            }
        }
        
        // Perform voting: count the 144-bit data excluding the address part
        int bit_idx = 0;
        for (int i = 0; i < BIT_LENGTH; i++) 
        {
            if (i < TARGET_START || i > TARGET_END) // Skip address bits
            {
                char bit = bit_stream[i];
                if (bit == '0') {
                    found_votes[address]->vote_0[bit_idx]++;
                } else if (bit == '1') {
                    found_votes[address]->vote_1[bit_idx]++;
                }
                bit_idx++;
            }
        }
        
        // Update total vote count
        found_votes[address]->total_votes++;
        already_found[address] = 1;
        // if (processed % 10000 == 0) 
        // {
        //     printf("Processed %d records\n", processed);
        // }
    }
    // printf("Processed %d records\n", processed);
    fclose(file);
    // printf("Valid read ratio: %lf\n",(double) vaild_read_num/processed);  
    FILE *found_file = fopen(new_found_codeword_file, "w");
    if (found_file == NULL) {
        printf("Cannot create file %s\n", new_found_codeword_file);
        goto cleanup;
    }
    
    for (uint32_t addr = 0; addr < MAX_ADDR; addr++) 
    {
        if (already_found[addr] && found_votes[addr] != NULL) 
        {
            fprintf(found_file, ">%u\n", addr); // Output address
            // Output voting data for each bit position
            for (int pos = 0; pos < RESULT_LENGTH; pos++) 
            {
                fprintf(found_file, "%d %d\n", found_votes[addr]->vote_0[pos], found_votes[addr]->vote_1[pos]);
            }
        }
    }
    
    fclose(found_file);
    // printf("Updated voting data saved to %s\n", new_found_codeword_file);

    // Create file to save address and sample count information
    FILE *addr_file = fopen(address_info, "w");
    if (addr_file == NULL) {
        printf("Cannot create address info file\n");
        goto cleanup;
    }
    
    // Record address and sample count
    fprintf(addr_file, "address\tcopies\n");
    for (uint32_t addr = 0; addr < MAX_ADDR; addr++) 
    {
        // Only record if found
        int found = ( already_found[addr] ) ? 1 : 0;
        fprintf(addr_file, "%u\t%d\n", addr, found);
    }
    fclose(addr_file);
    // printf("Address and sample count information saved to address_info.txt\n");
    
    // Count total erasure sequences (addresses not found)
    int total_erasure_seq = 0;
    for (uint32_t addr = 0; addr < MAX_ADDR; addr++) {
        if (!already_found[addr]) {
            total_erasure_seq++;
        }
    }

   
    printf("Total erased sequences: %d (%.3f)\n", total_erasure_seq, (double)total_erasure_seq/MAX_ADDR);
    

    // Create a 2D array to store voting results (MAX_ADDR rows, 144 columns)
    char **result_matrix = (char **)malloc(MAX_ADDR * sizeof(char *));
    if (result_matrix == NULL) {
        printf("Memory allocation failed for result matrix\n");
        return 1;
    }
    for (int i = 0; i < MAX_ADDR; i++) {
        result_matrix[i] = (char *)malloc(RESULT_LENGTH * sizeof(char));
        if (result_matrix[i] == NULL) {
            printf("Memory allocation failed for result matrix row %d\n", i);
            return 1;
        }
        // Initialize with 'E' (erasure)
        for (int j = 0; j < RESULT_LENGTH; j++) {
            result_matrix[i][j] = 'E';
        }
    }
    
    // Fill the result matrix with voting results
    for (uint32_t addr = 0; addr < MAX_ADDR; addr++) 
    {
        if (already_found[addr] && found_votes[addr] != NULL) 
        {
            // Generate 144-bit sequence based on voting results
            for (int pos = 0; pos < RESULT_LENGTH; pos++) 
            {
                // Majority voting decision
                if (found_votes[addr]->vote_1[pos] > found_votes[addr]->vote_0[pos]) {
                    result_matrix[addr][pos] = '1';
                } else if (found_votes[addr]->vote_0[pos] > found_votes[addr]->vote_1[pos]) {
                    result_matrix[addr][pos] = '0';
                } else {
                    result_matrix[addr][pos] = 'N';
                }
            }
        }
        // If not found, keep 'E' (already initialized)
    }
    
    // Calculate erasure errors for each RS codeword (12 groups, each with 12 columns)
    int max_erasure_errors = 0;
    int group_erasure_errors[12];  // Array to store erasure errors for each group
    int erasure_symbols = 0;

    for (int group = 0; group < 12; group++) 
    {
        int start_col = group * 12;
        int end_col = start_col + 12;
        
        // Count erasure errors in this group
        group_erasure_errors[group] = 0;
        for (uint32_t addr = 0; addr < MAX_ADDR; addr++) 
        {
            int has_erasure = 0;
            // Check if any of the 12 columns in this group has 'E' or 'N'
            for (int col = start_col; col < end_col; col++) 
            {
                if (result_matrix[addr][col] == 'E' || result_matrix[addr][col] == 'N') 
                {
                    has_erasure = 1;
                    break;
                }
            }
            if (has_erasure) {
                erasure_symbols++;
                group_erasure_errors[group]++;
            }
        }
        
        // Update maximum erasure errors
        if (group_erasure_errors[group] > max_erasure_errors) {
            max_erasure_errors = group_erasure_errors[group];
        }
    }

    // Output all 12 erasure errors in one line
    printf("Erasure errors of 12 RS codewords: ");
    for (int i = 0; i < 12; i++) {
        printf("%d ", group_erasure_errors[i]);
    }
    printf("\n");

    printf("Maximum erasure errors among all RS codewords: %d\n", max_erasure_errors);


    FILE *fw_total_erasure = fopen(total_erasure_record_file, "w");
    if (fw_total_erasure == NULL) 
    {
        printf("Cannot create file %s\n", total_erasure_record_file);
        return 1;
    }
    for (int i = 0; i < 12; i++) 
    {
        fprintf(fw_total_erasure,"%d ", group_erasure_errors[i]);
    }
    fclose(fw_total_erasure);


    FILE *fw_flag = fopen(ifdecoding_flag, "w");
    if (fw_flag == NULL) {
        printf("Cannot create file %s\n", ifdecoding_flag);
        // Free result matrix before return
        for (int i = 0; i < MAX_ADDR; i++) {
            free(result_matrix[i]);
        }
        free(result_matrix);
        return 1;
    }
   
    // Save results based on maximum erasure errors among all RS codewords
    if (max_erasure_errors < 95) 
    {
        fprintf(fw_flag, "yes\n");
        // printf("Maximum erasure errors (%d) < 100. Output: yes\n", max_erasure_errors);
    } 
    else 
    {
        fprintf(fw_flag, "no\n");
        // printf("Maximum erasure errors (%d) >= 100. Output: no\n", max_erasure_errors);
    }
    fclose(fw_flag);

    FILE *result_file = fopen(result_for_rs_decoding, "w");
    if (result_file == NULL) 
    {
        printf("Cannot create result file %s\n", result_for_rs_decoding);
        goto cleanup;
    }
    
    // Write decision results by address order using the result matrix
    for (uint32_t addr = 0; addr < MAX_ADDR; addr++) 
    {
        for (int pos = 0; pos < RESULT_LENGTH; pos++) 
        {
            fputc(result_matrix[addr][pos], result_file);
        }
        fputc('\n', result_file);
    }
    
    // Close result file
    fclose(result_file);
    

cleanup:
    // Free result matrix
    for (int i = 0; i < MAX_ADDR; i++) {
        free(result_matrix[i]);
    }
    free(result_matrix);
    
    // Free memory
    for (int i = 0; i < MAX_ADDR; i++) {
        if (found_votes != NULL && found_votes[i] != NULL)
        {
            free(found_votes[i]);
        }
    }
   
    if (found_votes != NULL) free(found_votes);
    if (already_found != NULL) free(already_found);
    
    return 0;
}



char dexor(char a, char b) // XOR operation
{
    //char output;
    if (a == '0' && b == '0')
    {
        // output='0';
        return '0';
    }
    if (a == '0' && b == '1')
    {
        // output='1';
        return '1';
    }
    if (a == '1' && b == '0')
    {
        // output='1';
        return '1';
    }
    if (a == '1' && b == '1')
    {
        // output='0';
        return '0';
    }
    return '0';
}

char mapping(char a, char b) // Mapping
{
    //char output;
    if (a == '0' && b == '0')
    {
        // output='0';
        return 'A';
    }
    if (a == '0' && b == '1')
    {
        // output='1';
        return 'T';
    }
    if (a == '1' && b == '0')
    {
        // output='1';
        return 'G';
    }
    if (a == '1' && b == '1')
    {
        // output='0';
        return 'C';
    }
    return 'N';
}

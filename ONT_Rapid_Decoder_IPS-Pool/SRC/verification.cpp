
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAXLEN  561600

int array[MAXLEN];
char buff[MAXLEN];

char write2[MAXLEN / 8];
int writepicture[MAXLEN / 8] = {0};

void fileReconstruct(const char* filename,const char* inputStream,int pictureSize)
{
	int i, j, t1 = 0;
	int len1 = 0;

	FILE* fp_input = NULL;

	FILE* fw_file = NULL;

	fp_input = fopen(inputStream, "rt"); 	
	if (fp_input == NULL)
	{
		printf("cannot open the file!\n");
		exit(0);
	}
	while (1)
	{
		int a = fscanf(fp_input, "%s", buff);
		if (feof(fp_input))
			break;
		printf("sourceBIT : %d bits\n", a);
	}


	fw_file = fopen(filename, "wb");
	if (fw_file == NULL)
	{
		printf("cannot open the png\n");
		exit(0);

	}

	len1 = strlen(buff);


	for (i = 0; i < len1; i++)
	{
		if (buff[i] == '0')
			array[i] = 0;
		else
			array[i] = 1;
	}

	t1 = 0;
	

	for (i = 0; i < (len1 / 8); i++)
	{
		
		writepicture[i] = 0;
		
	}
	for (i = 0; i < (len1 / 8); i++)
	{
		for (j = 8 * i; j < 8 * (i + 1); j++)
		{
			writepicture[t1] += (int)(array[j] * pow(2, (7 - (j - 8 * i))));
		}
		t1++;
	}

	for (i = 0; i < pictureSize; i++)
	{
		write2[i] = (char)writepicture[i];      
	}

	printf("%d\n", t1);

	fwrite(write2, sizeof(char), pictureSize, fw_file);

	fclose(fp_input);

	fclose(fw_file);
	
}
int main(int argc, char* argv[])
{
	
	// Check if command line arguments are provided
	if (argc != 4) {
		printf("Usage: %s <input_file> <output_file> <picture_size>\n", argv[0]);
		printf("Example: %s stream_1.txt output.png 69602\n", argv[0]);
		return 1;
	}

	char* inputStream = argv[1];
	char* filename = argv[2];
	int pictureSize = atoi(argv[3]);

	fileReconstruct(filename, inputStream, pictureSize);

	return 0;
}

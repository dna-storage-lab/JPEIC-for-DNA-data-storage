#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <float.h>
#include <limits.h>

#define MAX_RANDOM LONG_MAX    // Maximum value of random() 
#define informLen 561600
#define rsCodeNum 12
#define rsCodeLen 12
#define rs_k 3900
#define rs_n 4000
#define payloadLen 130


int ri;
int t;
int m, n, length, k, t2, d, red;
int init_zero;
int* p;
int numerr, errpos[600], errval[600], decerror;
int biterror, error;
char filename[40], name2[40];
int numera = 0;
int era[600]={0}, eraval[600]={0};
int* g, * rdata, * index_of, * alpha_to, * recd, * b;
	
char buff_scramble[561600] = {0};

void read_p(void);
void generate_gf(void);
void gen_poly(void);
void encode_rs(void);
//void bpsk_awgn(void);
void decode_rs(void);
//int weight(int word);

using namespace std;
char dexor1(int a, char b)
{
	//char output;
	if (a == 0 && b == '0')
	{
		// output='0';
		return 0;
	}
	if (a == 0 && b == '1')
	{
		// output='1';
		return 1;
	}
	if (a == 1 && b == '0')
	{
		// output='1';
		return 1;
	}
	if (a == 1 && b == '1')
	{
		// output='0';
		return 0;
	}
	return 0;
}
int main(int args, char** argv)
{
	if(args != 5)
	{
		printf("Usage: %s <LDPC_decoded_file> <bit_output_file> <success_decode_file>\n", argv[0]);
		return 1;
	}

	// printf("LDPC decoded file: %s\n", argv[1]);
    // printf("bit output file : %s\n", argv[2]);
	// printf("scramble sequence: %s\n", argv[3]);
	// printf("success decode file: %s\n", argv[4]);
    char *LDPC_decoded_file = argv[1];
    char *bit_output_file = argv[2];
	char *scramble_sequence = argv[3];
	char *success_decode_file = argv[4];

	FILE* fpQuery = NULL;
	fpQuery = fopen(LDPC_decoded_file, "rt"); //140bit
	if (!fpQuery)
	{
		printf("cannot open the sequence.txt file!\n");
		return 0;
	}

	FILE* fwOutstream = NULL;
	fwOutstream = fopen(bit_output_file, "wt");

	if (!fwOutstream)
	{
		printf("cannot open the out_txt file!\n");
		return 0;
	}

	FILE* fp_success_decode = NULL;
	fp_success_decode = fopen(success_decode_file, "wt");
	if (!fp_success_decode)
	{
		printf("cannot open the success decode file!\n");
		return 0;
	}


	int i = 0, j = 0;
	int pos_rs = 0;
	char buff[500] = { 0 };
	int reads_num = 0;
	int now_index = 0;
	//rs erase
	int ERROR[rsCodeNum][rs_n] = { 0 };
	
	for (int x = 0; x < 600; x++)
	{
		era[x] = 0;
	}
	numera = 0;

	int** rsCodeWithErr;
	rsCodeWithErr = (int**)malloc(sizeof(int*) * rsCodeNum);
	for (int i1 = 0; i1 < rsCodeNum; i1++)
	{
		rsCodeWithErr[i1] = (int*)malloc(sizeof(int) * rs_n);
		memset(rsCodeWithErr[i1], 0, sizeof(int) * rs_n);
	}

	int** codewordRS;
	codewordRS = (int**)malloc(sizeof(int*) * rs_k);
	for (int r = 0; r < rs_k; r++)
	{
		codewordRS[r] = (int*)malloc(sizeof(int) * rsCodeNum);
	}
	

	int** sequence11;
	sequence11 = (int**)malloc(sizeof(int*) * 4000);
	for (int l = 0; l < 4000; l++)
	{
		sequence11[l] = (int*)malloc(sizeof(int) * 144);
	}

	for (int x = 0; x < 4000; x++)
	{
		for (int y = 0; y < 144; y++)
		{
			sequence11[x][y] = 0;
		}
	}


	char** seq_after_ldpc;
	seq_after_ldpc = (char**)malloc(sizeof(char*) * rs_n);
	for (int l = 0; l < rs_n; l++)
	{
		seq_after_ldpc[l] = (char*)malloc(sizeof(char) * (rsCodeLen * rsCodeNum));
	}

	int** bit_after_ldpc;
	bit_after_ldpc = (int**)malloc(sizeof(int*) * rs_n);
	for (int l = 0; l < rs_n; l++)
	{
		bit_after_ldpc[l] = (int*)malloc(sizeof(int) * (rsCodeLen * rsCodeNum));
	}

	//=======================================================================//
	g = (int*)malloc(600 * sizeof(int));
	rdata = (int*)malloc(4096 * sizeof(int));
	index_of = (int*)malloc(4096 * sizeof(int));
	alpha_to = (int*)malloc(4096 * sizeof(int));
	recd = (int*)malloc(4096 * sizeof(int));
	b = (int*)malloc(600 * sizeof(int));
	p = (int*)malloc(600 * sizeof(int));

	m = rsCodeLen;
	length = 4095;
	red = 100;
	init_zero = 0;
	k = length - red;//3995
	t = red / 2;
	t2 = 2 * t;

	read_p();        // Read m */

	generate_gf();   // Construct the Galois Field GF(2^m) */

	gen_poly();      // Compute the generator polynomial of RS code */
	
	
	/*------------------------------------------------------------------------------------*/

	char** original;
	original = (char**)malloc(sizeof(char*) * rs_k);
	for (int r = 0; r < rs_k; r++)
	{
		original[r] = (char*)malloc(sizeof(char) * 141);
	}

	int seq_num = 0;
	char readWithoutPrimer[150] = { 0 };

	
	// Create a 1D array to store all sequences if it doesn't exist yet
	int* re_k = NULL;
	re_k = (int*)malloc(sizeof(int) * rs_k * 144);
	

	/*------------------------------------------------------*/
	reads_num = 0;

	while (reads_num < rs_n)
	{
		fscanf(fpQuery, "%s\n", readWithoutPrimer);
		
		reads_num++;

		for (int j = 0; j < 144; j++)
		{

			seq_after_ldpc[now_index][j] = readWithoutPrimer[j];

		}
		
		now_index++;

	}

	// printf("Finish reading %d sequences with 144bits\n",reads_num);

	for (i = 0; i < rs_n; i++)
	{
		
		int sum = 0;
		for (int n1 = 0; n1 < rsCodeNum; n1++)
		{
			for (int gg = 0; gg < rsCodeLen; gg++)
			{
				if(seq_after_ldpc[i][n1*rsCodeLen+gg] == 'E' || seq_after_ldpc[i][n1*rsCodeLen+gg] == 'N')
				{
					bit_after_ldpc[i][n1*rsCodeLen+gg] = 0;
					ERROR[n1][i] = 1;
				}
				else if(seq_after_ldpc[i][n1*rsCodeLen+gg] == '0')
				{
					bit_after_ldpc[i][n1*rsCodeLen+gg] = 0;
				}
				else if(seq_after_ldpc[i][n1*rsCodeLen+gg] == '1')
				{
					bit_after_ldpc[i][n1*rsCodeLen+gg] = 1;
				}	
			}
			
		}
		
		int sym_sum12 = 0;
		for (int g2 = 0; g2 < rsCodeNum; g2++)
		{
			sym_sum12 = 0;
			int base_index = 12 * g2;
			for (int gg = 0; gg < rsCodeLen; gg++)
			{
				sym_sum12 = (sym_sum12 << 1) + bit_after_ldpc[i][base_index + gg];
			}
			rsCodeWithErr[g2][i] = sym_sum12;
		}

	}
	
	numera = 0;
	int errorsymbol=0;
	
	bool is_success_decode = true;
	int record_flag=0;
	int sum_error = 0, sum_error1 = 0;
	for (i = 0; i < rsCodeNum; i++)
	{
		for (int x = 0; x < 600; x++)
		{
			era[x] = 0;
		}
		numera = 0;
		errorsymbol = 0;

		//erasure position

		for (j = 0; j < rs_n; j++)
		{
			
			if (ERROR[i][j] == 1)
			{
				era[numera] = 95+j;
				numera++;

			}
			
		}
		
		int error_num = 0;
		for (j = 0; j < 95; j++)
		{
			recd[j] = 0;
		}
		for (j = 0; j < rs_k; j++)
		{
			recd[j + 95] = rsCodeWithErr[i][j];
		}

		for (j = 3995; j < 4095; j++)
		{
			recd[j] = rsCodeWithErr[i][j - 95];
		}
		
		decode_rs();
		int CodeCheck[4096] = { 0 };
		error_num = 0;
		for (j = 0; j < rs_k; j++)
		{
			codewordRS[j][i] = recd[j + 95];
		}
		for (j = 0; j < 4095; j++)
		{
			CodeCheck[j] = *(recd + j);
			//printf("%d ", code[i1+251]);
		}
		//rs check after decoding
		int syn_error = 0;
		int s[200] = { 0 };
		for (int mm = 1; mm <= t2; mm++)
		{
			s[mm] = 0;
			for (int kk = 0; kk < length; kk++)
				if (CodeCheck[kk] != 0)
					s[mm] ^= alpha_to[(index_of[CodeCheck[kk]] + (mm + init_zero - 1) * kk) % n];
			// convert syndrome from vector form to log form  */
			if (s[mm] != 0)
				syn_error = 1;         // set flag if non-zero syndrome => error
				//
				// Note:    If the code is used only for ERROR DETECTION, then
				//          exit program here indicating the presence of errors.
				//
			s[mm] = index_of[s[mm]];
		}
		if(syn_error==1 && !record_flag)
		{
			is_success_decode = false;
			fprintf(fp_success_decode, "no");
			fclose(fp_success_decode);
			record_flag=1;
			// return 0;
		}

	}
	
	if(is_success_decode)
	{
		fprintf(fp_success_decode, "yes");
		fclose(fp_success_decode);
	}
	
	for (int jj = 0; jj < 3900; jj++)//
	{
		int sequence0[144] = { 0 };

		for (int aa = 0; aa < 12; aa++)
		{
			
			int value_rs = codewordRS[jj][aa];
			pos_rs = 11;
			while (value_rs > 0 && pos_rs >= 0)
			{
				sequence0[12 * aa + pos_rs] = value_rs % 2;
				value_rs = value_rs / 2;
				pos_rs--;
			}
		}

		for (i = 0; i < 144; i++)
		{
			sequence11[jj][i] = sequence0[i];
			
		}

		/*outfile<<endl;*/
	}
	// printf("===================RS decoding completed=====================\n");
	
	FILE* fp_scramble=NULL;

	fp_scramble = fopen(scramble_sequence, "r");
	if (fp_scramble == NULL)
	{
		printf("can not open scramble_sequence\n");
		exit(0);
	}
	int ii = 0;
	int len = 0;
	len = fread(buff_scramble, sizeof(char), 561600, fp_scramble);
	fclose(fp_scramble);
	// printf("scrambling code : %d bits\n", len);
	char write_buffer[3900 * 144 + 1]; 
	int buffer_index = 0;

	for (j = 0; j < 3900; j++)
	{
		for (i = 0; i < 144; i++)
		{
			sequence11[j][i] = dexor1(sequence11[j][i], buff_scramble[ii]);
			ii++;
			write_buffer[buffer_index++] = sequence11[j][i] + '0'; 
		}
	}
	write_buffer[buffer_index] = '\0';

	fprintf(fwOutstream, "%s", write_buffer);
	
	for (int ii1 = 0; ii1 < rsCodeNum; ii1++)
		free(rsCodeWithErr[ii1]);
	free(rsCodeWithErr);


	for (int rr2 = 0; rr2 < rs_k; rr2++)
		free(codewordRS[rr2]);
	free(codewordRS);


	for (int rr2 = 0; rr2 < rs_n; rr2++)
		free(seq_after_ldpc[rr2]);
	free(seq_after_ldpc);


	for (int ll = 0; ll < rs_k; ll++)
		free(original[ll]);
	free(original);

	// free(re_k);
	// free(rdata);
	// free(g);
	// free(index_of);
	// free(alpha_to);
	// free(recd);
	// free(b);
	// free(p);


	fclose(fwOutstream);
	fclose(fpQuery);
	
	return 0;
}

void read_p()
//      Read m, the degree of a primitive polynomial p(x) used to compute the
//      Galois field GF(2**m). Get precomputed coefficients p[] of p(x). Read
//      the code length.
{
	int i, ninf = 0;
	/*
	printf("\nSimulation of RS codes \n");
	printf("Copyright 2002 (c)Robert Morelos-Zaragoza. All rights reserved.\n\n");*/

	for (i = 1; i < m; i++)
		p[i] = 0;
	p[0] = p[m] = 1;

	if (m == 2)             p[1] = 1;
	else if (m == 3)        p[1] = 1;
	else if (m == 4)        p[3] = 1;
	// else if (m == 4)        p[1] = 1;  // Commented out to match example p. 68
	else if (m == 5)        p[2] = 1;
	else if (m == 6)        p[1] = 1;
	else if (m == 7)        p[1] = 1;
	else if (m == 8)        p[4] = p[5] = p[6] = 1;
	else if (m == 9)        p[4] = 1;
	else if (m == 10)       p[3] = 1;
	else if (m == 11)       p[2] = 1;
	else if (m == 12)       p[3] = p[4] = p[7] = 1;
	else if (m == 13)       p[1] = p[3] = p[4] = 1;
	else if (m == 14)       p[1] = p[11] = p[12] = 1;
	else if (m == 15)       p[1] = 1;
	else if (m == 16)       p[2] = p[3] = p[5] = 1;
	else if (m == 17)       p[3] = 1;
	else if (m == 18)       p[7] = 1;
	else if (m == 19)       p[1] = p[5] = p[6] = 1;
	else if (m == 20)       p[3] = 1;
	/*  printf("Primitive polynomial of GF(2^%d), (LSB first)   p(x) = ",m);*/

	n = 1;
	for (i = 0; i <= m; i++)
	{
		n *= 2;
		//     printf("%1d", p[i]);
	}
	//  printf("\n");

	n = n / 2 - 1;
}



void generate_gf()
// generate GF(2^m) from the irreducible polynomial p(X) in p[0]..p[m]
//
// lookup tables:  log->vector form           alpha_to[] contains j=alpha**i;
//                 vector form -> log form    index_of[j=alpha**i] = i
// alpha=2 is the primitive element of GF(2^m)
{

	register int i, mask;
	mask = 1;
	alpha_to[m] = 0;
	for (i = 0; i < m; i++)
	{
		alpha_to[i] = mask;
		index_of[alpha_to[i]] = i;
		if (p[i] != 0)
			alpha_to[m] ^= mask;
		mask <<= 1;
	}

	index_of[alpha_to[m]] = m;
	mask >>= 1;
	for (i = m + 1; i < n; i++)
	{
		if (alpha_to[i - 1] >= mask)
			alpha_to[i] = alpha_to[m] ^ ((alpha_to[i - 1] ^ mask) << 1);
		else alpha_to[i] = alpha_to[i - 1] << 1;
		index_of[alpha_to[i]] = i;
	}
	index_of[0] = -1;
	//#define PRINT_GF

#ifdef PRINT_GF/*
printf("Table of GF(%d):\n",n+1);
printf("----------------------\n");
printf("   i\tvector \tlog\n");
printf("----------------------\n");
for (i=0; i<=n; i++)
printf("%4d\t%4x\t%4d\n", i, alpha_to[i], index_of[i]);*/
#endif
}


void gen_poly()
// Compute the generator polynomial of the t-error correcting, length
// n=(2^m -1) Reed-Solomon code from the product of (X+alpha^i), for
// i = init_zero, init_zero + 1, ..., init_zero+length-k-1
{
	register int i, j;


	g[0] = alpha_to[init_zero];  //  <--- vector form of alpha^init_zero
	g[1] = 1;     // g(x) = (X+alpha^init_zero)
	for (i = 2; i <= length - k; i++)
	{
		g[i] = 1;
		for (j = i - 1; j > 0; j--)
			if (g[j] != 0)
				g[j] = g[j - 1] ^ alpha_to[(index_of[g[j]] + i + init_zero - 1) % n];
			else
				g[j] = g[j - 1];
		g[0] = alpha_to[(index_of[g[0]] + i + init_zero - 1) % n];
	}

	// convert g[] to log form for quicker encoding 
	for (i = 0; i <= length - k; i++)  g[i] = index_of[g[i]];


	//#define PRINT_POLY
#ifdef PRINT_POLY
	//printf("Generator polynomial (independent term first):\ng(x) = ");
	/*for (i=0; i<=length-k; i++) printf("%5d", g[i]);
	printf("\n");*/
#endif
}


void encode_rs()
// Compute the 2t parity symbols in b[0]..b[2*t-1]
// data[] is input and b[] is output in polynomial form.
// Encoding is done by using a feedback shift register with connections
// specified by the elements of g[].

{
	register int i, j;
	int feedback;


	for (i = 0; i < length - k; i++)
		b[i] = 0;
	for (i = k - 1; i >= 0; i--)
	{
		feedback = index_of[rdata[i] ^ b[length - k - 1]];
		if (feedback != -1)
		{
			for (j = length - k - 1; j > 0; j--)
				if (g[j] != -1)
					b[j] = b[j - 1] ^ alpha_to[(g[j] + feedback) % n];
				else
					b[j] = b[j - 1];
			b[0] = alpha_to[(g[0] + feedback) % n];
		}
		else
		{
			for (j = length - k - 1; j > 0; j--)
				b[j] = b[j - 1];
			b[0] = 0;
		}
	}
}


void decode_rs()
{
	register int i, j, u, q;
	int count = 0, syn_error = 0;
	int degphi, ell, temp = 0;


	int* d;
	d = (int*)malloc(600 * sizeof(int));
	int* l;
	l = (int*)malloc(600 * sizeof(int));
	int* u_lu;
	u_lu = (int*)malloc(600 * sizeof(int));
	int* s, * forney, * err, * omega, * phi, * phiprime;
	s = (int*)malloc(600 * sizeof(int));
	forney = (int*)malloc(600 * sizeof(int));
	err = (int*)malloc(4100 * sizeof(int));
	omega = (int*)malloc(4100 * sizeof(int));
	phi = (int*)malloc(4100 * sizeof(int));
	phiprime = (int*)malloc(4100 * sizeof(int));
	int* tau, * root, * loc, * z, * reg, * aux;
	tau = (int*)malloc(600 * sizeof(int));
	root = (int*)malloc(600 * sizeof(int));
	loc = (int*)malloc(600 * sizeof(int));
	z = (int*)malloc(600 * sizeof(int));
	reg = (int*)malloc(600 * sizeof(int));
	aux = (int*)malloc(600 * sizeof(int));
	/*
	int **elp = new int*[10536];
		for (i = 0; i < 10536; i++)
			elp[i] = new int[65538];

	*/


	int** elp;
	elp = (int**)malloc(sizeof(int*) * 600);
	for (i = 0; i < 600; i++)
		elp[i] = (int*)malloc(sizeof(int) * 600);

	// Compute the syndromes

 //#ifdef PRINT_SYNDROME
 //   printf("\ns =         0 ");
 //#endif
	for (i = 1; i <= t2; i++)
	{
		s[i] = 0;
		for (j = 0; j < length; j++)
			if (recd[j] != 0)
				s[i] ^= alpha_to[(index_of[recd[j]] + (i + init_zero - 1) * j) % n];
		// convert syndrome from vector form to log form  */
		if (s[i] != 0)
			syn_error = 1;         // set flag if non-zero syndrome => error
		  //
		  // Note:    If the code is used only for ERROR DETECTION, then
		  //          exit program here indicating the presence of errors.
		  //
		s[i] = index_of[s[i]];
		//#define PRINT_SYNDROME
#ifdef PRINT_SYNDROME
//   printf("%4d ", s[i]);
#endif
	}
	//printf("\nsrr =         0 ");
	// printf("INFO: check result%4d\n", syn_error);
	if (syn_error)       // if syndromes are nonzero then try to correct
	{

		s[0] = 0;  // S(x) = 1 + s_1x + ...

		 // TO HANDLE ERASURES

		if (numera)
			// if erasures are present, compute the erasure locator polynomial, tau(x)
		{
			for (i = 0; i <= t2; i++)
			{
				tau[i] = 0; aux[i] = 0;
			}

			aux[1] = alpha_to[era[0]];
			aux[0] = 1;       // (X + era[0]) 

			if (numera > 1)
				for (i = 0; i < numera; i++)
				{
					p[1] = era[i];
					p[0] = 0;
					for (j = 0; j < 2; j++)
						for (ell = 0; ell <= i; ell++)
							// Line below added 8/17/2003
							if ((p[j] != -1) && (aux[ell] != 0))
								tau[j + ell] ^= alpha_to[(p[j] + index_of[aux[ell]]) % n]; 
								// printf("%4d  %4d\n", ell, tau[ell]);
					if (i != (numera - 1))
						for (ell = 0; ell <= (i + 1); ell++)
						{
							aux[ell] = tau[ell];
							tau[ell] = 0;
						}
				}

			else {
				tau[0] = aux[0]; tau[1] = aux[1];
			}


			// Put in index (log) form
			for (i = 0; i <= numera; i++)
				tau[i] = index_of[tau[i]]; /* tau in log form */


#ifdef PRINT_SYNDROME
//printf("\ntau =    ");
//for (i=0; i<=numera; i++)
 // printf("%4d ", tau[i]);
//printf("\nforney = ");
#endif


		// Compute FORNEY modified syndrome:
		//            forney(x) = [ s(x) tau(x) + 1 ] mod x^{t2}

			for (i = 0; i <= n - k; i++) forney[i] = 0;
			for (i = 0; i <= n - k; i++)
				for (j = 0; j <= numera; j++)
					if (i + j <= (n - k)) // mod x^{n-k+1}
						if ((s[i] != -1) && (tau[j] != -1))
							forney[i + j] ^= alpha_to[(s[i] + tau[j]) % n];

			forney[0] ^= 1;

			for (i = 0; i <= n - k; i++)
				forney[i] = index_of[forney[i]];

#ifdef PRINT_SYNDROME
			/*for (i=0; i<=n-k; i++)
			  printf("%4d ", forney[i]);*/
#endif

		}

		else // No erasures
		{
			tau[0] = 0;
			for (i = 1; i <= n - k; i++) forney[i] = s[i];
		}

#ifdef PRINT_SYNDROME
		//printf("\n");
#endif

  // --------------------------------------------------------------
  //    THE BERLEKAMP-MASSEY ALGORITHM FOR ERRORS AND ERASURES
  // --------------------------------------------------------------

	  // initialize table entries
		d[0] = 0;                // log form
		d[1] = forney[numera + 1]; // log form
		elp[0][0] = 0;           // log form 
		elp[1][0] = 1;           // vector form 
		for (i = 1; i < t2; i++)
		{
			elp[0][i] = -1;   // log form
			elp[1][i] = 0;    // vector form 
		}
		l[0] = 0;
		l[1] = 0;
		u_lu[0] = -1;
		u_lu[1] = 0;
		u = 0;

		// if (numera < t2) 
		{  // If errors can be corrected
			do
			{
				u++;

				if (d[u] == -1)
				{
#ifdef PRINT_SYNDROME
					//printf("d[%d] is zero\n",u);
#endif


					l[u + 1] = l[u];

					for (i = 0; i <= l[u]; i++)
					{
						//	printf("elp[u][i]=%d,index_of[elp[u][i]]=%d",elp[u][i],index_of[elp[u][i]]);


						elp[u + 1][i] = elp[u][i];

						elp[u][i] = index_of[elp[u][i]];

					}

				}
				else
					// search for words with greatest u_lu[q] for which d[q]!=0
				{
					q = u - 1;
					while ((d[q] == -1) && (q > 0))    q--;
					// have found first non-zero d[q] 
					if (q > 0)
					{
						j = q;
						do
						{
							j--;
							if ((d[j] != -1) && (u_lu[q] < u_lu[j]))
								q = j;
						} while (j > 0);
					}

#ifdef PRINT_SYNDROME
					//printf("u = %4d, q = %4d, d[q] = %4d d[u] = %4d\n", u, q, d[q],d[u]);
#endif

			// have now found q such that d[u]!=0 and u_lu[q] is maximum 
			// store degree of new elp polynomial 
					if (l[u] > l[q] + u - q)  l[u + 1] = l[u];
					else  l[u + 1] = l[q] + u - q;

#ifdef PRINT_SYNDROME
					//printf("l[q] = %4d, l[u] = %4d\n", l[q], l[u]);
#endif

			// compute new elp(x) 
			// for (i=0; i<t2-numera; i++)    elp[u+1][i] = 0; 
					for (i = 0; i < t2; i++)    elp[u + 1][i] = 0;
					for (i = 0; i <= l[q]; i++)
						if (elp[q][i] != -1)
							elp[u + 1][i + u - q] = alpha_to[(d[u] + n - d[q] + elp[q][i]) % n];
					for (i = 0; i <= l[u]; i++)
					{
						elp[u + 1][i] ^= elp[u][i];
						elp[u][i] = index_of[elp[u][i]];
					}

#ifdef PRINT_SYNDROME
					//printf("l[u+1] = %4d, elp[u+1] = ", l[u+1]);
					/*for (i=0;  i<=l[u+1]; i++) printf("%4d ",index_of[elp[u+1][i]]);
					printf("\n");*/
#endif

				}


				u_lu[u + 1] = u - l[u + 1];
				// compute (u+1)th discrepancy 
				if (u < (t2 - numera)) // no discrepancy computed on last iteration 
				// --- if ( u < (l[u+1]+t-1-(numera/2)) ) 
				{
					// if (s[u+1]!=-1)
					if (forney[numera + u + 1] != -1)
						d[u + 1] = alpha_to[forney[numera + u + 1]];
					else
						d[u + 1] = 0;
#ifdef PRINT_SYNDROME
					//printf("discrepancy for u = %d: d[u+1] = %4d\n", u, index_of[d[u+1]]);
#endif

					for (i = 1; i <= l[u + 1]; i++)
						// if ((s[u+1-i]!=-1) && (elp[u+1][i]!=0))
						//   d[u+1] ^= alpha_to[(s[numera+u+1-i]
						if ((forney[numera + u + 1 - i] != -1) && (elp[u + 1][i] != 0)) {
							d[u + 1] ^= alpha_to[(forney[numera + u + 1 - i]
								+ index_of[elp[u + 1][i]]) % n];

#ifdef PRINT_SYNDROME
							//printf("i=%d, forney[%d] = %4d, d[u+1] = %4d\n",i,numera+u+1-i,
											 //  forney[numera+u+1-i],index_of[d[u+1]]);
#endif
						}
					d[u + 1] = index_of[d[u + 1]];     // put d[u+1] into index form 

#ifdef PRINT_SYNDROME
//printf("d[u+1] = %4d\n", d[u+1]);
#endif
				}
				//	printf("u=%d,(t2-numera)=%d,l[u+1]=%d,(t2-numera)/2)=%d\n",u,t2-numera,l[u+1],(t2-numera)/2);
			} while ((u < (t2 - numera)) && (l[u + 1] <= ((t2 - numera) / 2)));


		}

		// else // case of 2t erasures
		// {
		//   elp[1][0] = 0;
		//   count = 0;
		// }

		u++;
		//     printf("%d %d %d %d",l[u],t-numera/2,t,numera);
		//  if (l[u] <= t - numera / 2)         // can correct errors
		{

			// put elp into index form 
			for (i = 0; i <= l[u]; i++)   elp[u][i] = index_of[elp[u][i]];
			/*
			printf("\nBM algorithm, after %d iterations:\nsigma = ", (u-1));
			for (i=0; i<=l[u]; i++)   printf("%4d ", elp[u][i]);
			printf("\n");*/

			// find roots of the error location polynomial 
			for (i = 1; i <= l[u]; i++)
				reg[i] = elp[u][i];
			count = 0;
			for (i = 1; i <= n; i++)
			{
				q = 1;
				for (j = 1; j <= l[u]; j++)
					if (reg[j] != -1)
					{
						reg[j] = (reg[j] + j) % n;
						q ^= alpha_to[reg[j]];
					}
				if (!q)        // store root and error location number indices 
				{
					root[count] = i;
					loc[count] = n - i;

#ifdef PRINT_SYNDROME
					//printf("loc[%4d] = %4d\n", count, loc[count]);
#endif

					count++;
				}
			}

			if (count == l[u])    // no. roots = degree of elp hence <= t errors 
			{

				// Compute the errata evaluator polynomial, omega(x)

				forney[0] = 0;  // as a log, to construct 1+T(x)
				for (i = 0; i <= t2; i++)
					omega[i] = 0;
				for (i = 0; i <= t2; i++)
				{
					for (j = 0; j <= l[u]; j++)\
					{
						if (i + j <= t2) // mod x^{t2}
							if ((forney[i] != -1) && (elp[u][j] != -1))
								omega[i + j] ^= alpha_to[(forney[i] + elp[u][j]) % n];
					}
				}
				for (i = 0; i <= t2; i++)
					omega[i] = index_of[omega[i]];

#ifdef PRINT_SYNDROME
				/*printf("\nomega =    ");
				for (i=0; i<=t2; i++)   printf("%4d ", omega[i]);
				printf("\n");*/
#endif

				// Compute the errata locator polynomial, phi(x)

				degphi = numera + l[u];
				for (i = 0; i <= degphi; i++) phi[i] = 0;
				for (i = 0; i <= numera; i++)
					for (j = 0; j <= l[u]; j++)
						if ((tau[i] != -1) && (elp[u][j] != -1))
							phi[i + j] ^= alpha_to[(tau[i] + elp[u][j]) % n];
				for (i = 0; i <= degphi; i++)
					phi[i] = index_of[phi[i]];


#ifdef PRINT_SYNDROME
				/*printf("phi =      ");
				for (i=0; i<=degphi; i++)   printf("%4d ", phi[i]);
				printf("\n");*/
#endif

				// Compute the "derivative" of phi(x): phiprime

				for (i = 0; i <= degphi; i++) phiprime[i] = -1; // as a log
				for (i = 0; i <= degphi; i++)
					if (i % 2)  // Odd powers of phi(x) give terms in phiprime(x)
						phiprime[i - 1] = phi[i];

#ifdef PRINT_SYNDROME
				/*printf("phiprime = ");
				for (i=0; i<=degphi; i++)   printf("%4d ", phiprime[i]);
				printf("\n\n");*/
#endif

				//printf("**** count = %d\n", count);

				if (numera)
					// Add erasure positions to error locations
					for (i = 0; i < numera; i++) {
						loc[count + i] = era[i];
						root[count + i] = (n - era[i]) % n;
					}

				// evaluate errors at locations given by errata locations, loc[i]
				//                             for (i=0; i<l[u]; i++)    
				for (i = 0; i < degphi; i++)
				{
					// compute numerator of error term  
					err[loc[i]] = 0;
					for (j = 0; j <= t2; j++)
						if ((omega[j] != -1) && (root[i] != -1))
							err[loc[i]] ^= alpha_to[(omega[j] + j * root[i]) % n];

					// -------  The term loc[i]^{2-init_zero}
					if ((err[loc[i]] != 0) && (loc[i] != -1))
						err[loc[i]] = alpha_to[(index_of[err[loc[i]]]
							+ loc[i] * (2 - init_zero + n)) % n];
					if (err[loc[i]] != 0)
					{
						err[loc[i]] = index_of[err[loc[i]]];
						// compute denominator of error term 
						q = 0;
						for (j = 0; j <= degphi; j++)
							if ((phiprime[j] != -1) && (root[i] != -1))
								q ^= alpha_to[(phiprime[j] + j * root[i]) % n];

						// Division by q
						err[loc[i]] = alpha_to[(err[loc[i]] - index_of[q] + n) % n];

#ifdef PRINT_SYNDROME
						//printf("errata[%4d] = %4d (%4d) \n",loc[i],index_of[err[loc[i]]],err[loc[i]]);
#endif

						recd[loc[i]] ^= err[loc[i]];
					}
				}

				//    printf("\nCor =");
		  //         for (i=0; i<length; i++) {
		  //             printf("%4d ", index_of[recd[i]]);
		  //            }
		  //          printf("\n     ");
			 /*       for (i=0; i<length; i++) {
					   printf("%4d ", recd[i]);
					   }
					printf("\n");*/


			}

			else
				// no. roots != degree of elp => >t errors and cannot solve
				;

		}

		// else         // elp has degree has degree >t hence cannot solve 
		// 	;

	}

	else       // no non-zero syndromes => no errors: output received codeword 
		;

	for (i = 0; i < 600; i++)
		free(elp[i]);
	free(elp);
	free(d);
	free(l);
	free(u_lu);
	free(s);
	free(forney);
	free(err);
	free(omega);
	free(phi);
	free(phiprime);
	free(tau);
	free(root);
	free(loc);
	free(z);
	free(reg);
	free(aux);
	/*
	   for (i=0;i<10536;i++)
	   {
		delete []elp[i];
	   }
	  delete []elp;
	*/
}


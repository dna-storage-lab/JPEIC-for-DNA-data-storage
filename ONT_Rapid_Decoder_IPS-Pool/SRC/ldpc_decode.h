/* RCODE.H - Parity chk and gen matrix storage, and procedures to read them.*/

/* Copyright (c) 2000 by Radford M. Neal 
 *
 * Permission is granted for anyone to copy, use, or modify this program 
 * for purposes of research or education, provided this copyright notice 
 * is retained, and note is made of any changes that have been made. 
 *
 * This program is distributed without any warranty, express or implied.
 * As this program was written for research purposes only, it has not been
 * tested to the degree that would be advisable in any important application.
 * All use of this program is entirely at the user's own risk.
 */


#ifndef LDPC_DECODE_H
#define LDPC_DECODE_H

// #include "mod2sparse.h"
// #include "mod2dense.h"

class LDPCDecode {
private:
    /* VARIABLES HOLDING DATA */
    mod2sparse *H;          /* Parity check matrix */
    int M;                  /* Number of rows in parity check matrix */
    int N;                  /* Number of columns in parity check matrix */
    
    char type;              /* Type of generator matrix representation */
    int *cols;              /* Ordering of columns in generator matrix */
    
    mod2sparse *L, *U;      /* Sparse LU decomposition, if type=='s' */
    int *rows;              /* Ordering of rows in generator matrix (type 's') */
    
    mod2dense *G;           /* Dense or mixed representation of generator matrix */
    
    mod2sparse *Hm;
    int mM;
    int mN;
    
    mod2sparse *invHc;
    int cM;
    int cN;
    
    /* PRIVATE HELPER METHODS */
    void fftgf64(int* index, double* din, double* dout);
    void encodecom(int* indexg, int* src, int* code);
    void de2bi(int w, int din, char *dout);
    int bi2de(int w, int *din);
    void cnu(double *p, double *pifftdi, int *index, int *permutation, int *habc);
    void cnu3(double *p, double *pifftdi, int *index, int *permutation, int *habc);
    void vnu(double *c2v, double *prior, double *posterior, char *dec, double *v2c);
    void readGen(char *gen_file, int cols_only, int no_pchk_file);
    
public:
    /* CONSTRUCTOR AND DESTRUCTOR */
    LDPCDecode();
    ~LDPCDecode();
    
    /* PUBLIC INTERFACE METHODS */
    void readPchk(char *pchk_file);
    void readEncpara(char *pchk_Hm, char *pchk_Hc);
    void encode(mod2sparse *Hm,mod2sparse *invHc,char  *src,char  *code);
    void decode(char *dec, int* success, int* iteration, double *allf1, double *allf0, 
                int iternum, int *index, int *habc, int *checknode);
    
    /* GETTERS */
    mod2sparse* getH() const { return H; }
    mod2sparse* getHm() const { return Hm; }
    mod2sparse* getInvHc() const { return invHc; }
    int getM() const { return M; }
    int getN() const { return N; }
};

#endif // LDPC_DECODE_H
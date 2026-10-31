/* RCODE.C - Procedures to read parity check and generator matrices. */

/* Copyright (c) 2000, 2001 by Radford M. Neal 
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
#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "alloc.h"
#include "intio.h"
#include "open.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "check.h"

#include "ldpc_decode.h"

/* CONSTRUCTOR */
LDPCDecode::LDPCDecode() 
{
    H = NULL;
    M = 0;
    N = 0;
    type = 0;
    cols = NULL;
    L = NULL;
    U = NULL;
    rows = NULL;
    G = NULL;
    Hm = NULL;
    mM = 0;
    mN = 0;
    invHc = NULL;
    cM = 0;
    cN = 0;
}

/* DESTRUCTOR */
LDPCDecode::~LDPCDecode() 
{
    if (H) mod2sparse_free(H);
    if (L) mod2sparse_free(L);
    if (U) mod2sparse_free(U);
    if (G) mod2dense_free(G);
    if (Hm) mod2sparse_free(Hm);
    if (invHc) mod2sparse_free(invHc);
    if (cols) free(cols);
    if (rows) free(rows);
}

void LDPCDecode::readEncpara(char *pchk_Hm, char *pchk_Hc)
{
    FILE *f;

    f = open_file_std(pchk_Hm,"rb");
    if (f==NULL)
    { fprintf(stderr,"Can't open Hm file: %s\n",pchk_Hm);
        exit(1);
    }

    if (intio_read(f)!=('P'<<8)+0x80)
    { fprintf(stderr,"File %s doesn't contain a parity check matrix\n",pchk_Hm);
        exit(1);
    }

    Hm = (mod2sparse*)mod2sparse_read(f);

    if (Hm==0)
    { fprintf(stderr,"Error reading Hm from %s\n",pchk_Hm);
        exit(1);
    }

    mM = mod2sparse_rows(Hm);
    mN = mod2sparse_cols(Hm);

    fclose(f);

    f = open_file_std(pchk_Hc,"rb");
    if (f==NULL)
    { fprintf(stderr,"Can't open Hm file: %s\n",pchk_Hc);
        exit(1);
    }

    if (intio_read(f)!=('P'<<8)+0x80)
    { fprintf(stderr,"File %s doesn't contain a parity check matrix\n",pchk_Hc);
        exit(1);
    }

    invHc = (mod2sparse*)mod2sparse_read(f);

    if (invHc==0)
    { fprintf(stderr,"Error reading Hm from %s\n",pchk_Hc);
        exit(1);
    }

    cM = mod2sparse_rows(invHc);
    cN = mod2sparse_cols(invHc);

    fclose(f);
}

void LDPCDecode::readPchk(char *pchk_file)
{
    FILE *f;

    f = open_file_std(pchk_file,"rb");
    if (f==NULL)
    { fprintf(stderr,"Can't open parity check file: %s\n",pchk_file);
        exit(1);
    }

    if (intio_read(f)!=('P'<<8)+0x80)
    { fprintf(stderr,"File %s doesn't contain a parity check matrix\n",pchk_file);
        exit(1);
    }

    H = (mod2sparse*)mod2sparse_read(f);

    if (H==0)
    { fprintf(stderr,"Error reading parity check matrix from %s\n",pchk_file);
        exit(1);
    }

    M = mod2sparse_rows(H);
    N = mod2sparse_cols(H);

    fclose(f);
}

void LDPCDecode::readGen(char *gen_file, int cols_only, int no_pchk_file)
{
    int M2, N2;
    FILE *f;
    int i;

    f = open_file_std(gen_file,"rb");
    if (f==NULL)
    { fprintf(stderr,"Can't open generator matrix file: %s\n",gen_file);
        exit(1);
    }

    if (intio_read(f)!=('G'<<8)+0x80)
    { fprintf(stderr,"File %s doesn't contain a generator matrix\n",gen_file);
        exit(1);
    }

    if (fread (&type, 1, 1, f) != 1) goto error;

    M2 = intio_read(f);
    N2 = intio_read(f);

    if (feof(f) || ferror(f)) goto error;

    if (no_pchk_file)
    { M = M2;
        N = N2;
    }
    else 
    { if (M2!=M || N2!=N)
        { fprintf(stderr,
                  "Generator matrix and parity-check matrix are incompatible\n");
          exit(1);
        }
    }

    cols = (int *) chk_alloc (N, sizeof *cols);
    rows =(int *) chk_alloc (M, sizeof *rows);

    for (i = 0; i<N; i++)
    { 
        cols[i] = intio_read(f);
        if (feof(f) || ferror(f)) goto error;
    }

    if (!cols_only)
    {
        switch (type)
        {
          case 's':
          { 
            for (i = 0; i<M; i++)
            { rows[i] = intio_read(f);
              if (feof(f) || ferror(f)) goto error;
            }

            if ((L = (mod2sparse*)mod2sparse_read(f)) == 0) goto error;
            if ((U = (mod2sparse*)mod2sparse_read(f)) == 0) goto error;
      
            if (mod2sparse_rows(L)!=M || mod2sparse_cols(L)!=M) goto garbled;
            if (mod2sparse_rows(U)!=M || mod2sparse_cols(U)<M) goto garbled;
           
            break;
          }
      
          case 'd':
          {
            if ((G = mod2dense_read(f)) == 0) goto error;
      
            if (mod2dense_rows(G)!=M || mod2dense_cols(G)!=N-M) goto garbled;
      
            break;
          }
      
          case 'm':
          {
            if ((G = mod2dense_read(f)) == 0) goto error;
      
            if (mod2dense_rows(G)!=M || mod2dense_cols(G)!=M) goto garbled;
      
            break;
          }
      
          default: 
          { fprintf(stderr,
             "Unknown type of generator matrix in file %s\n",gen_file);
            exit(1);
          }
        }
    }
      
    fclose(f);
    return;

error:
    fprintf(stderr,"Error reading generator matrix from file %s\n",gen_file);
    exit(1);

garbled:
    fprintf(stderr,"Garbled generator matrix in file %s\n",gen_file);
    exit(1);
}

void LDPCDecode::fftgf64(int* index, double* din, double * dout)
{
    int i,j;
    for(i=0;i<64;i++)
    {
        *(dout+i) = 0.0;
        for(j=0;j<64;j++)
        {
            if(  *(index+j+64*i)  ==1  )		   
                *(dout+i) += *(din+j);
            else
                *(dout+i) -= *(din+j);
        }
    }	
}

void LDPCDecode::encodecom(int* indexg, int* src,int * code)
{
    int i,j;
    int t;
    for(i=0;i<18;i++)
    {
        *(code+i+6) = *(src+i) ;
    }
    for(i=0;i<6;i++)
    {
        t= 0;
        for(j=0;j<18;j++)
        {
            if(  *(indexg+j+18*i)  == 1  )		   
                t += *(src+j);
        }
        t= t %2;		
        *(code+i) =t; 
    }	
}

void LDPCDecode::de2bi(int w, int din, char *dout)
{
    int i;
    int t=din;
    for(i=0;i<w;i++)
    {
        if(t >= int (pow(float(2.0),w-i-1)  )   )
        {		
            *(dout+w-i-1) =1;
            t -= int (pow(float(2.0),w-i-1)  ) ;
        }
        else
           *(dout+w-i-1) =0;
    }
}

int LDPCDecode::bi2de(int w, int *din)
{
    int i;
    int dout=0;
    for(i=0;i<w;i++)
    {
        dout += (*(din+i))*(  int( pow(float(2.0),i) ) );  // pow  
    }
    return dout;
}

void LDPCDecode::cnu(double *p, double *pifftdi, int * index, int *permutation, int *habc)
{
    int j,o;
    double ffti[64];
    double ffto[64];
    double pi[64*4];
    double pifft[64*4];
    double pifftp[64*4];
    double pifftd[64*4];

     // transform of h_i(u)  
     for(j=0;j<64;j++)
     {
         *(pi+(*(habc+permutation[0]*64*2 +j))+0*64) =*(p+j+0*64);
         *(pi+(*(habc+permutation[1]*64*2 +j))+1*64) =*(p+j+1*64);
         *(pi+(*(habc+permutation[2]*64*2 +j))+2*64) =*(p+j+2*64);
         *(pi+(*(habc+permutation[3]*64*2 +j))+3*64) =*(p+j+3*64);
     }

     // FFT, indeed fast Hadamard transform
     for(j=0;j<4;j++)
     {
         for(o=0;o<64;o++)
         {				
            *(ffti+o) = *(pi+o+64*j);
         }
         fftgf64(index,ffti,ffto);
         for(o=0;o<64;o++)
         {				
            *(pifft+o+64*j)= *(ffto+o);
         }
     }

     // product, exactly dot product
     for(j=0;j<64;j++)
     {	
         *(pifftp+0*64+j)= (*(pifft+1*64+j))*(*(pifft+2*64+j))*(*(pifft+3*64+j)) ;
         *(pifftp+1*64+j)= (*(pifft+0*64+j))*(*(pifft+2*64+j))*(*(pifft+3*64+j)) ;
         *(pifftp+2*64+j)= (*(pifft+0*64+j))*(*(pifft+1*64+j))*(*(pifft+3*64+j)) ;
         *(pifftp+3*64+j)= (*(pifft+0*64+j))*(*(pifft+1*64+j))*(*(pifft+2*64+j)) ;
     }

   //  inversion of FFT (IFFT), here we use '-d' for IFFT
     for(j=0;j<4;j++)
     {
         for(o=0;o<64;o++)
         {				
            *(ffti+o) = *(pifftp+o+64*j);
         }
         fftgf64(index,ffti,ffto);

         for(o=0;o<64;o++)
         {				
            *(pifftd+o+64*j)= *(ffto+o)/64;
         }
     }

     // inversion of h_i(u)   
     for(j=0;j<64;j++)
     {
         *(pifftdi+(*(habc+permutation[0]*64*2  +64 +j)) +0*64) =*(pifftd+j+0*64);
         *(pifftdi+(*(habc+permutation[1]*64*2  +64 +j)) +1*64) =*(pifftd+j+1*64);
         *(pifftdi+(*(habc+permutation[2]*64*2  +64 +j)) +2*64) =*(pifftd+j+2*64);
         *(pifftdi+(*(habc+permutation[3]*64*2  +64 +j)) +3*64) =*(pifftd+j+3*64);
     }
}

void LDPCDecode::cnu3(double *p, double *pifftdi, int * index, int *permutation, int *habc)
{
    int j,o;
    double ffti[64];
    double ffto[64];
    double pi[64*4];
    double pifft[64*4];
    double pifftp[64*4];
    double pifftd[64*4];

     // transform of h_i(u)  
     for(j=0;j<64;j++)
     {
         *(pi+(*(habc+permutation[0]*64*2 +j))+0*64) =*(p+j+0*64);
         *(pi+(*(habc+permutation[1]*64*2 +j))+1*64) =*(p+j+1*64);
         *(pi+(*(habc+permutation[2]*64*2 +j))+2*64) =*(p+j+2*64);
     }

     // FFT, indeed fast Hadamard transform
     for(j=0;j<3;j++)
     {
         for(o=0;o<64;o++)
         {				
            *(ffti+o) = *(pi+o+64*j);
         }
         fftgf64(index,ffti,ffto);
         for(o=0;o<64;o++)
         {				
            *(pifft+o+64*j)= *(ffto+o);
         }
     }

     // product, exactly dot product
     for(j=0;j<64;j++)
     {	
         *(pifftp+0*64+j)= (*(pifft+1*64+j))*(*(pifft+2*64+j));
         *(pifftp+1*64+j)= (*(pifft+0*64+j))*(*(pifft+2*64+j));
         *(pifftp+2*64+j)= (*(pifft+0*64+j))*(*(pifft+1*64+j));
     }

   //  inversion of FFT (IFFT), here we use '-d' for IFFT
     for(j=0;j<3;j++)
     {
         for(o=0;o<64;o++)
         {				
            *(ffti+o) = *(pifftp+o+64*j);
         }
         fftgf64(index,ffti,ffto);

         for(o=0;o<64;o++)
         {				
            *(pifftd+o+64*j)= *(ffto+o)/64;
         }
     }

     // inversion of h_i(u)   
     for(j=0;j<64;j++)
     {
         *(pifftdi+(*(habc+permutation[0]*64*2  +64 +j)) +0*64) =*(pifftd+j+0*64);
         *(pifftdi+(*(habc+permutation[1]*64*2  +64 +j)) +1*64) =*(pifftd+j+1*64);
         *(pifftdi+(*(habc+permutation[2]*64*2  +64 +j)) +2*64) =*(pifftd+j+2*64);
     }
}

void LDPCDecode::vnu(double *c2v, double *prior, double *posterior, char *dec, double *v2c)
{
    int i,j;
    double pp1;
    int pos1;
    char decs[6];
    double pifftdiv[64];

    for(i=0;i<52;i++)
    {
         for(j=0;j<64;j++)
         {
             *(pifftdiv+j)    =  ( *(c2v+i*64  + j + 64*52) )  * (*(prior+i*64+j));						
         }
         // accumulation			
         pp1=0;
         for(j=0;j<64;j++)
         {
             pp1 += *(pifftdiv+j);
         }

         // uniform
         for(j=0;j<64;j++)
         {
             *(pifftdiv+j) /=pp1;
         }
                     
         for(j=0;j<64;j++)
         {
              *(v2c+i*64+j)  =  *(pifftdiv+j);				
         }


         //--------------------------------------------
         for(j=0;j<64;j++)
         {
             *(pifftdiv+j)    =  ( *(c2v+i*64  + j ) )  * (*(prior+i*64+j));						
         }
         // accumulation			
         pp1=0;
         for(j=0;j<64;j++)
         {
             pp1 += *(pifftdiv+j);
         }

         // uniform
         for(j=0;j<64;j++)
         {
             *(pifftdiv+j) /=pp1;
         }
                     
         for(j=0;j<64;j++)
         {
               *(v2c+i*64+j+ 64*52) =  *(pifftdiv+j);				
         }

         // variable node update, denoted as '-v'
         for(j=0;j<64;j++)
         {
             *(pifftdiv+j)=  *(c2v+i*64+j) * ( *(c2v+i*64 +j + 64*52) )  * (*(prior+i*64+j));				
         }

         // accumulation			
         pp1=0;
         for(j=0;j<64;j++)
         {
             pp1 += *(pifftdiv+j);
         }

         // uniform
         for(j=0;j<64;j++)
         {
             *(pifftdiv+j) /=pp1;
         }
                     
         for(j=0;j<64;j++)
         {
              (*(posterior+i*64+j))=  *(pifftdiv+j);				
         }

         pos1=0;
         pp1 = *(pifftdiv+0 );

         for(j=0;j<63;j++)
         {
             if( *(pifftdiv+j+1) > pp1  )
             {    
                 pp1 = *(pifftdiv+j+1+0*64);
                 pos1  = j+1;
             }
         }

         // converting the position of maximum probability to 
         // bit vector, here MSB is on the right
         // [LSB ... MSB]
         de2bi(6,pos1,decs);
         for(j=0;j<6;j++)
         {
             *(dec+j+6*i)  =  *(decs+j);
         }
    }
}

void LDPCDecode::decode(char *dec, int* success, int* iteration, double *allf1, double *allf0, 
                        int iternum, int *index, int *habc, int *checknode)// mod2sparse *H
{
    double prior[64*312/6];	
    double posterior[64*312/6];	
    double v2c[2*64*52];	
    double c2v[2*64*52];	
    double f1[24] = {0};
    double f0[24] = { 0 };

    int permutation[4];
    double p[64*4];
    double pifftdi[64*4];
    double pp1;
    int row,col,val, rowaddr,coladdr;
    char check_pt[156];
    int i,j,o;
    int num;
    char binseq[6] = {0};
    char* binseq_pt = binseq;

    for(i=0;i<312/6;i++)   // loop for each symbol
    {

        for(j=0;j<6;j++)
        {
            *(f0+j)  = *(allf0+j+i*6);
            *(f1+j)  = *(allf1+j+i*6);
        }
        // probability vector from symbol probability of each bit
        for(j=0;j<64;j++)
        {
            de2bi(6, j, binseq_pt);

            pp1=1;
            for(o=0;o<6;o++)
            {	 
                if(*(binseq_pt+o)==0)
                {
                    pp1 *= *(f0 +o) ;
                }
                else
                {
                    pp1 *= *(f1 +o) ;
                }				 			 
            }
            *(p+j) =pp1;   // only use part
        
        }
        // uniform
        pp1=0;
        for(j=0;j<64;j++)
        {
            pp1 += *(p+j);
        }
        for(j=0;j<64;j++)
        {
            *(p+j) /=pp1;
        }
        for(j=0;j<64;j++)
        {
            *(prior+64*i+j) = *(p+j) ;
        }
    }

    // c2v initilization with all zero
    for(i=0;i<2;i++)
    {
        for(j=0;j<312/6;j++)
        {
            for(o=0;o<64;o++)
            {
                *(c2v+ i*64*312/6 +j*64+o) = 1.0;
            }
        }
    }

    // vnu and decision
    // vnu with the extrinsic information
    vnu(c2v,prior,posterior,dec,v2c);

    success[0] = check(H,dec,check_pt)==0;
    // printf("Before Iteration : %d\n", success[0]);

    num = 0;
    while( success[0]== 0 && num < iternum)
    {
        i=0;
        {
            for (j=0;j<3;j++)
            {
                row = i*4*3 + 0 + j*3;
                col = i*4*3 + 1 + j*3;
                val = i*4*3 + 2 + j*3;
                rowaddr = checknode[row]-1;
                coladdr = checknode[col]-1;
            
                permutation[j] = checknode[val]-1;
            
                for(o=0;o<64;o++)
                {
                    *(p+ j*64 + o) = v2c[rowaddr* 52*64 + coladdr*64 + o];
                }
            } 
                
            cnu3(p,pifftdi,index,permutation,habc);

            for (j=0;j<3;j++)
            {
                row = i*4*3 + 0 + j*3;
                col = i*4*3 + 1 + j*3;
            
                rowaddr = checknode[row]-1;
                coladdr = checknode[col]-1;
            
                for(o=0;o<64;o++)
                {
                    c2v[rowaddr* 52*64 + coladdr*64 + o] =  *(pifftdi+ j*64 + o) ;
                }
            } 					 
        }

        for(i=1;i<26;i++)
        {
            // each component cnu
            for (j=0;j<4;j++)
            {
                row = i*4*3 + 0 + j*3;
                col = i*4*3 + 1 + j*3;
                val = i*4*3 + 2 + j*3;
                rowaddr = checknode[row]-1;
                coladdr = checknode[col]-1;
                permutation[j] = checknode[val]-1;
            
                for(o=0;o<64;o++)
                {
                    *(p+ j*64 + o) = v2c[rowaddr* 52*64 + coladdr*64 + o];
                }
            } 

            cnu(p,pifftdi,index,permutation,habc);

            // put in the memory for c2v
            // pifftdi  ->  c2v
            for (j=0;j<4;j++)
            {
                row = i*4*3 + 0 + j*3;
                col = i*4*3 + 1 + j*3;
            
                rowaddr = checknode[row]-1;
                coladdr = checknode[col]-1;
            
                for(o=0;o<64;o++)
                {
                    c2v[rowaddr* 52*64 + coladdr*64 + o] =  *(pifftdi+ j*64 + o) ;
                }
            } 				 
        }

    //vnu
        vnu(c2v,prior,posterior,dec,v2c);


        num++;
    
        success[0] = check((mod2sparse*)H,dec,check_pt)==0;  // decoding success sign

        // printf("iteratin num %d,  Success  %d \n\n",num, success[0]);
    }
        
    iteration[0] = num;                 // practical iteration number 
}



void LDPCDecode:: encode (mod2sparse *Hm,mod2sparse *invHc,char  *src,char  *code)
{
    int i;
    char temp[156];
    char checkbit[156];
    
    mod2sparse_mulvec(Hm, src,temp);
    mod2sparse_mulvec(invHc, temp,checkbit);
    
    for(i=0;i<156;i++)
    {
        *(code+i) = *(src+i);
    }
    for(i=0;i<156;i++)
    {
        *(code+i+156) = *(checkbit+i);
    }
}

#ifdef __cplusplus
}
#endif
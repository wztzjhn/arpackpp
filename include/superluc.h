/*
   ARPACK++ v1.2 2/20/2000
   c++ MKL_INTerface to ARPACK code.

   MODULE SuperLUc.h.
   Interface to SuperLU routines.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef SUPERLUC_H
#define SUPERLUC_H

#include "arch.h"
#include "arlspdef.h"
#include "arlsupm.h"
#include "arlcomp.h"

// gstrf.

inline void gstrf(superlu_options_t *options, SuperMatrix *A,
        MKL_INT relax, MKL_INT panel_size, MKL_INT *etree, void *work, MKL_INT lwork,
        MKL_INT *perm_c, MKL_INT *perm_r, SuperMatrix *L, SuperMatrix *U,
        SuperLUStat_t *stat, MKL_INT *info)
{
  if (A->Dtype == SLU_D) {       // calling the double precision routine.
    dGlobalLU_t Glu;
    dgstrf(options,A,relax,
           panel_size,etree,work,lwork,perm_c,perm_r,L,U,&Glu,stat,info);
  }
  else if (A->Dtype == SLU_S) {  // calling the single precision routine.
    sGlobalLU_t Glu;
    sgstrf(options,A,relax,
           panel_size,etree,work,lwork,perm_c,perm_r,L,U,&Glu,stat,info);
  }
  else if (A->Dtype == SLU_Z) {  // calling the double precision complex routine.
#ifdef ARCOMP_H
    zGlobalLU_t Glu;
    zgstrf(options,A,relax,
           panel_size,etree,work,lwork,perm_c,perm_r,L,U,&Glu,stat,info);
#endif
  }
  else {                      // calling the single precision complex routine.
#ifdef ARCOMP_H
    cGlobalLU_t Glu;
    cgstrf(options,A,relax,
           panel_size,etree,work,lwork,perm_c,perm_r,L,U,&Glu,stat,info);
#endif
  }

} // gstrf.


inline void gstrs(trans_t trans, SuperMatrix *L, SuperMatrix *U,
	          MKL_INT *perm_c, MKL_INT *perm_r, SuperMatrix *B, SuperLUStat_t* stat, MKL_INT *info)
{

  if (L->Dtype == SLU_D) {       // calling the double precision routine.
    dgstrs(trans,L,U,perm_c,perm_r,B,stat,info);
  }
  else if (L->Dtype == SLU_S) {  // calling the single precision routine.
    sgstrs(trans,L,U,perm_c,perm_r,B,stat,info);
  }
  else if (L->Dtype == SLU_Z) {  // calling the double precision complex routine.
#ifdef ARCOMP_H
    zgstrs(trans,L,U,perm_c,perm_r,B,stat,info);
#endif
  }
  else {                      // calling the single precision complex routine.
#ifdef ARCOMP_H
    cgstrs(trans,L,U,perm_c,perm_r,B,stat,info);
#endif
  }

} // gstrs.


// Create_CompCol_Matrix.

inline void Create_CompCol_Matrix(SuperMatrix* A, MKL_INT m, MKL_INT n, MKL_INT nnz,
                                  double* a, MKL_INT* irow, MKL_INT* pcol,
                                  Stype_t S, Mtype_t M)
{

  dCreate_CompCol_Matrix(A,m,n,nnz,a,irow,pcol,S,SLU_D,M);

} // Create_CompCol_Matrix (double).

inline void Create_CompCol_Matrix(SuperMatrix* A, MKL_INT m, MKL_INT n, MKL_INT nnz,
                                  float* a, MKL_INT* irow, MKL_INT* pcol,
                                  Stype_t S, Mtype_t M)
{

  sCreate_CompCol_Matrix(A,m,n,nnz,a,irow,pcol,S,SLU_S,M);

} // Create_CompCol_Matrix (float).

#ifdef ARCOMP_H

inline void Create_CompCol_Matrix(SuperMatrix* A, MKL_INT m, MKL_INT n, MKL_INT nnz,
                                  arcomplex<double>* a, MKL_INT* irow, MKL_INT* pcol,
                                  Stype_t S, Mtype_t M)
{

  zCreate_CompCol_Matrix(A,m,n,nnz,(ldcomplex*)a,irow,pcol,S,SLU_Z,M);

} // Create_CompCol_Matrix (complex<double>).

inline void Create_CompCol_Matrix(SuperMatrix* A, MKL_INT m, MKL_INT n, MKL_INT nnz,
                                  arcomplex<float>* a, MKL_INT* irow, MKL_INT* pcol,
                                  Stype_t S, Mtype_t M)
{

  cCreate_CompCol_Matrix(A,m,n,nnz,(lscomplex*)a,irow,pcol,S,SLU_C,M);

} // Create_CompCol_Matrix (complex<float>).

#endif // ARCOMP_H.


// Create_Dense_Matrix.

inline void Create_Dense_Matrix(SuperMatrix* A, MKL_INT m, MKL_INT n, double* x,
                                MKL_INT ldx, Stype_t S, Mtype_t M)
{

  dCreate_Dense_Matrix(A,m,n,x,ldx,S,SLU_D,M);

} // Create_Dense_Matrix (double).

inline void Create_Dense_Matrix(SuperMatrix* A, MKL_INT m, MKL_INT n, float* x,
                                MKL_INT ldx, Stype_t S, Mtype_t M)
{

  sCreate_Dense_Matrix(A,m,n,x,ldx,S,SLU_S,M);

} // Create_Dense_Matrix (float).

#ifdef ARCOMP_H

inline void Create_Dense_Matrix(SuperMatrix* A, MKL_INT m, MKL_INT n, arcomplex<double>* x,
                                MKL_INT ldx, Stype_t S, Mtype_t M)
{

  zCreate_Dense_Matrix(A,m,n,(ldcomplex*)x,ldx,S,SLU_Z,M);

} // Create_Dense_Matrix (complex<double>).

inline void Create_Dense_Matrix(SuperMatrix* A, MKL_INT m, MKL_INT n, arcomplex<float>* x,
                                MKL_INT ldx, Stype_t S, Mtype_t M)
{

  cCreate_Dense_Matrix(A,m,n,(lscomplex*)x,ldx,S,SLU_C,M);

} // Create_Dense_Matrix (complex<float>).

#endif // ARCOMP_H.

#endif // SUPERLUC_H

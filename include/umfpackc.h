/*
   ARPACK++ v1.2 2/20/2000
   c++ MKL_INTerface to ARPACK code.

   MODULE UMFPACKc.h.
   Interface to UMFPACK routines.

   Author of this class:
      Martin Reuter
      Date 2/28/2013
      
   Arpack++ Author:
      Francisco Gomes
      
   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef UMFPACKC_H
#define UMFPACKC_H

#ifdef __cplusplus
extern "C" {
#endif

#define UMFPACK_INFO 90
#define UMFPACK_CONTROL 20
#define UMFPACK_OK (0)
#define UMFPACK_A	(0)	/* Ax=b    */
#define UMFPACK_PRL 0			/* print level */

void umfpack_di_defaults
(
    double Control [UMFPACK_CONTROL]
) ;


int umfpack_di_symbolic
(
    MKL_INT n_row,
    MKL_INT n_col,
    const MKL_INT Ap [ ],
    const MKL_INT Ai [ ],
    const double Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_di_numeric
(
    const MKL_INT Ap [ ],
    const MKL_INT Ai [ ],
    const double Ax [ ],
    void *Symbolic,
    void **Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

void umfpack_di_free_symbolic
(
    void **Symbolic
) ;

void umfpack_di_free_numeric
(
    void **Numeric
) ;

int umfpack_di_triplet_to_col
(
    MKL_INT n_row,
    MKL_INT n_col,
    MKL_INT nz,
    const MKL_INT Ti [ ],
    const MKL_INT Tj [ ],
    const double Tx [ ],
    MKL_INT Ap [ ],
    MKL_INT Ai [ ],
    double Ax [ ],
    MKL_INT Map [ ]
) ;

int umfpack_di_solve
(
    MKL_INT sys,
    const MKL_INT Ap [ ],
    const MKL_INT Ai [ ],
    const double Ax [ ],
    double X [ ],
    const double B [ ],
    void *Numeric,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]
) ;

int umfpack_di_report_matrix
(
    MKL_INT n_row,
    MKL_INT n_col,
    const MKL_INT Ap [ ],
    const MKL_INT Ai [ ],
    const double Ax [ ],
    MKL_INT col_form,
    const double Control [UMFPACK_CONTROL]
) ;

#ifdef __cplusplus
  }
#endif

//#include "umfpack.h"
#include <fstream>

inline void Write_Triplet_Matrix(const std::string & fname, MKL_INT * tripi,
                                 MKL_INT * tripj, double* tripx, unsigned MKL_INT nnz)
{
  std::ofstream myfile; 
  myfile.open ( fname.c_str() );
	myfile.precision(20);
  for (unsigned MKL_INT i=0;i<nnz;i++)
  {
    myfile << tripi[i]+1 << " " << tripj[i]+1 << " " << tripx[i] << std::endl;
  }
  myfile.close();
}

/*inline void Write_Cholmod_Sparse_Matrix(const std::string & fname,
                             cholmod_sparse* A, cholmod_common *c)
{
  std::ofstream myfile; 
  myfile.open ( fname.c_str() );
  cholmod_triplet * T = cholmod_sparse_to_triplet(A,c);
  //std::cout << " [ " << std::endl;
	myfile.precision(20);
  for (unsigned MKL_INT i=0;i<T->nnz;i++)
  {
    myfile << ((MKL_INT*)T->i)[i]+1 << " " << ((MKL_INT*)T->j)[i]+1 << " " << ((double*)T->x)[i] << std::endl;
  }
  //std::cout << " ] " << std::endl;
  myfile.close();
  
  cholmod_free_triplet(&T,c);

}

// Create_Cholmod_Sparse_Matrix 
inline cholmod_sparse* Create_Cholmod_Sparse_Matrix(MKL_INT m, MKL_INT n, MKL_INT nnz,
      double* a, MKL_INT* irow, MKL_INT* pcol, char uplo, cholmod_common *c)
{
  
  cholmod_sparse* A = new cholmod_sparse;
  A->nrow = m;
  A->ncol = n;
  A->nzmax = nnz;
  A->p = pcol;
  A->i = irow;
  A->nz = NULL;
  A->x = a;
  A->z = NULL;
  if (uplo == 'L') A->stype = -1;
  else A->stype = 1;
  A->itype = CHOLMOD_INT;
  A->xtype = CHOLMOD_REAL; // real
  A->dtype = CHOLMOD_DOUBLE; // double
  A->sorted = 0;
  A->packed = 1;

  return A;  
  
  

  
} // Create_Cholmod_Sparse_Matrix (double).

// Create_Cholmod_Dense_Matrix (from Triplet)
inline cholmod_dense* Create_Cholmod_Dense_Matrix(MKL_INT m, MKL_INT n,
                                  double* a, cholmod_common *c)
{


  cholmod_dense* A = new cholmod_dense;
  A->nrow = m;
  A->ncol = n;
  A->nzmax = m*n;
  A->d = m;
  A->x = a;
  A->z = NULL;
  A->xtype = CHOLMOD_REAL; // real
  A->dtype = CHOLMOD_DOUBLE; // double

//  cholmod_dense* As = cholmod_copy_dense(A,c);
  
  return A;
  
} // Create_Cholmod_Dense_Matrix (double).

// Create_Cholmod_Dense_Matrix (from Triplet)
inline void Get_Cholmod_Dense_Data(cholmod_dense* A, MKL_INT n, double* a)
{
  memcpy(a,A->x,n*sizeof(double));
  
//  for (MKL_INT i = 0;i<n;i++)
//    a[i] = ((double*)A->x)[i];
  
} // Create_Cholmod_Dense_Matrix (double).

*/

#endif // UMFPACKC_H

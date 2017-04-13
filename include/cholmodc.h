/*
   ARPACK++ v1.2 2/20/2000
   c++ MKL_INTerface to ARPACK code.

   MODULE CHOLMODc.h.
   Interface to CHOLMOD routines.

   Author of this class:
      Martin Reuter
      Date 11/05/2012
      
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

#ifndef CHOLMODC_H
#define CHOLMODC_H

#include "cholmod.h"
#include <fstream>

inline void Write_Cholmod_Sparse_Matrix(const std::string & fname,
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
  
  
/*  double* hd  = new double[nnz];
  MKL_INT* hi     = new MKL_INT[nnz];
  MKL_INT* hp     = new MKL_INT[n+1];
  
  MKL_INT col,j;
  MKL_INT counter=0;
  MKL_INT counter2=0;
  for (col=0;col<n;col++) // column
  {
    hp[col] = counter2;
    for (j=pcol[col];j<pcol[col+1];j++)
    {
      MKL_INT & row = irow[counter];
      if ((uplo == 'L' && row >= col) ||(uplo == 'U' && row <= col))
      {
        hd[counter2] = a[counter];
        hi[counter2] = irow[counter];
        counter2++;
        //std::cout << " In : " << std::flush;
      }
      //else  std::cout << " Out : " << std::flush;

      //std::cout << row+1 << " " << col+1 << " " << a[counter] << std::endl;
      counter++;
    }
  
  }
  hp[n] = counter2;
  
  
  cholmod_sparse* A = new cholmod_sparse;
  A->nrow = m;
  A->ncol = n;
  A->nzmax = counter2;
  A->p = hp;
  A->i = hi;
  A->nz = NULL;
  A->x = hd;
  A->z = NULL;
  if (uplo == 'L') A->stype = -1;
  else A->stype = 1;
  A->itype = CHOLMOD_INT;
  A->xtype = CHOLMOD_REAL; // real
  A->dtype = CHOLMOD_DOUBLE; // double
  A->sorted = 0;
  A->packed = 1;
  
  //cholmod_sparse* As = cholmod_copy_sparse(A,c);

  return A;*/
  
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



#endif // CHOLMODC_H

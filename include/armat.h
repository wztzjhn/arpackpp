/*
   ARPACK++ v1.2 2/20/2000
   c++ MKL_INTerface to ARPACK code.

   MODULE ARMat.h
   Generic matrix template with a matrix-vector product.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARMAT_H
#define ARMAT_H

template<class ARTYPE>
class ARMatrix {

 protected:

  MKL_INT  m, n;    // Number of rows and columns.
  bool defined;
 
 public:

  ARMatrix() { defined = false; }
  // Short constructor.

  ARMatrix(MKL_INT nrows, MKL_INT ncols = 0)
  // Long constructor.
  {
    m = nrows;
    n = (ncols?ncols:nrows);
    defined = false;
  } // Constructor.

  virtual ~ARMatrix() { }
  // Destructor.

  MKL_INT nrows() { return m; }

  MKL_INT ncols() { return n; }

  bool IsDefined() { return defined; }

  virtual void MultMv(ARTYPE* v, ARTYPE* w) = 0;
  // Matrix-vector product: w = A*v.

}; // ARMatrix.

#endif // ARMAT_H


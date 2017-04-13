/*
   ARPACK++ v1.2 2/20/2000
   c++ MKL_INTerface to ARPACK code.

   MODULE ARRSComp.h.
   Arpack++ class ARrcCompStdEig definition.

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARRSCOMP_H
#define ARRSCOMP_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arerror.h"
#include "debug.h"
#include "arrseig.h"
#include "caupp.h"
#include "ceupp.h"

template<class ARFLOAT>
class ARrcCompStdEig: virtual public ARrcStdEig<ARFLOAT, arcomplex<ARFLOAT> > {

 protected:

 // a) Protected functions:

 // a.1) Memory control functions.

  void WorkspaceAllocate();
  // Allocates workspace for complex problems.


 // a.2) Functions that handle original FORTRAN ARPACK code.

  void Aupp();
  // Interface to FORTRAN subroutines CNAUPD and ZNAUPD.

  void Eupp();
  // Interface to FORTRAN subroutines CNEUPD and ZNEUPD.

 public:

 // b) Public functions:

 // b.1) Trace functions.

  void Trace(const MKL_INT digit = -5, const MKL_INT getv0 = 0, const MKL_INT aupd = 1,
             const MKL_INT aup2 = 0,  const MKL_INT aitr = 0,  const MKL_INT eigt = 0,
             const MKL_INT apps = 0,  const MKL_INT gets = 0,  const MKL_INT eupd = 0)
  { 
    cTraceOn(digit, getv0, aupd, aup2, aitr, eigt, apps, gets, eupd); 
  }
  // Turns on trace mode. 


 // b.2) Functions that perform all calculations in one step.

  MKL_INT Eigenvalues(arcomplex<ARFLOAT>* &EigValp, bool ivec = false,
                  bool ischur = false);
  // Overrides array EigValp with the eigenvalues of the problem.
  // Also calculates eigenvectors and Schur vectors if requested.

  MKL_INT EigenValVectors(arcomplex<ARFLOAT>* &EigVecp, 
                      arcomplex<ARFLOAT>* &EigValp, bool ischur = false);
  // Overrides array EigVecp sequentially with the eigenvectors of the
  // given eigen-problem. Also stores the eigenvalues in EigValp.
  // Calculates Schur vectors if requested.


 // b.3) Functions that return elements of vectors and matrices.

  arcomplex<ARFLOAT> Eigenvalue(MKL_INT i);
  // Provides i-eth eigenvalue.

  arcomplex<ARFLOAT> Eigenvector(MKL_INT i, MKL_INT j);
  // Provides element j of the i-eth eigenvector.


 // b.4) Functions that use STL vector class.

#ifdef STL_VECTOR_H

  vector<arcomplex<ARFLOAT> >* StlEigenvalues(bool ivec = false,
                                              bool ischur = false);
  // Calculates the eigenvalues and stores them in a single STL vector.
  // Also calculates eigenvectors and Schur vectors if requested.

  vector<arcomplex<ARFLOAT> >* StlEigenvector(MKL_INT i);
  // Returns the i-th eigenvector in a STL vector.

#endif // #ifdef STL_VECTOR_H.


 // b.5) Constructors and destructor.

  ARrcCompStdEig() { }
  // Short constructor.

  ARrcCompStdEig(MKL_INT np, MKL_INT nevp, const std::string& whichp = "LM",
                 MKL_INT ncvp = 0, ARFLOAT tolp = 0.0, MKL_INT maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARrcCompStdEig(MKL_INT np, MKL_INT nevp, arcomplex<ARFLOAT> sigma,
                 const std::string& whichp = "LM", MKL_INT ncvp = 0, ARFLOAT tolp = 0.0,
                 MKL_INT maxitp = 0, arcomplex<ARFLOAT>* residp = NULL,
                 bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARrcCompStdEig(const ARrcCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARrcCompStdEig() { }
  // Destructor.

 // c) Operators.

  ARrcCompStdEig& operator=(const ARrcCompStdEig& other);
  // Assignment operator.

}; // class ARrcCompStdEig.


// ------------------------------------------------------------------------ //
// ARrcCompStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARrcCompStdEig<ARFLOAT>::WorkspaceAllocate()
{

  this->lworkl  = this->ncv*(3*this->ncv+6);
  this->lworkv  = 2*this->ncv;
  this->lrwork  = this->ncv;
  this->workl   = new arcomplex<ARFLOAT>[this->lworkl+1];
  this->workv   = new arcomplex<ARFLOAT>[this->lworkv+1];
  this->rwork   = new ARFLOAT[this->lrwork+1];

} // WorkspaceAllocate.


template<class ARFLOAT>
inline void ARrcCompStdEig<ARFLOAT>::Aupp()
{

  caupp(this->ido, this->bmat, this->n, this->which, this->nev, this->tol, this->resid, this->ncv, this->V, this->n,
        this->iparam, this->ipntr, this->workd, this->workl, this->lworkl, this->rwork, this->info);

} // Aupp.


template<class ARFLOAT>
inline void ARrcCompStdEig<ARFLOAT>::Eupp()
{

  ceupp(this->rvec, this->HowMny, this->EigValR, this->EigVec, this->n, this->sigmaR, this->workv,
        this->bmat, this->n, this->which, this->nev, this->tol, this->resid, this->ncv, this->V, this->n, this->iparam,
        this->ipntr, this->workd, this->workl, this->lworkl, this->rwork, this->info);

} // Eupp.


template<class ARFLOAT>
MKL_INT ARrcCompStdEig<ARFLOAT>::
Eigenvalues(arcomplex<ARFLOAT>* &EigValp, bool ivec, bool ischur)
{

  if (this->ValuesOK) {                      // Eigenvalues are available .
    if (EigValp == NULL) {             // Moving eigenvalues.
      EigValp  = this->EigValR;
      this->EigValR  = NULL;
      this->newVal   = false;
      this->ValuesOK = false;
    }
    else {                             // Copying eigenvalues.
      copy(this->nconv,this->EigValR,1,EigValp,1);
    }
  }
  else {
    if (this->newVal) {
      delete[] this->EigValR;
      this->newVal = false;
    }
    if (EigValp == NULL) {
      try { EigValp = new arcomplex<ARFLOAT>[this->ValSize()]; }
      catch (ArpackError) { return 0; }
    }
    this->EigValR = EigValp;
    if (ivec) {                        // Finding eigenvalues and eigenvectors.
      this->nconv = this->FindEigenvectors(ischur);
    }
    else {                             // Finding eigenvalues only.
      this->nconv = this->FindEigenvalues();
    }
    this->EigValR = NULL;
  }
  return this->nconv;

} // Eigenvalues(EigValp, ivec, ischur).


template<class ARFLOAT>
MKL_INT ARrcCompStdEig<ARFLOAT>::
EigenValVectors(arcomplex<ARFLOAT>* &EigVecp, arcomplex<ARFLOAT>* &EigValp,
                bool ischur)
{

  if (this->ValuesOK) {                  // Eigenvalues are already available.
    this->nconv = Eigenvalues(EigValp, false);
    this->nconv = this->Eigenvectors(EigVecp, ischur);
  }
  else {                           // Eigenvalues and vectors are not available.
    if (this->newVec) {
      delete[] this->EigVec;
      this->newVec = false;
    }
    if (this->newVal) {
      delete[] this->EigValR;
      this->newVal = false;
    }  
    try {
      if (EigVecp == NULL) EigVecp = new arcomplex<ARFLOAT>[this->ValSize()*this->n];
      if (EigValp == NULL) EigValp = new arcomplex<ARFLOAT>[this->ValSize()];
    }
    catch (ArpackError) { return 0; }
    this->EigVec  = EigVecp;
    this->EigValR = EigValp;
    this->nconv   = this->FindEigenvectors(ischur);
    this->EigVec  = NULL;
    this->EigValR = NULL;
  }
  return this->nconv;

} // EigenValVectors(EigVecp, EigValp, ischur).


template<class ARFLOAT>
inline arcomplex<ARFLOAT> ARrcCompStdEig<ARFLOAT>::Eigenvalue(MKL_INT i)
// calcula e retorna um autovalor.

{

  // Returning i-eth eigenvalue.

  if (!this->ValuesOK) {
    throw ArpackError(ArpackError::VALUES_NOT_OK, "Eigenvalue(i)");
  }
  else if ((i>=this->nconv)||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvalue(i)");
  }
  return this->EigValR[i];

} // Eigenvalue(i).


template<class ARFLOAT>
inline arcomplex<ARFLOAT> ARrcCompStdEig<ARFLOAT>::
Eigenvector(MKL_INT i, MKL_INT j)
{

  // Returning element j of i-eth eigenvector.

  if (!this->VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "Eigenvector(i,j)");
  }
  else if ((i>=this->nconv)||(i<0)||(j>=this->n)||(j<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "Eigenvector(i,j)");
  }
  return this->EigVec[i*this->n+j];

} // Eigenvector(i,j).


#ifdef STL_VECTOR_H

template<class ARFLOAT>
inline vector<arcomplex<ARFLOAT> >* ARrcCompStdEig<ARFLOAT>::
StlEigenvalues(bool ivec, bool ischur)
{

  // Returning the eigenvalues in a STL vector.

  vector<arcomplex<ARFLOAT> >* ValR;
  arcomplex<ARFLOAT>*          ValPtr;

  try {
    ValR = new vector<arcomplex<ARFLOAT> >(ValSize());
  }
  catch (ArpackError) { return NULL; }
  ValPtr = ValR->begin();
  nconv = Eigenvalues(ValPtr, ivec, ischur);
  return ValR;

} // StlEigenvalues.


template<class ARFLOAT>
inline vector<arcomplex<ARFLOAT> >* ARrcCompStdEig<ARFLOAT>::
StlEigenvector(MKL_INT i)
{

  // Returning the i-th eigenvector in a STL vector.

  vector<arcomplex<ARFLOAT> >* Vec;

  if (!VectorsOK) {
    throw ArpackError(ArpackError::VECTORS_NOT_OK, "StlEigenvector(i)");
  }
  else if ((i>=ValSize())||(i<0)) {
    throw ArpackError(ArpackError::RANGE_ERROR, "StlEigenvector(i)");
  }
  try {
    Vec = new vector<arcomplex<ARFLOAT> >(&EigVec[i*n], &EigVec[(i+1)*n]);
  }
  catch (ArpackError) { return NULL; }
  return Vec;

} // StlEigenvector(i).

#endif // #ifdef STL_VECTOR_H.


template<class ARFLOAT>
inline ARrcCompStdEig<ARFLOAT>::
ARrcCompStdEig(MKL_INT np, MKL_INT nevp, const std::string& whichp, MKL_INT ncvp, ARFLOAT tolp,
               MKL_INT maxitp, arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->NoShift();
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARrcCompStdEig<ARFLOAT>::
ARrcCompStdEig(MKL_INT np, MKL_INT nevp, arcomplex<ARFLOAT> sigmap,
               const std::string& whichp, MKL_INT ncvp, ARFLOAT tolp, MKL_INT maxitp,
               arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->ChangeShift(sigmap);
  this->DefineParameters(np, nevp, whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARrcCompStdEig<ARFLOAT>& ARrcCompStdEig<ARFLOAT>::
operator=(const ARrcCompStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARRSCOMP_H


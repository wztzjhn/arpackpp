/*
   ARPACK++ v1.2 2/20/2000
   c++ MKL_INTerface to ARPACK code.

   MODULE ARUSComp.h.
   Arpack++ class ARluCompStdEig definition
   (umfpack version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUSCOMP_H
#define ARUSCOMP_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arscomp.h"
#include "arunsmat.h"
#include "arrseig.h"


template<class ARFLOAT>
class ARluCompStdEig:
  public virtual ARCompStdEig<ARFLOAT, 
                              ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> > {

 public:

 // a) Public functions:

 // a.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(arcomplex<ARFLOAT> sigmaRp);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(arcomplex<ARFLOAT> sigmap);

 // a.2) Constructors and destructor.

  ARluCompStdEig() { }
  // Short constructor.

  ARluCompStdEig(MKL_INT nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A, 
                 const std::string& whichp = "LM", MKL_INT ncvp = 0,
                 ARFLOAT tolp = 0.0, MKL_INT maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluCompStdEig(MKL_INT nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
                 arcomplex<ARFLOAT> sigma, const std::string& whichp = "LM",
                 MKL_INT ncvp = 0, ARFLOAT tolp = 0.0, MKL_INT maxitp = 0,
                 arcomplex<ARFLOAT>* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert mode).

  ARluCompStdEig(const ARluCompStdEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluCompStdEig() { }
  // Destructor.


 // b) Operators.

  ARluCompStdEig& operator=(const ARluCompStdEig& other);
  // Assignment operator.

}; // class ARluCompStdEig.


// ------------------------------------------------------------------------ //
// ARluCompStdEig member functions definition.                              //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARluCompStdEig<ARFLOAT>::
ChangeShift(arcomplex<ARFLOAT> sigmaRp)
{

  this->objOP->FactorAsI(sigmaRp);
  ARrcStdEig<ARFLOAT, arcomplex<ARFLOAT> >::ChangeShift(sigmaRp);

} // ChangeShift.


template<class ARFLOAT>
inline void ARluCompStdEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>, 
           ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetRegularMode(this->objOP, 
                   &ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultMv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARluCompStdEig<ARFLOAT>::
SetShiftInvertMode(arcomplex<ARFLOAT> sigmap)
{

  ARStdEig<ARFLOAT, arcomplex<ARFLOAT>, 
           ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT> >::
    SetShiftInvertMode(sigmap, this->objOP,
                       &ARumNonSymMatrix<arcomplex<ARFLOAT>,ARFLOAT>::MultInvv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline ARluCompStdEig<ARFLOAT>::
ARluCompStdEig(MKL_INT nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               const std::string& whichp, MKL_INT ncvp, ARFLOAT tolp,
               MKL_INT maxitp, arcomplex<ARFLOAT>* residp, bool ishiftp)

{

  this->NoShift();
  DefineParameters(A.ncols(), nevp, &A,
                   &ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultMv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARluCompStdEig<ARFLOAT>::
ARluCompStdEig(MKL_INT nevp, ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>& A,
               arcomplex<ARFLOAT> sigmap, const std::string& whichp, MKL_INT ncvp,
               ARFLOAT tolp, MKL_INT maxitp, arcomplex<ARFLOAT>* residp,
               bool ishiftp)

{

  DefineParameters(A.ncols(), nevp, &A, 
                   &ARumNonSymMatrix<arcomplex<ARFLOAT>, ARFLOAT>::MultInvv,
                   whichp, ncvp, tolp, maxitp, residp, ishiftp);
  ChangeShift(sigmap);

} // Long constructor (shift and invert mode).


template<class ARFLOAT>
ARluCompStdEig<ARFLOAT>& ARluCompStdEig<ARFLOAT>::
operator=(const ARluCompStdEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
    this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUSCOMP_H

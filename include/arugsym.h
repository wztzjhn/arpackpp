/*
   ARPACK++ v1.2 2/20/2000
   c++ MKL_INTerface to ARPACK code.

   MODULE ARUGSym.h.
   Arpack++ class ARluSymGenEig definition
   (UMFPACK version).

   ARPACK Authors
      Richard Lehoucq
      Danny Sorensen
      Chao Yang
      Kristi Maschhoff
      Dept. of Computational & Applied Mathematics
      Rice University
      Houston, Texas
*/

#ifndef ARUGSYM_H
#define ARUGSYM_H

#include <cstddef>
#include <string>
#include "arch.h"
#include "arusmat.h"
#include "aruspen.h"
#include "argsym.h"


template<class ARFLOAT>
class ARluSymGenEig:
  public virtual ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>,
                             ARumSymPencil<ARFLOAT> > {

 private:

 // a) Data structure used to store matrices.

  ARumSymPencil<ARFLOAT> Pencil;

 // b) Protected functions:

  virtual void Copy(const ARluSymGenEig& other);
  // Makes a deep copy of "other" over "this" object.
  // Old values are not deleted (this function is to be used
  // by the copy constructor and the assignment operator only).


 public:

 // c) Public functions:

 // c.1) Functions that allow changes in problem parameters.

  virtual void ChangeShift(ARFLOAT sigmap);

  virtual void SetRegularMode();

  virtual void SetShiftInvertMode(ARFLOAT sigmap);

  virtual void SetBucklingMode(ARFLOAT sigmap);

  virtual void SetCayleyMode(ARFLOAT sigmap);

 // c.2) Constructors and destructor.

  ARluSymGenEig() { }
  // Short constructor.

  ARluSymGenEig(MKL_INT nevp, ARumSymMatrix<ARFLOAT>& A,
                ARumSymMatrix<ARFLOAT>& B, const std::string& whichp = "LM",
                MKL_INT ncvp = 0, ARFLOAT tolp = 0.0, MKL_INT maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (regular mode).

  ARluSymGenEig(char InvertModep, MKL_INT nevp, ARumSymMatrix<ARFLOAT>& A,
                ARumSymMatrix<ARFLOAT>& B, ARFLOAT sigma, const std::string& whichp = "LM", 
                MKL_INT ncvp = 0, ARFLOAT tolp = 0.0, MKL_INT maxitp = 0,
                ARFLOAT* residp = NULL, bool ishiftp = true);
  // Long constructor (shift and invert, buckling and Cayley modes).

  ARluSymGenEig(const ARluSymGenEig& other) { Copy(other); }
  // Copy constructor.

  virtual ~ARluSymGenEig() { }
  // Destructor.

 // d) Operators.

  ARluSymGenEig& operator=(const ARluSymGenEig& other);
  // Assignment operator.

}; // class ARluSymGenEig.


// ------------------------------------------------------------------------ //
// ARluSymGenEig member functions definition.                               //
// ------------------------------------------------------------------------ //


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::
Copy(const ARluSymGenEig<ARFLOAT>& other)
{

  ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>,
              ARumSymPencil<ARFLOAT> >:: Copy(other);
  Pencil = other.Pencil;
  this->objOP  = &Pencil;
  this->objB   = &Pencil;
  this->objA   = &Pencil;

} // Copy.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::ChangeShift(ARFLOAT sigmap)
{

  this->objOP->FactorAsB(sigmap);
  ARrcSymGenEig<ARFLOAT>::ChangeShift(sigmap);

} // ChangeShift.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::SetRegularMode()
{

  ARStdEig<ARFLOAT, ARFLOAT, ARumSymPencil<ARFLOAT> >::
    SetRegularMode(&Pencil, &ARumSymPencil<ARFLOAT>::MultInvBAv);

} // SetRegularMode.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::
SetShiftInvertMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>, ARumSymPencil<ARFLOAT> >::
    SetShiftInvertMode(sigmap, &Pencil, &ARumSymPencil<ARFLOAT>::MultInvAsBv);
  this->ChangeMultBx(&Pencil, &ARumSymPencil<ARFLOAT>::MultBv);

} // SetShiftInvertMode.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::
SetBucklingMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>, ARumSymPencil<ARFLOAT> >::
    SetBucklingMode(sigmap, &Pencil, &ARumSymPencil<ARFLOAT>::MultInvAsBv);
  this->ChangeMultBx(&Pencil, &ARumSymPencil<ARFLOAT>::MultAv);

} // SetBucklingMode.


template<class ARFLOAT>
inline void ARluSymGenEig<ARFLOAT>::
SetCayleyMode(ARFLOAT sigmap)
{

  ARSymGenEig<ARFLOAT, ARumSymPencil<ARFLOAT>, ARumSymPencil<ARFLOAT> >::
    SetCayleyMode(sigmap, &Pencil, &ARumSymPencil<ARFLOAT>::MultInvAsBv,
                  &Pencil, &ARumSymPencil<ARFLOAT>::MultAv);
  this->ChangeMultBx(&Pencil, &ARumSymPencil<ARFLOAT>::MultBv);

} // SetCayleyMode.


template<class ARFLOAT>
inline ARluSymGenEig<ARFLOAT>::
ARluSymGenEig(MKL_INT nevp, ARumSymMatrix<ARFLOAT>& A,
              ARumSymMatrix<ARFLOAT>& B, const std::string& whichp, MKL_INT ncvp,
              ARFLOAT tolp, MKL_INT maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->InvertMode = 'S';
  this->NoShift();
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARumSymPencil<ARFLOAT>::MultInvBAv, &Pencil,
                   &ARumSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);

} // Long constructor (regular mode).


template<class ARFLOAT>
inline ARluSymGenEig<ARFLOAT>::
ARluSymGenEig(char InvertModep, MKL_INT nevp, ARumSymMatrix<ARFLOAT>& A,
              ARumSymMatrix<ARFLOAT>& B, ARFLOAT sigmap,
              const std::string& whichp, MKL_INT ncvp, ARFLOAT tolp,
              MKL_INT maxitp, ARFLOAT* residp, bool ishiftp)

{

  Pencil.DefineMatrices(A, B);
  this->DefineParameters(A.ncols(), nevp, &Pencil,
                   &ARumSymPencil<ARFLOAT>::MultInvAsBv, &Pencil,
                   &ARumSymPencil<ARFLOAT>::MultBv, whichp,
                   ncvp, tolp, maxitp, residp, ishiftp);
  this->InvertMode = this->CheckInvertMode(InvertModep);
  switch (this->InvertMode) {
  case 'B':  // Buckling mode.
    this->ChangeMultBx(&Pencil, &ARumSymPencil<ARFLOAT>::MultAv);
  case 'S':  // Shift and invert mode.
    this->ChangeShift(sigmap);
    break;  
  case 'C':  // Cayley mode.
    this->SetCayleyMode(sigmap);
  }

} // Long constructor (shift and invert, buckling and Cayley modes).


template<class ARFLOAT>
ARluSymGenEig<ARFLOAT>& ARluSymGenEig<ARFLOAT>::
operator=(const ARluSymGenEig<ARFLOAT>& other)
{

  if (this != &other) { // Stroustrup suggestion.
   this->ClearMem();
    Copy(other);
  }
  return *this;

} // operator=.


#endif // ARUGSYM_H

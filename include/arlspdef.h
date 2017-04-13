/*
   ARPACK++ v1.2 2/20/2000
   c++ MKL_INTerface to ARPACK code.

   MODULE ARLSpDef.h.
   ALTERED version of slu_sdefs.h slu_ddefs.h slu_cdefs.h slu_zdefs.h 
   (from SuperLU 3.0 package).
*/


#ifndef __SUPERLU_SP_DEFS /* allow multiple inclusions */
#define __SUPERLU_SP_DEFS

/*
 * File name:		sp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */
#include "arlnames.h"
#include "arlsupm.h"
#include "arlcomp.h"
#include "arlutil.h"
#ifdef _CRAY
#include <fortran.h>
#include <string.h>
#endif

/* Define my MKL_INTeger type MKL_INT_t */
typedef MKL_INT MKL_INT_t; /* default */

// /* No of marker arrays used in the symbolic factorization,
//    each of size n */
// #define NO_MARKER     3
// #define NUM_TEMPV(m,w,t,b)  ( MAX(m, (t + b)*w) )
// 
// typedef enum {LUSUP, UCOL, LSUB, USUB} MemType;
// typedef enum {HEAD, TAIL}              stack_end_t;
// typedef enum {SYSTEM, USER}            LU_space_t;

/*
 * Global data structures used in LU factorization -
 * 
 *   nsuper: #supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new 
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub 
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript MKL_INTerchange
 *	performed, similar MKL_INTerchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numercial LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 */
 
typedef struct {
    MKL_INT     *xsup;    /* supernode and column mapping */
    MKL_INT     *supno;   
    MKL_INT     *lsub;    /* compressed L subscripts */
    MKL_INT	    *xlsub;
    float  *lusup;   /* L supernodes */
    MKL_INT     *xlusup;
    float  *ucol;    /* U columns */
    MKL_INT     *usub;
    MKL_INT	    *xusub;
    MKL_INT     nzlmax;   /* current max size of lsub */
    MKL_INT     nzumax;   /*    "    "    "      ucol */
    MKL_INT     nzlumax;  /*    "    "    "     lusup */
    MKL_INT     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
    MKL_INT     num_expansions;
    ExpHeader *expanders; /* Array of pointers to 4 types of memory */
    LU_stack_t stack;     /* use user supplied memory */
} sGlobalLU_t;

typedef struct {
    MKL_INT     *xsup;    /* supernode and column mapping */
    MKL_INT     *supno;   
    MKL_INT     *lsub;    /* compressed L subscripts */
    MKL_INT	    *xlsub;
    double  *lusup;   /* L supernodes */
    MKL_INT     *xlusup;
    double  *ucol;    /* U columns */
    MKL_INT     *usub;
    MKL_INT	    *xusub;
    MKL_INT     nzlmax;   /* current max size of lsub */
    MKL_INT     nzumax;   /*    "    "    "      ucol */
    MKL_INT     nzlumax;  /*    "    "    "     lusup */
    MKL_INT     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
    MKL_INT     num_expansions;
    ExpHeader *expanders; /* Array of pointers to 4 types of memory */
    LU_stack_t stack;     /* use user supplied memory */
} dGlobalLU_t;

typedef struct {
    MKL_INT     *xsup;    /* supernode and column mapping */
    MKL_INT     *supno;   
    MKL_INT     *lsub;    /* compressed L subscripts */
    MKL_INT	    *xlsub;
    lscomplex  *lusup;   /* L supernodes */
    MKL_INT     *xlusup;
    lscomplex  *ucol;    /* U columns */
    MKL_INT     *usub;
    MKL_INT	    *xusub;
    MKL_INT     nzlmax;   /* current max size of lsub */
    MKL_INT     nzumax;   /*    "    "    "      ucol */
    MKL_INT     nzlumax;  /*    "    "    "     lusup */
    MKL_INT     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
    MKL_INT     num_expansions;
    ExpHeader *expanders; /* Array of pointers to 4 types of memory */
    LU_stack_t stack;     /* use user supplied memory */
} cGlobalLU_t;

typedef struct {
    MKL_INT     *xsup;    /* supernode and column mapping */
    MKL_INT     *supno;   
    MKL_INT     *lsub;    /* compressed L subscripts */
    MKL_INT	    *xlsub;
    ldcomplex  *lusup;   /* L supernodes */
    MKL_INT     *xlusup;
    ldcomplex  *ucol;    /* U columns */
    MKL_INT     *usub;
    MKL_INT	    *xusub;
    MKL_INT     nzlmax;   /* current max size of lsub */
    MKL_INT     nzumax;   /*    "    "    "      ucol */
    MKL_INT     nzlumax;  /*    "    "    "     lusup */
    MKL_INT     n;        /* number of columns in the matrix */
    LU_space_t MemModel; /* 0 - system malloc'd; 1 - user provided */
    MKL_INT     num_expansions;
    ExpHeader *expanders; /* Array of pointers to 4 types of memory */
    LU_stack_t stack;     /* use user supplied memory */
} zGlobalLU_t;

// typedef struct {
//     MKL_INT panel_size;
//     MKL_INT relax;
//     float diag_pivot_thresh;
//     float drop_tol;
// } sfactor_param_t;
// 
// typedef struct {
//     MKL_INT panel_size;
//     MKL_INT relax;
//     double diag_pivot_thresh;
//     double drop_tol;
// } dfactor_param_t;
//
//typedef struct {
//    float for_lu;
//    float total_needed;
//    MKL_INT   expansions;
//} mem_usage_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Driver routines */
extern void
sgssv(superlu_options_t *, SuperMatrix *, MKL_INT *, MKL_INT *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, MKL_INT *);
extern void
dgssv(superlu_options_t *, SuperMatrix *, MKL_INT *, MKL_INT *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, MKL_INT *);
extern void
cgssv(superlu_options_t *, SuperMatrix *, MKL_INT *, MKL_INT *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, MKL_INT *);
extern void
zgssv(superlu_options_t *, SuperMatrix *, MKL_INT *, MKL_INT *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, MKL_INT *);
extern void
sgssvx(superlu_options_t *, SuperMatrix *, MKL_INT *, MKL_INT *, MKL_INT *,
       char *, float *, float *, SuperMatrix *, SuperMatrix *,
       void *, MKL_INT, SuperMatrix *, SuperMatrix *,
       float *, float *, float *, float *,
       sGlobalLU_t *, mem_usage_t *, SuperLUStat_t *, MKL_INT *);
extern void
dgssvx(superlu_options_t *, SuperMatrix *, MKL_INT *, MKL_INT *, MKL_INT *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, MKL_INT, SuperMatrix *, SuperMatrix *,
       double *, double *, double *, double *,
       dGlobalLU_t *, mem_usage_t *, SuperLUStat_t *, MKL_INT *);
extern void
cgssvx(superlu_options_t *, SuperMatrix *, MKL_INT *, MKL_INT *, MKL_INT *,
       char *, float *, float *, SuperMatrix *, SuperMatrix *,
       void *, MKL_INT, SuperMatrix *, SuperMatrix *,
       float *, float *, float *, float *,
       cGlobalLU_t *, mem_usage_t *, SuperLUStat_t *, MKL_INT *);
extern void
zgssvx(superlu_options_t *, SuperMatrix *, MKL_INT *, MKL_INT *, MKL_INT *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, MKL_INT, SuperMatrix *, SuperMatrix *,
       double *, double *, double *, double *,
       zGlobalLU_t *, mem_usage_t *, SuperLUStat_t *, MKL_INT *);

/* Supernodal LU factor related */
extern void
sCreate_CompCol_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, float *,
		       MKL_INT *, MKL_INT *, Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_CompCol_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, double *,
		       MKL_INT *, MKL_INT *, Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_CompCol_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, lscomplex *,
		       MKL_INT *, MKL_INT *, Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_CompCol_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, ldcomplex *,
		       MKL_INT *, MKL_INT *, Stype_t, Dtype_t, Mtype_t);
extern void
sCreate_CompRow_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, float *,
		       MKL_INT *, MKL_INT *, Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_CompRow_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, double *,
		       MKL_INT *, MKL_INT *, Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_CompRow_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, lscomplex *,
		       MKL_INT *, MKL_INT *, Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_CompRow_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, ldcomplex *,
		       MKL_INT *, MKL_INT *, Stype_t, Dtype_t, Mtype_t);
extern void
sCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
dCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
cCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
zCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
sCreate_Dense_Matrix(SuperMatrix *, MKL_INT, MKL_INT, float *, MKL_INT,
		     Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_Dense_Matrix(SuperMatrix *, MKL_INT, MKL_INT, double *, MKL_INT,
		     Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_Dense_Matrix(SuperMatrix *, MKL_INT, MKL_INT, lscomplex *, MKL_INT,
		     Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_Dense_Matrix(SuperMatrix *, MKL_INT, MKL_INT, ldcomplex *, MKL_INT,
		     Stype_t, Dtype_t, Mtype_t);
extern void
sCreate_SuperNode_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, float *, 
		         MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_SuperNode_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, double *, 
		         MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_SuperNode_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, lscomplex *, 
		         MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
zCreate_SuperNode_Matrix(SuperMatrix *, MKL_INT, MKL_INT, MKL_INT, ldcomplex *, 
		         MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
sCopy_Dense_Matrix(MKL_INT, MKL_INT, float *, MKL_INT, float *, MKL_INT);
extern void
dCopy_Dense_Matrix(MKL_INT, MKL_INT, double *, MKL_INT, double *, MKL_INT);
extern void
cCopy_Dense_Matrix(MKL_INT, MKL_INT, lscomplex *, MKL_INT, lscomplex *, MKL_INT);
extern void
zCopy_Dense_Matrix(MKL_INT, MKL_INT, ldcomplex *, MKL_INT, ldcomplex *, MKL_INT);

// extern void    Destroy_SuperMatrix_Store(SuperMatrix *);
// extern void    Destroy_CompCol_Matrix(SuperMatrix *);
// extern void    Destroy_SuperNode_Matrix(SuperMatrix *);
// extern void    Destroy_CompCol_Permuted(SuperMatrix *);
// extern void    Destroy_Dense_Matrix(SuperMatrix *);
// extern void    get_perm_c(MKL_INT, SuperMatrix *, MKL_INT *);  
// extern void    sp_preorder (char*, SuperMatrix*, MKL_INT*, MKL_INT*, SuperMatrix*);
// //  extern void    countnz (const MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *, sGlobalLU_t *);
// //  extern void    fixupL (const MKL_INT, const MKL_INT *, sGlobalLU_t *);

extern void    sallocateA (MKL_INT, MKL_INT, float **, MKL_INT **, MKL_INT **);
extern void    dallocateA (MKL_INT, MKL_INT, double **, MKL_INT **, MKL_INT **);
extern void    callocateA (MKL_INT, MKL_INT, lscomplex **, MKL_INT **, MKL_INT **);
extern void    zallocateA (MKL_INT, MKL_INT, ldcomplex **, MKL_INT **, MKL_INT **);
extern void    sgstrf (superlu_options_t*, SuperMatrix*, 
                       MKL_INT, MKL_INT, MKL_INT*, void *, MKL_INT, MKL_INT *, MKL_INT *, 
                       SuperMatrix *, SuperMatrix *, sGlobalLU_t *, SuperLUStat_t*, MKL_INT *);
extern void    dgstrf (superlu_options_t*, SuperMatrix*,
                       MKL_INT, MKL_INT, MKL_INT*, void *, MKL_INT, MKL_INT *, MKL_INT *, 
                       SuperMatrix *, SuperMatrix *, dGlobalLU_t *, SuperLUStat_t*, MKL_INT *);
extern void    cgstrf (superlu_options_t*, SuperMatrix*,
                       MKL_INT, MKL_INT, MKL_INT*, void *, MKL_INT, MKL_INT *, MKL_INT *, 
                       SuperMatrix *, SuperMatrix *, cGlobalLU_t *, SuperLUStat_t*, MKL_INT *);
extern void    zgstrf (superlu_options_t*, SuperMatrix*,
                       MKL_INT, MKL_INT, MKL_INT*, void *, MKL_INT, MKL_INT *, MKL_INT *, 
                       SuperMatrix *, SuperMatrix *, zGlobalLU_t *, SuperLUStat_t*, MKL_INT *);
extern MKL_INT     ssnode_dfs (const MKL_INT, const MKL_INT, const MKL_INT *, const MKL_INT *,
			     const MKL_INT *, MKL_INT *, MKL_INT *, sGlobalLU_t *);
extern MKL_INT     dsnode_dfs (const MKL_INT, const MKL_INT, const MKL_INT *, const MKL_INT *,
			     const MKL_INT *, MKL_INT *, MKL_INT *, dGlobalLU_t *);
extern MKL_INT     csnode_dfs (const MKL_INT, const MKL_INT, const MKL_INT *, const MKL_INT *,
			     const MKL_INT *, MKL_INT *, MKL_INT *, cGlobalLU_t *);
extern MKL_INT     zsnode_dfs (const MKL_INT, const MKL_INT, const MKL_INT *, const MKL_INT *,
			     const MKL_INT *, MKL_INT *, MKL_INT *, zGlobalLU_t *);
extern MKL_INT     ssnode_bmod (const MKL_INT, const MKL_INT, const MKL_INT, float *,
                              float *, sGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     dsnode_bmod (const MKL_INT, const MKL_INT, const MKL_INT, double *,
                              double *, dGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     csnode_bmod (const MKL_INT, const MKL_INT, const MKL_INT, lscomplex *,
                              lscomplex *, cGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     zsnode_bmod (const MKL_INT, const MKL_INT, const MKL_INT, ldcomplex *,
                              ldcomplex *, zGlobalLU_t *, SuperLUStat_t*);
extern void    spanel_dfs (const MKL_INT, const MKL_INT, const MKL_INT, SuperMatrix *,
			   MKL_INT *, MKL_INT *, float *, MKL_INT *, MKL_INT *, MKL_INT *,
			   MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, sGlobalLU_t *);
extern void    dpanel_dfs (const MKL_INT, const MKL_INT, const MKL_INT, SuperMatrix *,
			   MKL_INT *, MKL_INT *, double *, MKL_INT *, MKL_INT *, MKL_INT *,
			   MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, dGlobalLU_t *);
extern void    cpanel_dfs (const MKL_INT, const MKL_INT, const MKL_INT, SuperMatrix *,
			   MKL_INT *, MKL_INT *, lscomplex *, MKL_INT *, MKL_INT *, MKL_INT *,
			   MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, cGlobalLU_t *);
extern void    zpanel_dfs (const MKL_INT, const MKL_INT, const MKL_INT, SuperMatrix *,
			   MKL_INT *, MKL_INT *, ldcomplex *, MKL_INT *, MKL_INT *, MKL_INT *,
			   MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, zGlobalLU_t *);
extern void    spanel_bmod (const MKL_INT, const MKL_INT, const MKL_INT, const MKL_INT,
                           float *, float *, MKL_INT *, MKL_INT *,
			   sGlobalLU_t *, SuperLUStat_t*);
extern void    dpanel_bmod (const MKL_INT, const MKL_INT, const MKL_INT, const MKL_INT,
                           double *, double *, MKL_INT *, MKL_INT *,
			   dGlobalLU_t *, SuperLUStat_t*);
extern void    cpanel_bmod (const MKL_INT, const MKL_INT, const MKL_INT, const MKL_INT,
                           lscomplex *, lscomplex *, MKL_INT *, MKL_INT *,
			   cGlobalLU_t *, SuperLUStat_t*);
extern void    zpanel_bmod (const MKL_INT, const MKL_INT, const MKL_INT, const MKL_INT,
                           ldcomplex *, ldcomplex *, MKL_INT *, MKL_INT *,
			   zGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     scolumn_dfs (const MKL_INT, const MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
			   MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, sGlobalLU_t *);
extern MKL_INT     dcolumn_dfs (const MKL_INT, const MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
			   MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, dGlobalLU_t *);
extern MKL_INT     ccolumn_dfs (const MKL_INT, const MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
			   MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, cGlobalLU_t *);
extern MKL_INT     zcolumn_dfs (const MKL_INT, const MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *,
			   MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, zGlobalLU_t *);
extern MKL_INT     scolumn_bmod (const MKL_INT, const MKL_INT, float *,
			   float *, MKL_INT *, MKL_INT *, MKL_INT, sGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     dcolumn_bmod (const MKL_INT, const MKL_INT, double *,
			   double *, MKL_INT *, MKL_INT *, MKL_INT, dGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     ccolumn_bmod (const MKL_INT, const MKL_INT, lscomplex *,
			   lscomplex *, MKL_INT *, MKL_INT *, MKL_INT, cGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     zcolumn_bmod (const MKL_INT, const MKL_INT, ldcomplex *,
			   ldcomplex *, MKL_INT *, MKL_INT *, MKL_INT, zGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     scopy_to_ucol (MKL_INT, MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *,
                              float *, sGlobalLU_t *);         
extern MKL_INT     dcopy_to_ucol (MKL_INT, MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *,
                              double *, dGlobalLU_t *);         
extern MKL_INT     ccopy_to_ucol (MKL_INT, MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *,
                              lscomplex *, cGlobalLU_t *);         
extern MKL_INT     zcopy_to_ucol (MKL_INT, MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *,
                              ldcomplex *, zGlobalLU_t *);         
extern MKL_INT     spivotL (const MKL_INT, const float, MKL_INT *, MKL_INT *, 
                              MKL_INT *, MKL_INT *, MKL_INT *, sGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     dpivotL (const MKL_INT, const double, MKL_INT *, MKL_INT *, 
                              MKL_INT *, MKL_INT *, MKL_INT *, dGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     cpivotL (const MKL_INT, const float, MKL_INT *, MKL_INT *, 
                              MKL_INT *, MKL_INT *, MKL_INT *, cGlobalLU_t *, SuperLUStat_t*);
extern MKL_INT     zpivotL (const MKL_INT, const double, MKL_INT *, MKL_INT *, 
                              MKL_INT *, MKL_INT *, MKL_INT *, zGlobalLU_t *, SuperLUStat_t*);
extern void    spruneL (const MKL_INT, const MKL_INT *, const MKL_INT, const MKL_INT,
			     const MKL_INT *, const MKL_INT *, MKL_INT *, sGlobalLU_t *);
extern void    dpruneL (const MKL_INT, const MKL_INT *, const MKL_INT, const MKL_INT,
			     const MKL_INT *, const MKL_INT *, MKL_INT *, dGlobalLU_t *);
extern void    cpruneL (const MKL_INT, const MKL_INT *, const MKL_INT, const MKL_INT,
			     const MKL_INT *, const MKL_INT *, MKL_INT *, cGlobalLU_t *);
extern void    zpruneL (const MKL_INT, const MKL_INT *, const MKL_INT, const MKL_INT,
			     const MKL_INT *, const MKL_INT *, MKL_INT *, zGlobalLU_t *);
extern void    sreadmt (MKL_INT *, MKL_INT *, MKL_INT *, float **, MKL_INT **, MKL_INT **);
extern void    dreadmt (MKL_INT *, MKL_INT *, MKL_INT *, double **, MKL_INT **, MKL_INT **);
extern void    creadmt (MKL_INT *, MKL_INT *, MKL_INT *, lscomplex **, MKL_INT **, MKL_INT **);
extern void    zreadmt (MKL_INT *, MKL_INT *, MKL_INT *, ldcomplex **, MKL_INT **, MKL_INT **);
extern void    sGenXtrue (MKL_INT, MKL_INT, float *, MKL_INT);
extern void    dGenXtrue (MKL_INT, MKL_INT, double *, MKL_INT);
extern void    cGenXtrue (MKL_INT, MKL_INT, lscomplex *, MKL_INT);
extern void    zGenXtrue (MKL_INT, MKL_INT, ldcomplex *, MKL_INT);
extern void    sFillRHS (trans_t, MKL_INT, float *, MKL_INT, SuperMatrix *,
			SuperMatrix *);
extern void    dFillRHS (trans_t, MKL_INT, double *, MKL_INT, SuperMatrix *,
			SuperMatrix *);
extern void    cFillRHS (trans_t, MKL_INT, lscomplex *, MKL_INT, SuperMatrix *,
			SuperMatrix *);
extern void    zFillRHS (trans_t, MKL_INT, ldcomplex *, MKL_INT, SuperMatrix *,
			SuperMatrix *);
extern void    sgstrs (trans_t, SuperMatrix *, SuperMatrix *, MKL_INT *, MKL_INT *,
                        SuperMatrix *, SuperLUStat_t*, MKL_INT *);
extern void    dgstrs (trans_t, SuperMatrix *, SuperMatrix *, MKL_INT *, MKL_INT *,
                        SuperMatrix *, SuperLUStat_t*, MKL_INT *);
extern void    cgstrs (trans_t, SuperMatrix *, SuperMatrix *, MKL_INT *, MKL_INT *,
                        SuperMatrix *, SuperLUStat_t*, MKL_INT *);
extern void    zgstrs (trans_t, SuperMatrix *, SuperMatrix *, MKL_INT *, MKL_INT *,
                        SuperMatrix *, SuperLUStat_t*, MKL_INT *);


/* Driver related */

extern void    sgsequ (SuperMatrix *, float *, float *, float *,
			     float *, float *, MKL_INT *);
extern void    dgsequ (SuperMatrix *, double *, double *, double *,
			     double *, double *, MKL_INT *);
extern void    cgsequ (SuperMatrix *, float *, float *, float *,
			     float *, float *, MKL_INT *);
extern void    zgsequ (SuperMatrix *, double *, double *, double *,
			     double *, double *, MKL_INT *);
extern void    slaqgs (SuperMatrix *, float *, float *, float,
                             float, float, char *);
extern void    dlaqgs (SuperMatrix *, double *, double *, double,
                             double, double, char *);
extern void    claqgs (SuperMatrix *, float *, float *, float,
                             float, float, char *);
extern void    zlaqgs (SuperMatrix *, double *, double *, double,
                             double, double, char *);
extern void    sgscon (char *, SuperMatrix *, SuperMatrix *, 
		         float, float *, SuperLUStat_t*, MKL_INT *);
extern void    dgscon (char *, SuperMatrix *, SuperMatrix *, 
		         double, double *, SuperLUStat_t*, MKL_INT *);
extern void    cgscon (char *, SuperMatrix *, SuperMatrix *, 
		         float, float *, SuperLUStat_t*, MKL_INT *);
extern void    zgscon (char *, SuperMatrix *, SuperMatrix *, 
		         double, double *, SuperLUStat_t*, MKL_INT *);

extern float   sPivotGrowth(MKL_INT, SuperMatrix *, MKL_INT *, 
                            SuperMatrix *, SuperMatrix *);
extern double  dPivotGrowth(MKL_INT, SuperMatrix *, MKL_INT *, 
                            SuperMatrix *, SuperMatrix *);
extern float   cPivotGrowth(MKL_INT, SuperMatrix *, MKL_INT *, 
                            SuperMatrix *, SuperMatrix *);
extern double  zPivotGrowth(MKL_INT, SuperMatrix *, MKL_INT *, 
                            SuperMatrix *, SuperMatrix *);
extern void    sgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, MKL_INT *, MKL_INT *, char *, float *, 
                       float *, SuperMatrix *, SuperMatrix *,
                       float *, float *, SuperLUStat_t*, MKL_INT *);
extern void    dgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, MKL_INT *, MKL_INT *, char *, double *, 
                       double *, SuperMatrix *, SuperMatrix *,
                       double *, double *, SuperLUStat_t*, MKL_INT *);
extern void    cgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, MKL_INT *, MKL_INT *, char *, float *, 
                       float *, SuperMatrix *, SuperMatrix *,
                       float *, float *, SuperLUStat_t*, MKL_INT *);
extern void    zgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, MKL_INT *, MKL_INT *, char *, double *, 
                       double *, SuperMatrix *, SuperMatrix *,
                       double *, double *, SuperLUStat_t*, MKL_INT *);

extern MKL_INT     sp_strsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, float *, SuperLUStat_t*, MKL_INT *);
extern MKL_INT     sp_dtrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, double *, SuperLUStat_t*, MKL_INT *);
extern MKL_INT     sp_ctrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, lscomplex *, SuperLUStat_t*, MKL_INT *);
extern MKL_INT     sp_ztrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, ldcomplex *, SuperLUStat_t*, MKL_INT *);
extern MKL_INT     sp_sgemv (char *, float, SuperMatrix *, float *,
			int, float, float *, MKL_INT);
extern MKL_INT     sp_dgemv (char *, double, SuperMatrix *, double *,
			int, double, double *, MKL_INT);
extern MKL_INT     sp_cgemv (char *, lscomplex, SuperMatrix *, lscomplex *,
			int, lscomplex, lscomplex *, MKL_INT);
extern MKL_INT     sp_zgemv (char *, ldcomplex, SuperMatrix *, ldcomplex *,
			int, ldcomplex, ldcomplex *, MKL_INT);

extern MKL_INT     sp_sgemm (char *, char *, MKL_INT, MKL_INT, MKL_INT, float,
			SuperMatrix *, float *, MKL_INT, float, 
			float *, MKL_INT);
extern MKL_INT     sp_dgemm (char *, char *, MKL_INT, MKL_INT, MKL_INT, double,
			SuperMatrix *, double *, MKL_INT, double, 
			double *, MKL_INT);
extern MKL_INT     sp_cgemm (char *, char *, MKL_INT, MKL_INT, MKL_INT, lscomplex,
			SuperMatrix *, lscomplex *, MKL_INT, lscomplex, 
			lscomplex *, MKL_INT);
extern MKL_INT     sp_zgemm (char *, char *, MKL_INT, MKL_INT, MKL_INT, ldcomplex,
			SuperMatrix *, ldcomplex *, MKL_INT, ldcomplex, 
			ldcomplex *, MKL_INT);

/* Memory-related */
extern MKL_INT     sLUMemInit (fact_t, void *, MKL_INT, MKL_INT, MKL_INT, MKL_INT, MKL_INT,
			     float, SuperMatrix *, SuperMatrix *,
			     sGlobalLU_t *, MKL_INT **, float **);
extern MKL_INT     dLUMemInit (fact_t, void *, MKL_INT, MKL_INT, MKL_INT, MKL_INT, MKL_INT,
			     double, SuperMatrix *, SuperMatrix *,
			     dGlobalLU_t *, MKL_INT **, double **);
extern MKL_INT     cLUMemInit (fact_t, void *, MKL_INT, MKL_INT, MKL_INT, MKL_INT, MKL_INT,
			     float, SuperMatrix *, SuperMatrix *,
			     cGlobalLU_t *, MKL_INT **, lscomplex **);
extern MKL_INT     zLUMemInit (fact_t, void *, MKL_INT, MKL_INT, MKL_INT, MKL_INT, MKL_INT,
			     double, SuperMatrix *, SuperMatrix *,
			     zGlobalLU_t *, MKL_INT **, ldcomplex **);
extern void    sSetRWork (MKL_INT, MKL_INT, float *, float **, float **);
extern void    dSetRWork (MKL_INT, MKL_INT, double *, double **, double **);
extern void    cSetRWork (MKL_INT, MKL_INT, lscomplex *, lscomplex **, lscomplex **);
extern void    zSetRWork (MKL_INT, MKL_INT, ldcomplex *, ldcomplex **, ldcomplex **);
extern void    sLUWorkFree (MKL_INT *, float *, sGlobalLU_t *);
extern void    dLUWorkFree (MKL_INT *, double *, dGlobalLU_t *);
extern void    cLUWorkFree (MKL_INT *, lscomplex *, cGlobalLU_t *);
extern void    zLUWorkFree (MKL_INT *, ldcomplex *, zGlobalLU_t *);
extern MKL_INT     sLUMemXpand (MKL_INT, MKL_INT, MemType, MKL_INT *, sGlobalLU_t *);
extern MKL_INT     dLUMemXpand (MKL_INT, MKL_INT, MemType, MKL_INT *, dGlobalLU_t *);
extern MKL_INT     cLUMemXpand (MKL_INT, MKL_INT, MemType, MKL_INT *, cGlobalLU_t *);
extern MKL_INT     zLUMemXpand (MKL_INT, MKL_INT, MemType, MKL_INT *, zGlobalLU_t *);

extern float  *floatMalloc(MKL_INT);
extern double  *doubleMalloc(MKL_INT);
extern lscomplex  *complexMalloc(MKL_INT);
extern ldcomplex  *doublecomplexMalloc(MKL_INT);
extern float  *floatCalloc(MKL_INT);
extern double  *doubleCalloc(MKL_INT);
extern lscomplex  *complexCalloc(MKL_INT);
extern ldcomplex  *doublecomplexCalloc(MKL_INT);
extern MKL_INT     smemory_usage(const MKL_INT, const MKL_INT, const MKL_INT, const MKL_INT);
extern MKL_INT     dmemory_usage(const MKL_INT, const MKL_INT, const MKL_INT, const MKL_INT);
extern MKL_INT     cmemory_usage(const MKL_INT, const MKL_INT, const MKL_INT, const MKL_INT);
extern MKL_INT     zmemory_usage(const MKL_INT, const MKL_INT, const MKL_INT, const MKL_INT);
extern MKL_INT     sQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern MKL_INT     dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern MKL_INT     cQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
extern MKL_INT     zQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/* Auxiliary routines */
extern void    sreadhb(FILE *, MKL_INT *, MKL_INT *, MKL_INT *, float **, MKL_INT **, MKL_INT **);
extern void    dreadhb(FILE *, MKL_INT *, MKL_INT *, MKL_INT *, double **, MKL_INT **, MKL_INT **);
extern void    creadhb(FILE *, MKL_INT *, MKL_INT *, MKL_INT *, lscomplex **, MKL_INT **, MKL_INT **);
extern void    zreadhb(FILE *, MKL_INT *, MKL_INT *, MKL_INT *, ldcomplex **, MKL_INT **, MKL_INT **);
extern void    sCompRow_to_CompCol(MKL_INT, MKL_INT, MKL_INT, float*, MKL_INT*, MKL_INT*,
		                   float **, MKL_INT **, MKL_INT **);
extern void    dCompRow_to_CompCol(MKL_INT, MKL_INT, MKL_INT, double*, MKL_INT*, MKL_INT*,
		                   double **, MKL_INT **, MKL_INT **);
extern void    cCompRow_to_CompCol(MKL_INT, MKL_INT, MKL_INT, lscomplex*, MKL_INT*, MKL_INT*,
		                   lscomplex **, MKL_INT **, MKL_INT **);
extern void    zCompRow_to_CompCol(MKL_INT, MKL_INT, MKL_INT, ldcomplex*, MKL_INT*, MKL_INT*,
		                   ldcomplex **, MKL_INT **, MKL_INT **);
extern void    sfill (float *, MKL_INT, float);
extern void    dfill (double *, MKL_INT, double);
extern void    cfill (lscomplex *, MKL_INT, lscomplex);
extern void    zfill (ldcomplex *, MKL_INT, ldcomplex);
extern void    sinf_norm_error (MKL_INT, SuperMatrix *, float *);
extern void    dinf_norm_error (MKL_INT, SuperMatrix *, double *);
extern void    cinf_norm_error (MKL_INT, SuperMatrix *, lscomplex *);
extern void    zinf_norm_error (MKL_INT, SuperMatrix *, ldcomplex *);
//  extern void    PrintPerf (SuperMatrix *, SuperMatrix *, mem_usage_t *,
//                           float, float, float *, float *, char *);

/* Routines for debugging */
extern void    sPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    dPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    cPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    zPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    sPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    dPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    cPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    zPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    sPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    dPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    cPrint_Dense_Matrix(char *, SuperMatrix *);
extern void    zPrint_Dense_Matrix(char *, SuperMatrix *);
//     extern void    print_lu_col(char *, MKL_INT, MKL_INT, MKL_INT *, sGlobalLU_t *);
//     extern void    check_tempv(MKL_INT, float *);

/* Reordering routine */


#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_SP_DEFS */


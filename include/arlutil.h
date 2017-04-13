/*
   ARPACK++ v1.2 2/20/2000
   c++ MKL_INTerface to ARPACK code.

   MODULE ARLUtil.h.
   Unaltered copy of util.h (from SuperLU package)
   and superlu_enum_consts.h
*/



#ifndef __SUPERLU_ENUM_CONSTS /* allow multiple inclusions */
#define __SUPERLU_ENUM_CONSTS

/***********************************************************************
 * Enumerate types
 ***********************************************************************/
typedef enum {NO, YES}                                          yes_no_t;
typedef enum {DOFACT, SamePattern, SamePattern_SameRowPerm, FACTORED} fact_t;
typedef enum {NOROWPERM, LargeDiag, MY_PERMR}                   rowperm_t;
typedef enum {NATURAL, MMD_ATA, MMD_AT_PLUS_A, COLAMD,
	      METIS_AT_PLUS_A, PARMETIS, ZOLTAN, MY_PERMC}      colperm_t;
typedef enum {NOTRANS, TRANS, CONJ}                             trans_t;
typedef enum {NOEQUIL, ROW, COL, BOTH}                          DiagScale_t;
typedef enum {NOREFINE, SLU_SINGLE=1, SLU_DOUBLE, SLU_EXTRA}    IterRefine_t;
typedef enum {LUSUP, UCOL, LSUB, USUB, LLVL, ULVL}              MemType;
typedef enum {HEAD, TAIL}                                       stack_end_t;
typedef enum {SYSTEM, USER}                                     LU_space_t;
typedef enum {ONE_NORM, TWO_NORM, INF_NORM}			norm_t;
typedef enum {SILU, SMILU_1, SMILU_2, SMILU_3}			milu_t;
#if 0
typedef enum {NODROP		= 0x0000,
	      DROP_BASIC	= 0x0001, /* ILU(tau) */
	      DROP_PROWS	= 0x0002, /* ILUTP: keep p maximum rows */
	      DROP_COLUMN	= 0x0004, /* ILUTP: for j-th column, 
					     p = gamma * nnz(A(:,j)) */
	      DROP_AREA 	= 0x0008, /* ILUTP: for j-th column, use
					     nnz(F(:,1:j)) / nnz(A(:,1:j))
					     to limit memory growth  */
	      DROP_SECONDARY	= 0x000E, /* PROWS | COLUMN | AREA */
	      DROP_DYNAMIC	= 0x0010,
	      DROP_INTERP	= 0x0100}			rule_t;
#endif


/* 
 * The following enumerate type is used by the statistics variable 
 * to keep track of flop count and time spent at various stages.
 *
 * Note that not all of the fields are disjoint.
 */
typedef enum {
    COLPERM, /* find a column ordering that minimizes fills */
    ROWPERM, /* find a row ordering maximizes diagonal. */
    RELAX,   /* find artificial supernodes */
    ETREE,   /* compute column etree */
    EQUIL,   /* equilibrate the original matrix */
    SYMBFAC, /* symbolic factorization. */
    DIST,    /* distribute matrix. */
    FACT,    /* perform LU factorization */
    COMM,    /* communication for factorization */
    SOL_COMM,/* communication for solve */
    RCOND,   /* estimate reciprocal condition number */
    SOLVE,   /* forward and back solves */
    REFINE,  /* perform iterative refinement */
    TRSV,    /* fraction of FACT spent in xTRSV */
    GEMV,    /* fraction of FACT spent in xGEMV */
    FERR,    /* estimate error bounds after iterative refinement */
    NPHASES  /* total number of phases */
} PhaseType;


#endif /* __SUPERLU_ENUM_CONSTS */



#ifndef __SUPERLU_UTIL /* allow multiple inclusions */
#define __SUPERLU_UTIL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*
#ifndef __STDC__
#include <malloc.h>
#endif
*/
#include <assert.h>

/***********************************************************************
 * Macros
 ***********************************************************************/
#define FIRSTCOL_OF_SNODE(i)	(xsup[i])
/* No of marker arrays used in the symbolic factorization,
   each of size n */
#define NO_MARKER     3
#define NUM_TEMPV(m,w,t,b)  ( SUPERLU_MAX(m, (t + b)*w) )

#ifndef USER_ABORT
#define USER_ABORT(msg) superlu_abort_and_exit(msg)
#endif

#define ABORT(err_msg) \
 { char msg[256];\
   sprintf(msg,"%s at line %d in file %s\n",err_msg,__LINE__, __FILE__);\
   USER_ABORT(msg); }


#ifndef USER_MALLOC
#if 1
#define USER_MALLOC(size) superlu_malloc(size)
#else
/* The following may check out some uninitialized data */
#define USER_MALLOC(size) memset (superlu_malloc(size), '\x0F', size)
#endif
#endif

#define SUPERLU_MALLOC(size) USER_MALLOC(size)

#ifndef USER_FREE
#define USER_FREE(addr) superlu_free(addr)
#endif

#define SUPERLU_FREE(addr) USER_FREE(addr)

#define CHECK_MALLOC(where) {                 \
    extern MKL_INT superlu_malloc_total;        \
    printf("%s: malloc_total %d Bytes\n",     \
	   where, superlu_malloc_total); \
}

#define SUPERLU_MAX(x, y) 	( (x) > (y) ? (x) : (y) )
#define SUPERLU_MIN(x, y) 	( (x) < (y) ? (x) : (y) )

/*********************************************************
 * Macros used for easy access of sparse matrix entries. *
 *********************************************************/
#define L_SUB_START(col)     ( Lstore->rowind_colptr[col] )
#define L_SUB(ptr)           ( Lstore->rowind[ptr] )
#define L_NZ_START(col)      ( Lstore->nzval_colptr[col] )
#define L_FST_SUPC(superno)  ( Lstore->sup_to_col[superno] )
#define U_NZ_START(col)      ( Ustore->colptr[col] )
#define U_SUB(ptr)           ( Ustore->rowind[ptr] )


/***********************************************************************
 * Constants 
 ***********************************************************************/
#define EMPTY	(-1)
/*#define NO	(-1)*/
#define FALSE	0
#define TRUE	1

#define NO_MEMTYPE  4      /* 0: lusup;
			      1: ucol;
			      2: lsub;
			      3: usub */

#define GluIntArray(n)   (5 * (n) + 5)

/* Dropping rules */
#define  NODROP	        ( 0x0000 )
#define	 DROP_BASIC	( 0x0001 )  /* ILU(tau) */
#define  DROP_PROWS	( 0x0002 )  /* ILUTP: keep p maximum rows */
#define  DROP_COLUMN	( 0x0004 )  /* ILUTP: for j-th column, 
				              p = gamma * nnz(A(:,j)) */
#define  DROP_AREA 	( 0x0008 )  /* ILUTP: for j-th column, use
 		 			      nnz(F(:,1:j)) / nnz(A(:,1:j))
					      to limit memory growth  */
#define  DROP_SECONDARY	( 0x000E )  /* PROWS | COLUMN | AREA */
#define  DROP_DYNAMIC	( 0x0010 )  /* adaptive tau */
#define  DROP_INTERP	( 0x0100 )  /* use MKL_INTerpolation */


#if 1
#define MILU_ALPHA (1.0e-2) /* multiple of drop_sum to be added to diagonal */
#else
#define MILU_ALPHA  1.0 /* multiple of drop_sum to be added to diagonal */
#endif


/***********************************************************************
 * Type definitions
 ***********************************************************************/
typedef float    flops_t;
typedef unsigned char Logical;

/* 
 *-- This contains the options used to control the solution process.
 *
 * Fact   (fact_t)
 *        Specifies whether or not the factored form of the matrix
 *        A is supplied on entry, and if not, how the matrix A should
 *        be factorizaed.
 *        = DOFACT: The matrix A will be factorized from scratch, and the
 *             factors will be stored in L and U.
 *        = SamePattern: The matrix A will be factorized assuming
 *             that a factorization of a matrix with the same sparsity
 *             pattern was performed prior to this one. Therefore, this
 *             factorization will reuse column permutation vector 
 *             ScalePermstruct->perm_c and the column elimination tree
 *             LUstruct->etree.
 *        = SamePattern_SameRowPerm: The matrix A will be factorized
 *             assuming that a factorization of a matrix with the same
 *             sparsity	pattern and similar numerical values was performed
 *             prior to this one. Therefore, this factorization will reuse
 *             both row and column scaling factors R and C, both row and
 *             column permutation vectors perm_r and perm_c, and the
 *             data structure set up from the previous symbolic factorization.
 *        = FACTORED: On entry, L, U, perm_r and perm_c contain the 
 *              factored form of A. If DiagScale is not NOEQUIL, the matrix
 *              A has been equilibrated with scaling factors R and C.
 *
 * Equil  (yes_no_t)
 *        Specifies whether to equilibrate the system (scale A's row and
 *        columns to have unit norm).
 *
 * ColPerm (colperm_t)
 *        Specifies what type of column permutation to use to reduce fill.
 *        = NATURAL: use the natural ordering 
 *        = MMD_ATA: use minimum degree ordering on structure of A'*A
 *        = MMD_AT_PLUS_A: use minimum degree ordering on structure of A'+A
 *        = COLAMD: use approximate minimum degree column ordering
 *        = MY_PERMC: use the ordering specified by the user
 *         
 * Trans  (trans_t)
 *        Specifies the form of the system of equations:
 *        = NOTRANS: A * X = B        (No transpose)
 *        = TRANS:   A**T * X = B     (Transpose)
 *        = CONJ:    A**H * X = B     (Transpose)
 *
 * IterRefine (IterRefine_t)
 *        Specifies whether to perform iterative refinement.
 *        = NO: no iterative refinement
 *        = SLU_SINGLE: perform iterative refinement in single precision
 *        = SLU_DOUBLE: perform iterative refinement in double precision
 *        = SLU_EXTRA: perform iterative refinement in extra precision
 *
 * DiagPivotThresh (double, in [0.0, 1.0]) (only for sequential SuperLU)
 *        Specifies the threshold used for a diagonal entry to be an
 *        acceptable pivot.
 *
 * SymmetricMode (yest_no_t)
 *        Specifies whether to use symmetric mode. Symmetric mode gives 
 *        preference to diagonal pivots, and uses an (A'+A)-based column
 *        permutation algorithm.
 *
 * PivotGrowth (yes_no_t)
 *        Specifies whether to compute the reciprocal pivot growth.
 *
 * ConditionNumber (ues_no_t)
 *        Specifies whether to compute the reciprocal condition number.
 *
 * RowPerm (rowperm_t) (only for SuperLU_DIST or ILU)
 *        Specifies whether to permute rows of the original matrix.
 *        = NO: not to permute the rows
 *        = LargeDiag: make the diagonal large relative to the off-diagonal
 *        = MY_PERMR: use the permutation given by the user
 *
 * ILU_DropRule (MKL_INT)
 *        Specifies the dropping rule:
 *	  = DROP_BASIC:   Basic dropping rule, supernodal based ILUTP(tau).
 *	  = DROP_PROWS:   Supernodal based ILUTP(p,tau), p = gamma * nnz(A)/n.
 *	  = DROP_COLUMN:  Variant of ILUTP(p,tau), for j-th column,
 *			      p = gamma * nnz(A(:,j)).
 *	  = DROP_AREA:    Variation of ILUTP, for j-th column, use
 *			      nnz(F(:,1:j)) / nnz(A(:,1:j)) to control memory.
 *	  = DROP_DYNAMIC: Modify the threshold tau during factorizaion:
 *			  If nnz(L(:,1:j)) / nnz(A(:,1:j)) > gamma
 *				  tau_L(j) := MIN(tau_0, tau_L(j-1) * 2);
 *			  Otherwise
 *				  tau_L(j) := MAX(tau_0, tau_L(j-1) / 2);
 *			  tau_U(j) uses the similar rule.
 *			  NOTE: the thresholds used by L and U are separate.
 *	  = DROP_INTERP:  Compute the second dropping threshold by
 *	                  MKL_INTerpolation instead of sorting (default).
 *  		          In this case, the actual fill ratio is not
 *			  guaranteed to be smaller than gamma.
 *   	  Note: DROP_PROWS, DROP_COLUMN and DROP_AREA are mutually exclusive.
 *	  ( Default: DROP_BASIC | DROP_AREA )
 *
 * ILU_DropTol (double)
 *        numerical threshold for dropping.
 *
 * ILU_FillFactor (double) 
 *        Gamma in the secondary dropping.
 *
 * ILU_Norm (norm_t)
 *        Specify which norm to use to measure the row size in a
 *        supernode: infinity-norm, 1-norm, or 2-norm.
 *
 * ILU_FillTol (double)
 *        numerical threshold for zero pivot perturbation.
 *
 * ILU_MILU (milu_t)
 *        Specifies which version of MILU to use.
 *
 * ILU_MILU_Dim (double) 
 *        Dimension of the PDE if available.
 *
 * ReplaceTinyPivot (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether to replace the tiny diagonals by
 *        sqrt(epsilon)*||A|| during LU factorization.
 *
 * SolveInitialized (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether the initialization has been performed to the
 *        triangular solve.
 *
 * RefineInitialized (yes_no_t) (only for SuperLU_DIST)
 *        Specifies whether the initialization has been performed to the
 *        sparse matrix-vector multiplication routine needed in iterative
 *        refinement.
 *
 * PrintStat (yes_no_t)
 *        Specifies whether to print the solver's statistics.
 */
typedef struct {
    fact_t        Fact;
    yes_no_t      Equil;
    colperm_t     ColPerm;
    trans_t       Trans;
    IterRefine_t  IterRefine;
    double        DiagPivotThresh;
    yes_no_t      SymmetricMode;
    yes_no_t      PivotGrowth;
    yes_no_t      ConditionNumber;
    rowperm_t     RowPerm;
    MKL_INT 	  ILU_DropRule;
    double	  ILU_DropTol;    /* threshold for dropping */
    double	  ILU_FillFactor; /* gamma in the secondary dropping */
    norm_t	  ILU_Norm;       /* infinity-norm, 1-norm, or 2-norm */
    double	  ILU_FillTol;    /* threshold for zero pivot perturbation */
    milu_t	  ILU_MILU;
    double	  ILU_MILU_Dim;   /* Dimension of PDE (if available) */
    yes_no_t      ParSymbFact;
    yes_no_t      ReplaceTinyPivot; /* used in SuperLU_DIST */
    yes_no_t      SolveInitialized;
    yes_no_t      RefineInitialized;
    yes_no_t      PrintStat;
    MKL_INT           nnzL, nnzU;      /* used to store nnzs for now       */
    MKL_INT           num_lookaheads;  /* num of levels in look-ahead      */
    yes_no_t      lookahead_etree; /* use etree computed from the
				      serial symbolic factorization */
    yes_no_t      SymPattern;      /* symmetric factorization          */
} superlu_options_t;

/*! \brief Headers for 4 types of dynamatically managed memory */
typedef struct e_node {
    MKL_INT size;      /* length of the memory that has been used */
    void *mem;     /* pointer to the new malloc'd store */
} ExpHeader;

typedef struct {
    MKL_INT  size;
    MKL_INT  used;
    MKL_INT  top1;  /* grow upward, relative to &array[0] */
    MKL_INT  top2;  /* grow downward */
    void *array;
} LU_stack_t;

typedef struct {
    MKL_INT     *panel_histo; /* histogram of panel size distribution */
    double  *utime;       /* running time at various phases */
    flops_t *ops;         /* operation count at various phases */
    MKL_INT     TinyPivots;   /* number of tiny pivots */
    MKL_INT     RefineSteps;  /* number of iterative refinement steps */
    MKL_INT     expansions;   /* number of memory expansions */
} SuperLUStat_t;

typedef struct {
    float for_lu;
    float total_needed;
} mem_usage_t;


/***********************************************************************
 * Prototypes
 ***********************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

extern void    Destroy_SuperMatrix_Store(SuperMatrix *);
extern void    Destroy_CompCol_Matrix(SuperMatrix *);
extern void    Destroy_CompRow_Matrix(SuperMatrix *);
extern void    Destroy_SuperNode_Matrix(SuperMatrix *);
extern void    Destroy_CompCol_Permuted(SuperMatrix *);
extern void    Destroy_Dense_Matrix(SuperMatrix *);
extern void    get_perm_c(MKL_INT, SuperMatrix *, MKL_INT *);
extern void    set_default_options(superlu_options_t *options);
extern void    ilu_set_default_options(superlu_options_t *options);
extern void    sp_preorder (superlu_options_t *, SuperMatrix*, MKL_INT*, MKL_INT*,
			    SuperMatrix*);
extern void    superlu_abort_and_exit(char*);
extern void    *superlu_malloc (size_t);
extern MKL_INT     *intMalloc (MKL_INT);
extern MKL_INT     *intCalloc (MKL_INT);
extern void    superlu_free (void*);
extern void    SetIWork (MKL_INT, MKL_INT, MKL_INT, MKL_INT *, MKL_INT **, MKL_INT **, MKL_INT **,
                         MKL_INT **, MKL_INT **, MKL_INT **, MKL_INT **);
extern MKL_INT     sp_coletree (MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT, MKL_INT, MKL_INT *);
extern void    relax_snode (const MKL_INT, MKL_INT *, const MKL_INT, MKL_INT *, MKL_INT *);
extern void    heap_relax_snode (const MKL_INT, MKL_INT *, const MKL_INT, MKL_INT *, MKL_INT *);
extern MKL_INT     mark_relax(MKL_INT, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT *);
extern void    ilu_relax_snode (const MKL_INT, MKL_INT *, const MKL_INT, MKL_INT *,
				int *, MKL_INT *);
extern void    ilu_heap_relax_snode (const MKL_INT, MKL_INT *, const MKL_INT, MKL_INT *,
				     MKL_INT *, MKL_INT*);
extern void    resetrep_col (const MKL_INT, const MKL_INT *, MKL_INT *);
extern MKL_INT     spcoletree (MKL_INT *, MKL_INT *, MKL_INT *, MKL_INT, MKL_INT, MKL_INT *);
extern MKL_INT     *TreePostorder (MKL_INT, MKL_INT *);
extern double  SuperLU_timer_ ();
extern MKL_INT     sp_ienv (MKL_INT);
extern MKL_INT     lsame_ (char *, char *);
extern MKL_INT     xerbla_ (char *, MKL_INT *);
extern void    ifill (MKL_INT *, MKL_INT, MKL_INT);
extern void    snode_profile (MKL_INT, MKL_INT *);
extern void    super_stats (MKL_INT, MKL_INT *);
extern void    check_repfnz(MKL_INT, MKL_INT, MKL_INT, MKL_INT *);
extern void    PrintSumm (char *, MKL_INT, MKL_INT, MKL_INT);
extern void    StatInit(SuperLUStat_t *);
extern void    StatPrint (SuperLUStat_t *);
extern void    StatFree(SuperLUStat_t *);
extern void    print_panel_seg(MKL_INT, MKL_INT, MKL_INT, MKL_INT, MKL_INT *, MKL_INT *);
extern MKL_INT     print_int_vec(char *,MKL_INT, MKL_INT *);
extern MKL_INT     slu_PrintInt10(char *, MKL_INT, MKL_INT *);

#ifdef __cplusplus
  }
#endif

#endif /* __SUPERLU_UTIL */

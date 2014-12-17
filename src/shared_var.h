#ifndef Shared_Var_h
#define Shared_Var_h

#define MPI_COMM_WORLD ((MPI_Comm)0x44000000)

typedef int MPI_Comm;

typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)0x4c000101)
#define MPI_SIGNED_CHAR    ((MPI_Datatype)0x4c000118)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)0x4c000102)
#define MPI_BYTE           ((MPI_Datatype)0x4c00010d)
#define MPI_WCHAR          ((MPI_Datatype)0x4c00040e)
#define MPI_SHORT          ((MPI_Datatype)0x4c000203)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)0x4c000204)
#define MPI_INT            ((MPI_Datatype)0x4c000405)
#define MPI_UNSIGNED       ((MPI_Datatype)0x4c000406)
#define MPI_LONG           ((MPI_Datatype)0x4c000807)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)0x4c000808)
#define MPI_FLOAT          ((MPI_Datatype)0x4c00040a)
#define MPI_DOUBLE         ((MPI_Datatype)0x4c00080b)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)0x4c00100c)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)0x4c000809)
#define MPI_UNSIGNED_LONG_LONG ((MPI_Datatype)0x4c000819)
#define MPI_LONG_LONG      MPI_LONG_LONG_INT


extern "C" 
{
    int MPI_Init(int *, char ***);
    int MPI_Finalize(void);
    int MPI_Dims_create(int, int, int *);
    int MPI_Barrier( MPI_Comm comm );
    void blacs_pinfo_ ( int *mypnum, int *nprocs );
    void blacs_get_ ( int *ConTxt, int *what, int *val );
    void blacs_gridinit_ ( int *ConTxt, char *order, int *nprow, int *npcol );
    void blacs_gridexit_ ( int *ConTxt );
    void blacs_pcoord_ ( int *ConTxt, int *nodenum, int *prow, int *pcol );
    void descinit_ ( int*, int*, int*, int*, int*, int*, int*, int*, int*, int* );
    void pdpotrf_ ( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
    void pdpotri_ ( char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info );
    void pdsymm_( char *side, char *uplo, int *m, int *n, double *alpha, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb,
                  int *descb, double *beta, double *c, int *ic, int *jc, int *descc );
    void pddot_( int *n, double *dot, double *x, int *ix, int *jx, int *descx, int *incx, double *y, int *iy, int *jy, int *descy, int *incy );
    void dgesd2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rdest, int *cdest );
    void dgerv2d_ ( int *ConTxt, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc );
    void dgemm_ ( const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *lda, const double *b, const int *ldb, const double *beta, double *c, const int *ldc );
}

typedef struct MPI_Status {
    int count;
    int cancelled;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
    
} MPI_Status;


class CSRdouble;
class ParDiSO;

extern "C"{
    int MPI_Send(void*, int, MPI_Datatype, int, int, MPI_Comm);
    int MPI_Recv(void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *);
    int MPI_Get_count(MPI_Status*, MPI_Datatype, int*);
    void blacs_barrier_ ( int*, char* );
    int blacs_pnum_ ( int *ConTxt, int *prow, int *pcol );
    double MPI_Wtime();
}

void printdense ( int m, int n, double *mat, const char *filename );
int read_in_BD ( int * DESCD, double * Dmat, CSRdouble& BT_i, CSRdouble& B_j, CSRdouble& Btsparse ) ;
int read_input ( char* filename ) ;
int make_Sij_sparse_parallel (CSRdouble& A, CSRdouble& BT_i, CSRdouble& B_j, double* T_ij, int lld_Tij );
int make_Sij_parallel_denseB(CSRdouble& A, CSRdouble& BT_i, CSRdouble& B_j, double * T_ij, int lld_T, double* AB_sol) ;
void dense2CSR ( double *mat, int m, int n, CSRdouble& A );
void dense2CSR_sub ( double *mat, int m, int n, int lld_mat, CSRdouble& A, int startrow, int startcol ) ;
void CSR2dense ( CSRdouble& matrix, double* T_ij ) ;
void CSR2dense_lld ( CSRdouble& matrix,double *dense, int lld_dense ) ;

void mult_colsA_colsC ( CSRdouble& A, double *B, int lld_B, int Acolstart, int Acolstop, int Ccolstart, int Ccolstop,
                        CSRdouble& C, bool trans );


#endif


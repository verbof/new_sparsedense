#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <algorithm>

#include "shared_var.h"
#include "config.hpp"
#include "CSRdouble.hpp"
#include "IO.hpp"
#include "ParDiSO.hpp"
#include "Solver.hpp"
#include "SelInvProcess.hpp"
#include "smat.h"
#include "timing.hpp"
#include "Options.hpp"
#include "Keywords.hpp"

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
}


void computeDiagonalUsingMPI(int* argc, char*** argv)
{

    Options& options = Options::reference();

    SelInvProcess process(*argc, argv);


    process.readMatrices(options.blocksize, options.dimD, options.dimA, (options.fileA).c_str(), (options.fileB).c_str(), (options.fileD).c_str());

    if (options.printDebug)
    {
        process.printDebugInfo(options.blocksize);
    }

    process.makeSchur();
    process.factorizeSchur();

    process.extractDiagInvD();
    
    if (process.iam == 0)
    {
        process.invertA();
    }


    process.callBlacsBarrier();

    process.computeDiagInvA();

    process.saveDiagonal(options.output_prefix);


}

void computeDiagonalUsingParDiSO()
{
    // Nothing, at the moment...
}


int main(int argc, char **argv)
{
    if (argc != 3)
    {
        cout << "usage: " << argv[0] << "mode MPI_grid_size inputfilename" << endl;
        cout << "examp: " << argv[0] << "MPI  5 simplecase" << endl;
        cout << "examp: " << argv[0] << "PARDISO 4 simplecase" << endl;

        exit(1);
    }

    string mode = argv[1];

    // convert mode to uppercase
    std::transform(mode.begin(), mode.end(), mode.begin(), (int(*)(int))toupper);

    int N                 = (int) strtol(argv[2], NULL, 0);
    const char* inputfile = argv[3];

    Options& options = Options::reference();
    Keywords keywords(inputfile, &options);
    keywords.read();

    options.print();

    if (mode == "MPI")
    {
        computeDiagonalUsingMPI(&argc, &argv);
    }
    else if (mode == "PARDISO")
    {
        computeDiagonalUsingParDiSO();
    }
    else
    {
        cout << "your choice: \"" << mode << "\" is not supported, ... yet" << endl;
    }

    return 0;
}


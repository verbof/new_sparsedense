#include "SelInvProcess.hpp"


SelInvProcess::~SelInvProcess()
{

}


SelInvProcess::SelInvProcess(int argc, char*** argv)
{

    //Initialise MPI and some MPI-variables
    info = MPI_Init(&argc, argv);
    if (info != 0)
    {
        printf("Error in MPI initialisation: %d\n", info);
        return;
    }

    int* position = new int[2];
    int* dims     = new int[2];

    //BLACS is the interface used by PBLAS and ScaLAPACK on top of MPI

    blacs_pinfo_(&iam, &size);          //determine the number of processes involved

    if ((iam == 0) && (floor(sqrt(size))*round(sqrt(size))) != size)
    {
        cout << "*** Number of procesors ( " << size << " ) is not an exact square! ***" << endl;
        exit(1);
    }

    info = MPI_Dims_create(size, 2, dims);  //determine the best 2D cartesian grid with the number of processes
    if (info != 0)
    {
        printf ( "Error in MPI creation of dimensions: %d\n", info );
        return;
    }

    I = dims[0];
    J = dims[1];

    int i_negone = -1, i_zero = 0;

    blacs_get_(&i_negone, &i_zero, &ICTXT2D);

    blacs_gridinit_(&ICTXT2D, "R", &M, &N);     //Initialisation of the BLACS process grid, which is referenced as ICTXT2D

    if (M != N)
    {
        cout << "*** Processors grid is not square. Dimensions: M=" << M << ", N=" << N << "\n" << endl;
        exit(1);

    }

    //The rank (iam) of the process is mapped to a 2D grid: position= (process row, process column)
    blacs_pcoord_(&ICTXT2D, &iam, &I, &J);

    if (I == -1)
    {
        printf("Error in processor grid\n");
        return;
    }

    callBlacsBarrier();

    if (I == 0 && J == 0)
    {
        cout << "Processor grid made by " << size << " (" << M << "x" << N << ") processors! \n" << endl;
    }


    solver = new Solver(iam, ICTXT2D, options->isBsparse);

}

void SelInvProcess::callBlacsBarrier()
{
    char all = 'A';
    blacs_barrier_(&ICTXT2D, &all);
}


void SelInvProcess::readMatrices()
{

    int blocksize = options->blocksize;
    int dimD      = options->dimD;
    int dimA      = options->dimA;


    //Define number of blocks needed to store a complete column/row of D
    //Dblocks = dimD % blocksize == 0 ? dimD/blocksize : dimD/blocksize + 1;

    Dblocks = ceil(dimD/(double)blocksize);


    //Define the number of rowblocks (and columnblocks) needed by the current process to store its part of the dense matrix D
    //Drows = (Dblocks-I) % M == 0 ? (Dblocks-I)/M : (Dblocks-I)/M + 1;
    //Dcols = (Dblocks-J) % N == 0 ? (Dblocks-J)/N : (Dblocks-J)/N + 1;

    Drows = ceil((Dblocks-I)/(double)M);
    Dcols = Dcols < 1 ? 1 : Dcols;

    Dcols = ceil((Dblocks-J)/(double)N);
    Drows = Drows < 1 ? 1 : Drows;

    //Define the local leading dimension of D (keeping in mind that matrices are always stored column-wise)
    lld_D = Drows * blocksize;

    //Initialise the descriptor of the dense distributed matrix
    DESCD = (int*) malloc(DLEN_*sizeof(int));
    if (DESCD == NULL)
    {
        printf("Unable to allocate memory for descriptor for C. \n");
        return;
    }

    // D with dimensions (dimD,dimD) is distributed over all processes in ICTXT2D, with the first element in process (0,0)
    // D is distributed into blocks of size (blocksize,blocksize), having a local leading dimension lld_D in this specific process

    int i_zero = 0;

    descinit_(DESCD, &dimD, &dimD, &blocksize, &blocksize, &i_zero, &i_zero, &ICTXT2D, &lld_D, &info);
    if (info != 0)
    {
        printf("Descriptor of matrix C returns info: %d. \n", info);
        return;
    }

    // Allocate the space necessary to store the part of D that is held into memory of this process.
    D_ij.allocate(Drows*blocksize, Dcols*blocksize);

    read_in_BD(D_ij.data);

    callBlacsBarrier();

    if (iam == 0)
        printf("Matrices B & D read in. \n");


    //Now every process has to read in the sparse matrix A


    A.loadFromFileSym(options->fileA.c_str());
    A.matrixType = SYMMETRIC;

    assert(A.nrows == dimA);
    assert(A.ncols == dimA);

    callBlacsBarrier();

}


void SelInvProcess::read_in_BD(double* Dmat)
{

    FILE *fD;
    int ni, i, j;
    int nstrips;

    int blocksize = options->blocksize;
    int dimD      = options->dimD;
    int dimA      = options->dimA;

    MPI_Status status;

    // pcol= * ( position+1 );

    fD = fopen(options->fileD.c_str(), "rb");

    if (fD == NULL)
    {
        printf("Error opening file! \n");
        return; 
    }

    nstrips = dimD % (blocksize*N) == 0 ?  dimD/(blocksize*N) : dimD/(blocksize*N) + 1;
           // dimD % ( blocksize ** ( dims + 1 ) ) == 0 ?  dimD / ( blocksize ** ( dims + 1 ) ) : ( dimD / ( blocksize ** ( dims + 1 ) ) ) + 1;

    // Set up of matrix D and B per strip of T'

    for (ni = 0; ni < nstrips; ni++)
    {
        if ((Dblocks-1) % M == I  &&  dimD % blocksize != 0)
        {
            if (ni == 0)
            {
                info = fseek(fD, (long) (J * blocksize * dimD * sizeof(double)), SEEK_SET);
                    // fseek ( fD, ( long ) ( pcol * blocksize * ( dimD ) * sizeof ( double ) ),SEEK_SET );
                if (info != 0)
                {
                    printf("Error in setting correct begin position for reading D file \n \t Processor (%d,%d), error: %d. \n", I, J, info);
                    return;
                }
            }
            else
            {
                info = fseek(fD, (long) ((N-1) * blocksize * dimD * sizeof(double)), SEEK_CUR);
                    // fseek ( fD, ( long ) ( blocksize * ( * ( dims + 1 ) - 1 ) * ( dimD ) * sizeof ( double ) ), SEEK_CUR );
                if (info != 0)
                {
                    printf("Error in setting correct begin position for reading D file \n \t Processor (%d,%d), error: %d. \n", I, J, info);
                    return;
                }
            }
            for (i = 0; i < blocksize; i++)
            {
                info = fseek(fD, (long) (I * blocksize * sizeof(double)), SEEK_CUR);
                    // fseek ( fD, ( long ) ( blocksize **position * sizeof ( double ) ), SEEK_CUR );
                if ( info != 0 )
                {
                    printf("Error in setting correct begin position for reading D file \n \t Processor (%d,%d), error: %d. \n", I, J, info);
                    return;
                }
                for (j = 0; j < Drows-1; j++)
                {
                    fread(Dmat + i*Drows*blocksize + j*blocksize + ni*blocksize*Drows*blocksize, sizeof(double), blocksize, fD);
                    info = fseek(fD, (long) ((M-1) * blocksize * sizeof(double)), SEEK_CUR);
                        
                    //fread ( Dmat + i * Drows * blocksize + j * blocksize + ni * blocksize * Drows * blocksize, sizeof ( double ), blocksize, fD );
                        // fseek ( fD, ( long ) ( ( ( *dims ) - 1 ) * blocksize * sizeof ( double ) ), SEEK_CUR );
                    
                    if ( info != 0 )
                    {
                        printf("Error in setting correct begin position for reading D file \n \t Processor (%d,%d), error: %d. \n", I, J, info);
                        return;
                    }
                }
                fread(Dmat + i*Drows*blocksize + j*blocksize + ni*blocksize*Drows*blocksize, sizeof(double), dimD % blocksize, fD);
                // fread ( Dmat + i * Drows * blocksize + j * blocksize + ni * blocksize * Drows * blocksize, sizeof ( double ), dimD % blocksize, fD );
            }
            //Normal read-in of the strips of T from a binary file (each time blocksize elements are read in)
        }
        else
        {
            if (ni == 0)
            {
                info = fseek(fD, (long) (J * blocksize * dimD * sizeof(double)), SEEK_SET);
                    // fseek ( fD, ( long ) ( pcol * blocksize * ( dimD ) * sizeof ( double ) ), SEEK_SET );

                if (info != 0)
                {
                    printf("Error in setting correct begin position for reading D file \n \t Processor (%d,%d), error: %d. \n", I, J, info);
                    return;
                }
            }
            else
            {
                info = fseek(fD, (long) (blocksize * (N-1) * dimD * sizeof(double)), SEEK_CUR);
                    // fseek ( fD, ( long ) ( blocksize * ( * ( dims + 1 ) - 1 ) * ( dimD ) * sizeof ( double ) ), SEEK_CUR );

                if (info != 0)
                {
                    printf("Error in setting correct begin position for reading D file \n \t Processor (%d,%d), error: %d. \n", I, J, info);
                    return;
                }
            }
            for (i = 0; i < blocksize; i++)
            {
                info = fseek (fD, (long) (I * blocksize * sizeof(double)), SEEK_CUR);
                    // fseek ( fD, ( long ) ( blocksize **position * sizeof ( double ) ), SEEK_CUR );

                if (info != 0)
                {
                    printf("Error in setting correct begin position for reading D file \n \t Processor (%d,%d), error: %d. \n", I, J, info);
                    return;
                }

                for (j = 0; j < Drows-1; j++)
                {
                    fread(Dmat + i*Drows*blocksize + j*blocksize + ni*blocksize*Drows*blocksize, sizeof(double), blocksize, fD);
                    info = fseek(fD, (long) ((M-1) * blocksize * sizeof(double)), SEEK_CUR);

                    // fread ( Dmat + i * Drows * blocksize + j * blocksize + ni * blocksize * Drows * blocksize, sizeof ( double ), blocksize, fD );
                        // fseek ( fD, ( long ) ( ( * ( dims ) - 1 ) * blocksize * sizeof ( double ) ), SEEK_CUR );

                    if (info != 0)
                    {
                    printf("Error in setting correct begin position for reading D file \n \t Processor (%d,%d), error: %d. \n", I, J, info);
                    return;
                    }
                }

                fread(Dmat + i*Drows*blocksize + j*blocksize + ni*blocksize*Drows*blocksize, sizeof(double), blocksize, fD);
                info = fseek(fD, (long) ((dimD - blocksize * ((Drows-1)*M + I + 1)) * sizeof(double)), SEEK_CUR);

                // fread ( Dmat + i * Drows * blocksize + j * blocksize + ni * blocksize * Drows * blocksize, sizeof ( double ), blocksize, fD );
                    // fseek ( fD, ( long ) ( ( dimD - blocksize * ( ( Drows - 1 ) **dims + *position + 1 ) ) * sizeof ( double ) ), SEEK_CUR );
                if (info != 0)
                {
                    printf("Error in setting correct begin position for reading D file \n \t Processor (%d,%d), error: %d. \n", I, J, info);
                    return;
                }
            }
        }

        callBlacsBarrier();
    }

    fclose(fD);

    // End of read-in

    // Matrix B
    if (iam != 0)
    {

        // Each process receives the necessary BT_i and B_j
        // Blocking sends are used, which is why the order of the receives is critical depending on the coordinates of the process

        int nonzeroes;
        int count;

        if (I >= J)
        {

            MPI_Recv(&nonzeroes, 1, MPI_INT, 0, iam, MPI_COMM_WORLD, &status);
            
            BT_i.allocate(blocksize * Drows, dimA, nonzeroes);
            
            MPI_Recv(&(BT_i.pRows[0]), blocksize*Drows + 1, MPI_INT, 0, iam + size, MPI_COMM_WORLD, &status);


            MPI_Get_count(&status, MPI_INT, &count);

            BT_i.nrows = count - 1;
            
            MPI_Recv(&(BT_i.pCols[0]), nonzeroes, MPI_INT,    0, iam + 2*size, MPI_COMM_WORLD, &status);
            MPI_Recv(&(BT_i.pData[0]), nonzeroes, MPI_DOUBLE, 0, iam + 3*size, MPI_COMM_WORLD, &status);
            MPI_Recv(&nonzeroes,       1,         MPI_INT,    0, iam + 4*size, MPI_COMM_WORLD, &status);

            B_j.allocate(blocksize*Dcols, dimA, nonzeroes);

            MPI_Recv(&(B_j.pRows[0]), blocksize*Dcols + 1, MPI_INT, 0, iam + 5*size, MPI_COMM_WORLD, &status);

            MPI_Get_count(&status, MPI_INT, &count);
            
            B_j.nrows = count - 1;
            
            MPI_Recv(&(B_j.pCols[0]), nonzeroes, MPI_INT,    0, iam + 6*size, MPI_COMM_WORLD, &status);
            MPI_Recv(&(B_j.pData[0]), nonzeroes, MPI_DOUBLE, 0, iam + 7*size, MPI_COMM_WORLD, &status);

            //Actually BT_j is sent, so it still needs to be transposed

            B_j.transposeIt(1);
        }
        else
        {
            MPI_Recv(&nonzeroes, 1, MPI_INT, 0, iam + 4*size, MPI_COMM_WORLD, &status);

            B_j.allocate(blocksize*Dcols, dimA, nonzeroes);

            MPI_Recv(&(B_j.pRows[0]), blocksize*Dcols + 1, MPI_INT, 0, iam + 5*size, MPI_COMM_WORLD, &status);

            MPI_Get_count(&status, MPI_INT, &count);

            B_j.nrows = count - 1;

            MPI_Recv(&(B_j.pCols[0]), nonzeroes, MPI_INT,    0, iam + 6*size, MPI_COMM_WORLD, &status);
            MPI_Recv(&(B_j.pData[0]), nonzeroes, MPI_DOUBLE, 0, iam + 7*size, MPI_COMM_WORLD, &status);

            B_j.transposeIt (1);

            MPI_Recv (&nonzeroes, 1, MPI_INT, 0, iam, MPI_COMM_WORLD, &status);

            BT_i.allocate(blocksize*Drows, dimA, nonzeroes);

            MPI_Recv(&(BT_i.pRows[0]), blocksize*Drows + 1, MPI_INT, 0, iam + size, MPI_COMM_WORLD, &status);

            MPI_Get_count(&status, MPI_INT, &count);

            BT_i.nrows = count - 1;

            MPI_Recv(&(BT_i.pCols[0]), nonzeroes, MPI_INT,    0, iam + 2*size, MPI_COMM_WORLD, &status);
            MPI_Recv(&(BT_i.pData[0]), nonzeroes, MPI_DOUBLE, 0, iam + 3*size, MPI_COMM_WORLD, &status);
        }
    }
    else
    {

        Btsparse.loadFromFile(options->fileB.c_str());

        assert(Btsparse.nrows == dimA);
        assert(Btsparse.ncols == dimD);

        Btsparse.transposeIt(1);

        // For each process row i BT_i is created which is also sent to processes in column i to become B_j.
        for (int rowproc = M-1; rowproc >= 0; rowproc--)
        {
            BT_i.ncols    = Btsparse.ncols;
            BT_i.nrows    = 0;
            BT_i.nonzeros = 0;

            int Drows_rowproc;

            if (rowproc != 0)
            {
                Drows_rowproc = (Dblocks-rowproc) % M == 0 ? (Dblocks-rowproc)/M : (Dblocks-rowproc)/M + 1;
                Drows_rowproc = Drows_rowproc < 1 ? 1 : Drows_rowproc;
            }
            else
            {
                Drows_rowproc = Drows;
            }

            for (i = 0; i < Drows_rowproc; i++)
            {
                //Each process in row i can hold several blocks of contiguous rows of D for which we need the corresponding rows of B_T
                // Therefore we use the function extendrows to create BT_i (see src/tools.cpp)

                BT_i.extendrows(Btsparse, (i*M + rowproc) * blocksize, blocksize);

            }

            for (int colproc = (rowproc == 0); colproc < N; colproc++)      // int colproc = (rowproc == 0 ? 1 : 0);
            {
                // int* curpos;
                int  rankproc;

                rankproc = blacs_pnum_(&ICTXT2D, &rowproc, &colproc);

                MPI_Send(&(BT_i.nonzeros), 1,              MPI_INT,    rankproc, rankproc,          MPI_COMM_WORLD);
                MPI_Send(&(BT_i.pRows[0]), BT_i.nrows + 1, MPI_INT,    rankproc, rankproc + size,   MPI_COMM_WORLD);
                MPI_Send(&(BT_i.pCols[0]), BT_i.nonzeros,  MPI_INT,    rankproc, rankproc + 2*size, MPI_COMM_WORLD);
                MPI_Send(&(BT_i.pData[0]), BT_i.nonzeros,  MPI_DOUBLE, rankproc, rankproc + 3*size, MPI_COMM_WORLD);

                //printf("BT_i's sent to processor %d\n",rankproc);

                rankproc = blacs_pnum_(&ICTXT2D, &colproc, &rowproc);

                MPI_Send(&(BT_i.nonzeros), 1,              MPI_INT,    rankproc, rankproc + 4*size, MPI_COMM_WORLD);
                MPI_Send(&(BT_i.pRows[0]), BT_i.nrows + 1, MPI_INT,    rankproc, rankproc + 5*size, MPI_COMM_WORLD);
                MPI_Send(&(BT_i.pCols[0]), BT_i.nonzeros,  MPI_INT,    rankproc, rankproc + 6*size, MPI_COMM_WORLD);
                MPI_Send(&(BT_i.pData[0]), BT_i.nonzeros,  MPI_DOUBLE, rankproc, rankproc + 7*size, MPI_COMM_WORLD);

                //printf("B_j's sent to processor %d\n",rankproc);
            }
        }

        B_j.make(BT_i.nrows, BT_i.ncols, BT_i.nonzeros, BT_i.pRows, BT_i.pCols, BT_i.pData);

        B_j.transposeIt(1);

    }
}






void SelInvProcess::makeSchur()
{
    DESCAB_sol = (int*) malloc(DLEN_*sizeof(int));     // AB_sol will contain the solution of A*X=B, distributed across the process rows. 
                                                            // Processes in the same process row possess the same part of AB_sol
    int blocksize = options->blocksize;
    int dimD      = options->dimD;
    int dimA      = options->dimA;
    int i_zero    = 0;


    if (DESCAB_sol == NULL)
    {
        printf("Unable to allocate memory for descriptor for AB_sol\n");
        return;
    }

    descinit_(DESCAB_sol, &dimA, &dimD, &dimA, &blocksize, &i_zero, &i_zero, &ICTXT2D, &dimA, &info ); 
                                                            // AB_sol (dimA, dimD) is distributed across all processes in ICTXT2D starting 
                                                            // from process (0,0) into blocks of size (dimA, blocksize)
    
    if (info != 0)
    {
        printf("Descriptor of matrix C returns info: %d\n", info);
        return;
    }

    // AB_sol = (double*) calloc(dimA*Dcols*blocksize, sizeof(double));
    AB_sol.allocate(dimA, Dcols*blocksize);

    callBlacsBarrier();

    // Each process calculates the Schur complement of the part of D at its disposal. (see src/schur.cpp)
    // The solution of A * X = B_j is stored in AB_sol (= A^{-1} * B_j)

    solver->make_Sij_parallel_denseB(A, BT_i, B_j, D_ij.data, lld_D, AB_sol.data);


    callBlacsBarrier();

}


void SelInvProcess::factorizeSchur()
{
    int dimD      = options->dimD;
    int i_one     = 1;

 

    pdpotrf_("U", &dimD, D_ij.data, &i_one, &i_one, DESCD, &info); //The Schur complement is factorised (by ScaLAPACK)

    if (info != 0)
    {
        printf("Cholesky decomposition of D was unsuccessful, error returned: %d\n", info);
        return;
    }

    callBlacsBarrier();

}

void SelInvProcess::invertSchur()
{

    int dimD      = options->dimD;
    int i_one     = 1;

    pdpotri_("U", &dimD, D_ij.data, &i_one, &i_one, DESCD, &info); //The Schur complement is inverteded (by ScaLAPACK)

    if (info != 0)
    {
        printf("Inverse of D was unsuccessful, error returned: %d\n", info);
        return;
    }

    callBlacsBarrier();
}


void SelInvProcess::extractDiagInvD()
{
    int blocksize = options->blocksize;
    int dimD      = options->dimD;
    int dimA      = options->dimA;
    int i_one     = 1;
    int i_zero    = 0;

    
    //InvD_T_Block = (double*) calloc(dimA + Dblocks*blocksize, sizeof(double));
    InvD_T_Block = new double[dimA + Dblocks*blocksize];

    //Diagonal elements of the (1,1) block of C^-1 are still distributed and here they are gathered in InvD_T_Block in the root process.
    if (I == J)    // if (*position == pcol)
    {
        for (int i = 0; i < dimD; i++)
        {
            if (J == (i / blocksize) % M)
            {
                int Dpos                    = i % blocksize + ((i / blocksize) / M) * blocksize;
                InvD_T_Block[dimA + i] = D_ij.data[Dpos + lld_D*Dpos];                             //*(InvD_T_Block + i + dimA) = *(D_ij + Dpos + lld_D*Dpos);
            }
        }

        for (int i = 0, j = 0; i < Dblocks; i++, j++)
        {
            if (j == M)
            {
                j = 0;
            }

            if (j == I)
            {
                dgesd2d_(&ICTXT2D, &blocksize, &i_one, &InvD_T_Block[dimA + i*blocksize], &blocksize, &i_zero, &i_zero);
            }

            if (I == 0)
            {
                dgerv2d_(&ICTXT2D, &blocksize, &i_one, &InvD_T_Block[dimA + i*blocksize], &blocksize, &j,      &j);
            }
        }
    }

}

void SelInvProcess::invertA()
{
    /*int pardiso_message_level = 1; int pardiso_mtype=-2; ParDiSO pardiso ( pardiso_mtype, pardiso_message_level );*/

    int nops  = 1;                          // int number_of_processors = 1;
    char* var = getenv("OMP_NUM_THREADS");
    
    if (var != NULL)
        sscanf(var, "%d", &nops);
    else
    {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }

    pardiso_var.iparm[2]   = 2;
    pardiso_var.iparm[3]   = nops;
    pardiso_var.iparm[8]   = 0;
    pardiso_var.iparm[11]  = 1;
    pardiso_var.iparm[13]  = 0;
    pardiso_var.iparm[28]  = 0;

    //This function calculates the factorisation of A once again so this might be optimized.

    pardiso_var.findInverseOfA(A);

    printf("Processor %d inverted matrix A! \n", iam);

    callBlacsBarrier();

}


void SelInvProcess::computeDiagInvA()
{

    int blocksize = options->blocksize;
    int dimD      = options->dimD;
    int dimA      = options->dimA;
    int i_one     = 1;
    int i_zero    = 0;
    double d_one  = 1.0;
    double d_zero = 0.0;


    // To minimize memory usage, and because only the diagonal elements of the inverse are needed, X' * S is calculated row by row
    // the diagonal element is calculated as the dot product of this row and the corresponding column of X. (X is solution of AX=B)
    //XSrow     = (double*) calloc(Dcols*blocksize, sizeof(double));
    XSrow.allocate(Dcols, blocksize);
    DESCXSROW = (int*)    malloc(DLEN_*sizeof(int));

    if (DESCXSROW == NULL)
    {
        printf("Unable to allocate memory for descriptor for AB_sol\n");
        return;
    }

    // XSrow (1,dimD) is distributed acrros processes of ICTXT2D starting from process (0,0) into blocks of size (1,blocksize)
    descinit_(DESCXSROW, &i_one, &dimD, &i_one, &blocksize, &i_zero, &i_zero, &ICTXT2D, &i_one, &info);
    if (info != 0)
    {
        printf("Descriptor of matrix C returns info: %d\n", info);
        return;
    }


    //Calculating diagonal elements 1 by 1 of the (0,0)-block of C^-1.

    for (int i = 1; i <= dimA; i++)
    {
        pdsymm_("R", "U", &i_one, &dimD, &d_one, D_ij.data, &i_one, &i_one, DESCD, AB_sol.data, &i, &i_one, DESCAB_sol, &d_zero, XSrow.data, &i_one, &i_one, DESCXSROW);
        pddot_(&dimD, &InvD_T_Block[i-1], AB_sol.data, &i, &i_one, DESCAB_sol, &dimA, XSrow.data, &i_one, &i_one, DESCXSROW, &i_one);
    }

    if (iam == 0)
    {
        for (int i = 0; i < dimA; i++)
        {
                       int j = A.pRows[i];
            InvD_T_Block[i] += A.pData[j];
        }
        
    }
    

}


void SelInvProcess::saveDiagonal()
{
    int dimD        = options->dimD;
    int dimA        = options->dimA;
    string filename = options->output_prefix;

    if (filename.at(filename.size()-1) != '/')
        filename += "/";

    filename += "diag_inverse_C_parallel.txt";
    
    printdense(dimA + dimD, 1, InvD_T_Block, filename.c_str());
}


void SelInvProcess::finalizeMPIandBLACS()
{
    callBlacsBarrier();

    blacs_gridexit_(&ICTXT2D);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

}


void SelInvProcess::printDebugInfo()
{

    callBlacsBarrier();

    CSRdouble Dmat, Dblock, Csparse;

    MPI_Status status;

    int dimD      = options  -> dimD;
    int blocksize = options  -> blocksize;
    int Drows     = D_ij.nrows / blocksize;
    int Dcols     = D_ij.ncols / blocksize;

    Dblock.nrows  = Dblocks   * blocksize;
    Dblock.ncols  = Dblocks   * blocksize;

    Dblock.allocate(Dblocks*blocksize, Dblocks*blocksize, 0);
    Dmat.allocate(0, 0, 0);

    for (int i = 0; i < Drows; ++i)
    {
        for (int j = 0; j < Dcols; ++j)
        {
            dense2CSR_sub(&(D_ij.data[i*blocksize + j*lld_D*blocksize]), blocksize, blocksize, lld_D, Dblock, (M*i + I)*blocksize, (N*j + J)*blocksize);
            //dense2CSR_sub(D + i * blocksize + j * lld_D * blocksize, blocksize, blocksize, lld_D, Dblock, ( * ( dims) * i + *position ) *blocksize, ( * ( dims + 1 ) * j + pcol ) *blocksize);

            if (Dblock.nonzeros > 0)
            {
                if (Dmat.nonzeros == 0)
                {
                    Dmat.make2(Dblock.nrows, Dblock.ncols, Dblock.nonzeros, Dblock.pRows, Dblock.pCols, Dblock.pData);
                }
                else
                {
                    Dmat.addBCSR(Dblock);
                }
            }
            Dblock.clear();
        }
    }

    callBlacsBarrier();

    if (iam != 0)
    {
        //Each process other than root sends its Dmat to the root process.
        MPI_Send (&(Dmat.nonzeros), 1,              MPI_INT,    0, iam,            MPI_COMM_WORLD);
        MPI_Send (&(Dmat.pRows[0]), Dmat.nrows + 1, MPI_INT,    0, iam + size,     MPI_COMM_WORLD);
        MPI_Send (&(Dmat.pCols[0]), Dmat.nonzeros,  MPI_INT,    0, iam + 2 * size, MPI_COMM_WORLD);
        MPI_Send (&(Dmat.pData[0]), Dmat.nonzeros,  MPI_DOUBLE, 0, iam + 3 * size, MPI_COMM_WORLD);
        Dmat.clear();
        Btsparse.clear();
    }
    else
    {
        for (int i = 1; i < size; i++)
        {
            // The root process receives parts of Dmat sequentially from all processes and directly adds them together.
            int nonzeroes;

            MPI_Recv(&nonzeroes, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            // MPI_Get_count(&status, MPI_INT, &count); printf("Process 0 received %d elements of process %d\n",count,i);

            if (nonzeroes > 0)
            {
                //printf("Nonzeroes : %d\n ",nonzeroes);

                Dblock.allocate(Dblocks*blocksize, Dblocks*blocksize, nonzeroes);

                MPI_Recv(&(Dblock.pRows[0]), Dblocks * blocksize + 1, MPI_INT,    i, i + size,     MPI_COMM_WORLD, &status);
                MPI_Recv(&(Dblock.pCols[0]), nonzeroes,               MPI_INT,    i, i + 2 * size, MPI_COMM_WORLD, &status);
                MPI_Recv(&(Dblock.pData[0]), nonzeroes,               MPI_DOUBLE, i, i + 3 * size, MPI_COMM_WORLD, &status);

                // MPI_Get_count(&status, MPI_DOUBLE, &count); printf("Process 0 received %d elements of process %d\n",count,i);
                // MPI_Get_count(&status, MPI_INT, &count);    printf("Process 0 received %d elements of process %d\n",count,i);
                // MPI_Get_count(&status, MPI_INT, &count);    printf("Process 0 received %d elements of process %d\n",count,i);

                Dmat.addBCSR(Dblock);
                Dblock.clear();
            }
        }

        //Dmat.writeToFile("D_sparse.csr");
        
        Btsparse.transposeIt(1);
        
        Dmat.reduceSymmetric();
        Dmat.nrows = dimD;
        Dmat.ncols = dimD;
        Dmat.pRows = (int*) realloc(Dmat.pRows, (dimD+1)*sizeof(int));
        
        solver->create2x2SymBlockMatrix(A, Btsparse, Dmat, Csparse);
        

        string filenameC = options->output_prefix;

        if (filenameC.at(filenameC.size()-1) != '/')
            filenameC += "/";

        filenameC += "debug_C.csr";

        Csparse.writeToFile(filenameC.c_str());


        Btsparse.clear();
         Csparse.clear();
            Dmat.clear();
        
        /* 
        if (filenameC != NULL)
            free(filenameC.c_str());
        filenameC = NULL;
        */
        

    }

    callBlacsBarrier();

}



#ifndef Solver_hpp
#define Solver_hpp


#include "CSRdouble.hpp"
#include "ParDiSO.hpp"
#include "timing.hpp"
#include "shared_var.h"
#include <cassert>
#include <cstring>

class Solver
{
    public:

    private:
        int  iam;
        int  ICTXT2D;
        bool isBsparse;

    public:

        ~Solver();
        Solver(int process_rank, int context, bool isBsparse_bool);

        void solveSystem(CSRdouble& A, double* X, double* B, int pardiso_mtype, int number_of_rhs);
        void solveSystem(CSRdouble& A, double* X, double* B, int pardiso_mtype, int number_of_rhs, ParDiSO& solver);
        void calculateSchurComplement(CSRdouble& A, int pardiso_mtype, CSRdouble& S);
        void errorReport(int number_of_rhs, CSRdouble& A, double* x, double* b);
        void generateRhs(int number_of_rhs, int n, double* B);
        void writeSolution(int number_of_solution_vectors, int n, double* X);


        void create1x2BlockMatrix(CSRdouble& A, CSRdouble& B, CSRdouble& C);
        void create2x2BlockMatrix(CSRdouble& A, CSRdouble& B, CSRdouble& C, CSRdouble& D, CSRdouble& W);

        void create2x2SymBlockMatrix(CSRdouble& A, CSRdouble& B, CSRdouble& T, CSRdouble& C);
        void create2x2SymBlockMatrix_denseT(CSRdouble& A, CSRdouble& B, double* T, CSRdouble& C);
        void create2x2BlockMatrix_denseT_lldT(CSRdouble& A, CSRdouble& B_i_T, CSRdouble& B_j, double* T, int lld_T, CSRdouble& C);

        void createAugmentedMatrix(CSRdouble& A, CSRdouble& B, CSRdouble& C, CSRdouble& K);
        void make_Sij_sparse_parallel(CSRdouble& A, CSRdouble& BT_i, CSRdouble& B_j, double * T_ij, int lld_T);
        void make_Sij_parallel_denseB(CSRdouble& A, CSRdouble& BT_i, CSRdouble& B_j, double * T_ij, int lld_T, double * AB_sol_out);
        void solveManyRhsUsingSchurComplement(CSRdouble& A, int nrhs, int pardiso_mtype);
        void makeRhsCSRdoubleMatrix(int number_of_rhs, int n, double* B, CSRdouble& Bmat);
        void makeIdentity(int n, CSRdouble& I);
        void fillSymmetric(CSRdouble* pmatrix);
        void norm(int n, double* x, double* y);


};

#endif


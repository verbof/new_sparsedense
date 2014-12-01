#ifndef SelInvProcess_hpp
#define SelInvProcess_hpp

#include "CSRdouble.hpp"
#include "Matrix2D.hpp"
#include "ParDiSO.hpp"
#include "Solver.hpp"

class SelInvProcess
{
    public:
        int M, N;
        int I, J;
        int iam;
        int size;
        int ICTXT2D;



    private:
        int       Dblocks;
        int       lld_D;
        int       DESCD;
        int       DESCAB_sol;
        int       DESCXSROW;

        int       DLEN_;

        ParDiSO   pardiso_var;
        CSRdouble A;
        CSRdouble BT_i;
        CSRdouble B_j;
        CSRdouble Btsparse;

        Solver solver;

        Matrix2D<double> D_ij;
        Matrix2D<double> AB_sol;
        Matrix2D<double> InvD_T_Block;
        Matrix2D<double> XSrow;


    public:
        ~SelInvProcess();
        SelInvProcess(int argc, char*** argv);

        void initialize(int i, int j);
        void callBlacsBarrier();
        void readMatrices(int blocksize, int dims, int Ddim, const char* fileA, const char* fileB, const char* fileD);
        void makeSchur();
        void factorizeSchur();
        void extractDiagInvD();
        void invertA();
        void computeDiagInvA();
        void saveDiagonal(string output_prefix);
        void finalizeMPIandBLACS();
        void printDebugInfo(int blocksize);
};

#endif


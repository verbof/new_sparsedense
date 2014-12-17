#ifndef SelInvProcess_hpp
#define SelInvProcess_hpp

#include "CSRdouble.hpp"
#include "Matrix2D.hpp"
#include "ParDiSO.hpp"
#include "Options.hpp"
#include "Solver.hpp"
#include "RealMath.hpp"
#include "shared_var.h"
#include <cassert>




class SelInvProcess
{
    public:
        int M, N;
        int I, J;
        int iam;
        int size;
        int ICTXT2D;



    private:

        static const int DLEN_ = 9;

        int*      DESCD;
        int*      DESCXSROW;
        int*      DESCAB_sol;
        int       Drows;
        int       Dcols;
        int       Dblocks;
        int       lld_D;

        int       info;

        ParDiSO   pardiso_var;
        Solver*   solver;
        Options*  options;

        CSRdouble A;
        CSRdouble BT_i;
        CSRdouble B_j;
        CSRdouble Btsparse;

        Matrix2D<double> D_ij;
        Matrix2D<double> AB_sol;
        Matrix2D<double> XSrow;

        double* InvD_T_Block;

    public:
        ~SelInvProcess();
        SelInvProcess(int argc, char*** argv);

        void initialize(int i, int j);
        void callBlacsBarrier();
        void readMatrices();
        void read_in_BD(double* Dmat);
        void makeSchur();
        void factorizeSchur();
        void invertSchur();
        void extractDiagInvD();
        void invertA();
        void computeDiagInvA();
        void saveDiagonal();
        void finalizeMPIandBLACS();
        void printDebugInfo();
};

#endif


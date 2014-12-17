#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include <algorithm>

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


void computeDiagonalUsingMPI(int* argc, char*** argv)
{

    SelInvProcess process(*argc, argv);

    Options& options = Options::reference();

    process.readMatrices();

    if (options.printDebug)
    {
        process.printDebugInfo();
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

    process.saveDiagonal();


}

void computeDiagonalUsingParDiSO()
{
    // Nothing, at the moment...
}




int read_input ( string filename, Options& options ) {
    std::ifstream inputfile ( filename.c_str() );
    string line;

    int lambda=0;
    int blocksize=64;


    while (std::getline(inputfile, line)) 
    {

        if ( line == "#DimensionD" ) 
        {
            std::getline ( inputfile,line );
            options.dimD = atoi(line.c_str());
        }
        else if ( line == "#DimensionA" ) 
        {
            std::getline ( inputfile,line );
            options.dimA = atoi ( line.c_str() );
        } 
        else if ( line == "#SparseFileA" ) 
        {
            std::getline ( inputfile,line );
            options.fileA = line;
        } 
        else if ( line == "#DenseFileD" ) 
        {
            std::getline ( inputfile,line );
            options.fileD = line;
        } 
        else if ( line == "#SparseFileB" ) 
        {
            std::getline ( inputfile,line );
            options.fileB = line;
        } 
        else if ( line == "#BlockSize" ) 
        {
            std::getline ( inputfile,line );
            options.blocksize = atoi ( line.c_str() );
        } 
        else if ( line == "#Bsparse" ) 
        {
            std::getline ( inputfile,line );
            options.isBsparse = (bool) atoi ( line.c_str() );
        } 
        else if ( line == "#OutputFileSparseC" ) 
        {
            std::getline ( inputfile,line );
        } 
        else if ( line[0] == '/' || line.size() == 0 ) 
        {
            // Nothing!
        }
        else 
        {
            printf ( "Unknown parameter in inputfile, the following line was ignored: \n" );
            printf ( "%s\n",line.c_str() );
        }
    }
    return 0;
}




int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "usage: " << argv[0] << "mode MPI_grid_size inputfilename" << endl;
        cout << "examp: " << argv[0] << "MPI  5 simplecase" << endl;
        cout << "examp: " << argv[0] << "PARDISO 4 simplecase" << endl;

        exit(1);
    }

    string mode = argv[1];


    // convert mode to uppercase
    std::transform(mode.begin(), mode.end(), mode.begin(), (int(*)(int))toupper);

    const char* inputfile = argv[3];

    cout << "Input file: " << inputfile << endl;

    Options& options = Options::reference();
    //Keywords keywords(inputfile, &options);
    //keywords.read();

    int a = read_input(inputfile, options);

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


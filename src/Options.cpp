// include files
// =============
#include "config.hpp"
#include "Options.hpp"
// =============

Options* Options::pOptions_ = 0;

// keyword MATRIX_SIZES
  int    Options::dimD;
  int    Options::dimA;
 
  // keyword MATRIX_FILES
  string Options::fileA;
  string Options::fileB;
  string Options::fileD;

  // keyword B_MATRIX_TYPE
  bool Options::isBsparse;

  // keyword BLOCK_SIZE
  int  Options::blocksize;

  // keyword OUTPUT
  string Options::output_prefix;

  // keyword DEBUG
  bool Options::printDebug;

Options& Options::reference()
{
  if (pOptions_ == 0)
  {
    static Options global_options;
    pOptions_ = &global_options;
    return global_options;
  }
  else 
  {
    return *pOptions_;
  }
  
}




void Options::print()
{
  cout << "\n\n\n";
  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  cout << "@@@                         PARAMETERS                        @@@\n";
  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  
  cout << "\n";
  
  cout << "dimA             #" << setw(15) << dimA << std::endl;
  cout << "dimD             #" << setw(15) << dimD << std::endl;
  cout << endl;
  cout << endl;

  cout << "fileA             #" << setw(15) << fileA << std::endl;
  cout << "fileB             #" << setw(15) << fileB << std::endl;
  cout << "fileD             #" << setw(15) << fileD << std::endl;
  cout << endl;
  cout << endl;

  cout << "isBsparse         #" << setw(15) << isBsparse  << std::endl;
  cout << "printDebug        #" << setw(15) << printDebug << std::endl;
  cout << endl;
  cout << endl;

  cout << "blocksize         #" << setw(15) << blocksize << "\n";
  cout << "output_prefix     #" << setw(15) << output_prefix << "\n";
  cout << endl;
  cout << endl;
  
  cout << "\n";
  cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
  cout << "\n\n\n";
}



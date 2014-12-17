// include files
// =================
#include "config.hpp"
#include "Keywords.hpp"
#include "Options.hpp"


Keywords::Keywords(const char* filename, Options* poptions)
{
  pOptions_ = poptions;

  init(filename);
}



void Keywords::read()
{
  KeywordsBase::read();
}


bool Keywords::readKeyword_(const string& keyword)
{
  if (keyword == "MATRIX_SIZES")
    MATRIX_SIZESKeyword();
  else if (keyword == "MATRIX_FILES")
    MATRIX_FILESKeyword();
  else if (keyword == "B_MATRIX_TYPE")
    B_MATRIX_TYPEKeyword();
  else if (keyword == "BLOCK_SIZE")
    BLOCK_SIZEKeyword();
  else if (keyword == "OUTPUT")
    OUTPUTKeyword();
  else if (keyword == "DEBUG")
    DEBUGKeyword();
  else
    return false;

  return true;
}



void Keywords::MATRIX_SIZESKeyword()
{
  parse >> Options::dimA >> Options::dimD;
}

void Keywords::MATRIX_FILESKeyword()
{
  parse >> Options::fileA >> Options::fileB >> Options::fileD;
}

void Keywords::B_MATRIX_TYPEKeyword()
{
  string btype;
  parse >> btype;

  std::transform(btype.begin(), btype.end(), btype.begin(), (int(*)(int))toupper);   
  Options::isBsparse = (btype == "YES");

}

void Keywords::BLOCK_SIZEKeyword()
{
  parse >> Options::blocksize;
}


void Keywords::DEBUGKeyword()
{
  string debugout;
  parse >> debugout;
  
  std::transform(debugout.begin(), debugout.end(), debugout.begin(), (int(*)(int))toupper);   
  Options::printDebug = (debugout == "YES");

}

void Keywords::OUTPUTKeyword()
{
  parse >> Options::output_prefix;
}

#ifndef Options_hpp

#define Options_hpp

// include files
// =============
#include "config.hpp"
#include "Matrix2D.hpp"
// =============


class Options
{
private:
  static Options* pOptions_;
  Options() {}
  ~Options() {}
  Options(const Options&);
  Options& operator=(const Options&);
  
public:
 // keyword MATRIX_SIZES
  static int    dimD;
  static int    dimA;
 
  // keyword MATRIX_FILES
  static string fileA;
  static string fileB;
  static string fileD;

  // keyword B_MATRIX_TYPE
  static bool isBsparse;

  // keyword BLOCK_SIZE
  static int  blocksize;


// keyword OUTPUT
  static bool   printDebug;
  static string output_prefix;


public:
  static Options& reference();
  static void   print();
};

#endif


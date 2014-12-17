#ifndef Keywords_hpp
#define Keywords_hpp

#include "KeywordsBase.hpp"
#include <algorithm>

class Options;

/// class for handling all keywords in the input file
class Keywords : public KeywordsBase
{
  public:
    ~Keywords(){};
    Keywords(const char* filename, Options* poptions);
    void read();

    void MATRIX_SIZESKeyword();
    void MATRIX_FILESKeyword();
    void B_MATRIX_TYPEKeyword();
    void BLOCK_SIZEKeyword();
    void OUTPUTKeyword();
    void DEBUGKeyword();
  
  private:
    Options* pOptions_;
    bool readKeyword_(const string& keyword);
};


#endif 


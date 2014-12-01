#ifndef Parse_hpp
#define Parse_hpp

// include files
// =============
#include "config.hpp"
#include "Matrix2D.hpp"
// =============

class Parse
{
private:
  fstream* pFile_;

public:
  ~Parse();
  Parse();

  void setCurrentStream(fstream* pstream);

  template<typename T>
  Parse& operator>>(T& s);

  Parse& operator>>(vector<int>& v);
  Parse& operator>>(vector<size_t>& v);
  Parse& operator>>(vector<double>& v);
  Parse& operator>>(vector<string>& v);
  Parse& operator>>(vector<vector<string> >& vv);
  Parse& operator>>(Matrix2D<double>& A);
  
  void omitEmptySpace(fstream& fin);
  void eatCommentsAndSpaces(fstream& fin);

  void skipComments();
  void skipCommentsAndSpaces();
  void readAWholeLine(string& line);
  void readUntilEndCharacter(string& line);
  void readUntilEndCharacterMatrices(string& line);
  void readTableOfData(vector<string>& vlines); 
  void readTableOfDataMatrices(vector<string>& vlines); 

  void tokenizeString(const string& str, vector<string>& tokens, const string& delims=", \t");
  void expand(string& line, vector<string>& ventries);
};



template<typename T>
Parse& Parse::operator>>(T& s)
{
  *pFile_ >> s;
      
  return *this;
}


#endif

#ifndef KEYWORDS_BASE_HPP
#define KEYWORDS_BASE_HPP

// include files
// =============
#include "Parse.hpp"
#include <stack>
#include <string>

using std::stack;
using std::string;
// =============

/// abstract class for Keywords data reader
class KeywordsBase
{
public:
  ~KeywordsBase();
  virtual void read();
  static std::string parseRelativePath(const char* filename);
  static std::string parseBaseFileNameWithExt(const char *filename);  

public:
  string               relative_path;

protected:
  KeywordsBase(){}

  virtual void init(const char* filename);

private:
  virtual void readKeyword(const string& keyword); 
  virtual bool readKeyword_(); 
  void INCLUDEkeyword();

protected:
  string               currentFileName_;
  Parse                parse;
  stack<fstream*>      filesToBeParsed_;
};



#endif // KEYWORDS_BASE_HPP

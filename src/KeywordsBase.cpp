// include files
// =================
#include "KeywordsBase.hpp"
#include <sstream>
#include <algorithm>

KeywordsBase::~KeywordsBase()
{
  while (!filesToBeParsed_.empty())
  {
    delete filesToBeParsed_.top();
    filesToBeParsed_.pop();
  }
}

void KeywordsBase::init(const char* filename)
{
  fstream* pfile = new fstream(filename, ios::in);
  if (!pfile->is_open())
  {
    delete pfile;
    cout << "Could not open file " << filename << " for input!" << std::endl;
    exit(-1);
  }
  filesToBeParsed_.push(pfile);
  relative_path = parseRelativePath(filename);
}



std::string KeywordsBase::parseRelativePath(const char* filename)
{
  string rel_path_temp = filename;
  int rel_path_pos = -1;
  for (int i = rel_path_temp.length() - 1; i >= 0; i--)
  {
    if (rel_path_temp[i] == '\\' || rel_path_temp[i] == '/')
    {
      rel_path_pos = i;
      break;
    }
  }
  if (rel_path_pos < 0)
    return "";
  rel_path_temp.resize(rel_path_pos+1);
  return rel_path_temp;
}


/*! \brief get base name from full path of file
  \example result of parseBaseFileNameWithExt("c:\\work\\1.txt") is "1.txt"   
  \param [in] filename - name of file for parsing*/      
std::string KeywordsBase::parseBaseFileNameWithExt(const char *filename)
{
  string file_name_with_ext = filename;
  int rel_path_pos = -1;
  for (int i = file_name_with_ext.length() - 1; i >= 0; i--)
  {
    if (file_name_with_ext[i] == '\\' || file_name_with_ext[i] == '/')
    {
      rel_path_pos = i;
      break;
    }
  }
  if (rel_path_pos < 0)
    return "";
  file_name_with_ext = file_name_with_ext.substr(rel_path_pos+1);
  return file_name_with_ext;
}




void KeywordsBase::read()
{
  string oneline;
  string string_read;

  cout << "##################\n";
  while (!filesToBeParsed_.empty())
  {
    fstream* pcurrent_file = filesToBeParsed_.top();
    parse.setCurrentStream(pcurrent_file);

    parse.eatCommentsAndSpaces(*pcurrent_file);
    parse.readAWholeLine(oneline);
    string_read.clear(); 

    // remove spaces before the keyword and after if any exist
    for (int i = 0; i < (int)oneline.length(); i++)
    {
      if (oneline[i] != ' ')
        string_read += oneline[i];
    }
    parse.eatCommentsAndSpaces(*pcurrent_file);

    if (string_read == "END" || string_read == "end" || pcurrent_file->eof())
    {
      delete pcurrent_file;
      filesToBeParsed_.pop();

      cout << "Finished with file\n";
      cout << "##################\n";
      // cout << filesToBeParsed_.empty() << "\n";
    }
    else 
      readKeyword(string_read);
  }
}


void KeywordsBase::readKeyword(const string& kwd)
{
  string keyword = kwd;
  std::transform(keyword.begin(), keyword.end(), keyword.begin(), (int(*)(int))toupper);   

  if (!keyword.empty())
  {
    cout << "reading keyword ..." << keyword << "...\n";

    if (keyword == "INCLUDE")
      INCLUDEkeyword();

    else if( !readKeyword_(keyword) )
      cout << "****** keyword not found ******\n";
  }
}

void KeywordsBase::INCLUDEkeyword()
{
  string filename_temp, filename;
  parse >> filename_temp;

  filename         = relative_path + filename_temp;
  currentFileName_ = filename;

  cout << "filename = " << filename << "\n";
  fstream* pnew_file = new fstream(filename.c_str(), ios::in);
  if (pnew_file->fail())
  {
    cout << "failed to open file: " << filename << "\n";
    exit(1);
  }
  filesToBeParsed_.push(pnew_file);
}



bool KeywordsBase::readKeyword_(const string& keyword)
{
  return true;
}



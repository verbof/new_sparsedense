// include files
// =================
#include "Parse.hpp"
// =================
Parse::~Parse()
{}

Parse::Parse()
{
  pFile_ = 0;
}


void Parse::setCurrentStream(fstream* pstream)
{
  pFile_ = pstream;
}


void 
Parse::omitEmptySpace(fstream& fin) 
{
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // @ omitEmptySpace_ removes all kinds of empty space characters from the @
  // @ input stream (spaces, tabs and new lines) from the input stream until@
  // @ it finds a character which is non of the above.                      @
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // * 
  // * Parameters
  // * ----------
  // * fin       the stream from which empty space will be removed.
  // *
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    

  char temp;
  do {
    temp = fin.get();
    if (fin.eof())
    {
      // cout << "Omit Empty space EOF FOUND\n";
      return;
    }
  } while ( temp == ' ' || 
            temp == '\t' || 
            temp == '\n' || 
            temp == char(10) || 
            temp == char(13));

  fin.putback(temp);
}


void 
Parse::eatCommentsAndSpaces(fstream& fin)
{
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // @ eatCommentsAndSpaces_ removes all kinds of empty space characters from@
  // @ the input stream (spaces, tabs and new lines) as well as any comments @
  // @ until it finds a character which is non of the above.                 @
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  // *
  // * Any line which starts by one of the following characters is 
  // * considered to be a comment and it is removed from the input steram.
  // * Characters indicating comments may be: #, - /
  // * 
  // * Parameters
  // * ----------
  // *
  // * fin       the stream from which empty space and comments will be
  // *           removed.
  // *
  // *
  // * Class methods used
  // * ------------------
  // *
  // * Parse::omitEmptySpace_(fstream&)    remove all kinds of empty space
  // * 
  // @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  char character;
  char second_character;
  bool foundcomment;
  bool filefinished;

  do
  {
    // omit all spaces
    omitEmptySpace(fin);

    character = fin.peek();

    fin.get();
    if (fin.eof())
    {
      // cout << "EAT COMMENTS AND SPACES EOF FOUND\n";
      return;
    }

    // try second character 
    second_character = fin.peek();
    
    // put the character back
    fin.putback(character);

    filefinished = fin.eof();

    foundcomment = (character == '#'  || 
                    (character == '-' && second_character == '-')) && !filefinished;

    if (foundcomment)
    {
      // if this line is a comment then read it until you find enter
      bool finished;
      do
      {
        character=fin.get();

        if (fin.eof())
        {
          cout << "IN LOOP OMIT COMMENTS AND SPACES EOF FOUND\n";
        }

        finished =  character == '\n'
          || character == char(10)
          || character == char(13);

      } while (!finished);
    }

    // cout << character;		       
  } while( foundcomment  );
}


void Parse::skipComments()
{
  string oneline;
  char temp = pFile_->peek();
  while(temp == '#' || temp == '\n' || temp == char(10) || temp == char(13))
  {
    readAWholeLine(oneline);
    temp = pFile_->peek();
  }
}

void Parse::readAWholeLine(string& line)
{
  line.clear();
  char ch  = pFile_->get();

  while( ch != '\n' && ch != char(10) && ch != char(13) && ch != '/' && ch != char(-1) 
         && line != "end" && line != "END" )
  {
    line += ch;
    ch  = pFile_->get();
    if (pFile_->eof())
    {
      cout << "HERE EOF FOUND\n";
    }

  }
}


void Parse::readUntilEndCharacter(string& line)
{
  line.clear();
  char ch  = pFile_->get();
  bool quotes_havent_been_opened = true;
  char character;
  while( ch != '/' || !quotes_havent_been_opened )
  {
    if ( !(ch == '\n' || ch == char(10) || ch == char(13)) )
    {
      if (ch == '\"' || ch == '\'')
      {
        quotes_havent_been_opened = !quotes_havent_been_opened;
      }
      else
      {
        if (ch == '#'  || (ch == '-' && pFile_->peek() == '-'))
        {
          // if this line is a comment then read it until you find enter
          bool finished;
          do
          {
            character = pFile_->get();
            finished =  character == '\n'
            || character == char(10)
            || character == char(13);

          } while (!finished);
          line += " ";
        }
        else
          line += ch;
      }
    }
    else
      line += " ";

    ch  = pFile_->get();
  }
}



void Parse::readUntilEndCharacterMatrices(string& line)
{
  line.clear();
  char ch  = pFile_->get();
  bool quotes_havent_been_opened = true;
  char character;
  while( (ch != '/' && ch != '\n') || !quotes_havent_been_opened )
  {
    if ( !(ch == '\n' || ch == char(10) || ch == char(13)) )
    {
      if (ch == '\"' || ch == '\'')
      {
        quotes_havent_been_opened = !quotes_havent_been_opened;
      }
      else
      {
        if (ch == '#'  || (ch == '-' && pFile_->peek() == '-'))
        {
          // if this line is a comment then read it until you find enter
          bool finished;
          do
          {
            character = pFile_->get();
            finished =  character == '\n'
                     || character == char(10)
                     || character == char(13);

          } while (!finished);
          line += " ";
        }
        else
          line += ch;
      }
    }
    else
      line += " ";

    ch  = pFile_->get();
  }
}

// COMPDAT
// 1 2 3 4 5 str/
// 4 5 6 7 8 upr/
//                /

void Parse::readTableOfData(vector<string>& vlines)
{
  bool finished;
  vlines.clear();

  char ch  = pFile_->peek();
  finished = ch == '/';
  if(finished) pFile_->get();	// to eat '/'

  string line, temp;
  while( !finished )
  {    
    line.clear();
    // getline (*pFile_, line, '/');
    readUntilEndCharacter(line);
    eatCommentsAndSpaces(*pFile_);
    int spaces = 0;
    int line_length = line.size();
    for (int i = 0; i < line_length; i++)
    {
      if (line[i] == ' ' || line[i] == '\t')
        spaces++;
      else
        break;
    }

    line.erase (line.begin(), line.begin() + spaces);

    finished = line.empty();

    if (!finished)
      vlines.push_back(line);
  }
}


void Parse::readTableOfDataMatrices(vector<string>& vlines)
{
  bool finished;
  vlines.clear();

  char ch  = pFile_->peek();
  finished = ch == '/';
  if(finished) pFile_->get();	// to eat '/'

  string line, temp;
  while( !finished )
  {    
    line.clear();
    // getline (*pFile_, line, '/');
    readUntilEndCharacterMatrices(line);
    eatCommentsAndSpaces(*pFile_);
    int spaces = 0;
    int line_length = line.size();
    for (int i = 0; i < line_length; i++)
    {
      if (line[i] == ' ' || line[i] == '\t')
        spaces++;
      else
        break;
    }

    line.erase (line.begin(), line.begin() + spaces);

    finished = line.empty();

    if (!finished)
      vlines.push_back(line);
  }
}



/// \brief convert input string into vector of string tokens
///
/// \note consecutive delimiters will be treated as single delimiter
/// \note delimiters are _not_ included in return data
///
/// \param input string to be parsed
/// \param delims list of delimiters.

void Parse::tokenizeString(const string& str, vector<string>& tokens, const string& delims)
{
  // output vector
  tokens.clear();

  // Skip delims at beginning, find start of first token
  string::size_type lastPos = str.find_first_not_of(delims, 0);

  while (lastPos != string::npos)  // check lastPos is enough (if lastPos is string::npos, so is pos)
  {
    // Find next delimiter at end of token
    const string::size_type pos = str.find_first_of(delims, lastPos);

    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // Skip delims.  Note the "not_of". this is beginning of token
    lastPos = str.find_first_not_of(delims, pos);
  }

  // In this case we have find something like * or 5* and we should replace it with * 
  // or 5 * (since 5 has already been pushed, we only need to push *)
  if (delims == "*" && str[str.length() - 1] == '*')
    tokens.push_back("*");
}


void Parse::expand(string& line, vector<string>& ventries)
{
  // this will be the expanded array of strings
  ventries.clear();

  vector<string> vspace_tokens, vnumbers;
  // this reads the input string line and expands its to tokens
  // seperated with ',' or ' ', or tab (by default)
  tokenizeString(line, vspace_tokens);
  ventries.reserve(vspace_tokens.size());

  // we want to search now each of these tokens because some of them
  // may be of the form 5*1234 which means that 5 of them follow with
  // value 1234
  for (size_t i = 0; i < vspace_tokens.size(); i++)
  {
    // each token in vspace_tokens is splitted to two numbers if a
    // character '*' exists or a single number if no character '*'
    // exists
    tokenizeString(vspace_tokens[i], vnumbers, "*");

    if (vnumbers.size() == 1)
      ventries.push_back(vnumbers[0]);
    else if (vnumbers.size() == 2)
    {
      // here it has found a token of the form number1*number2 which
      // is supposed to mean, number1 numbers equal to number2
      int number_of_identical_entries = int( strtol(vnumbers[0].c_str(), NULL, 0) );
      for (int j = 0; j < number_of_identical_entries; j++)
        ventries.push_back(vnumbers[1]);
    }
  }
}




Parse& Parse::operator>>(vector<double>& v)
{
  eatCommentsAndSpaces(*pFile_);
  string line;
 
  readUntilEndCharacter(line);

  vector<string> tokens;
  expand(line, tokens);

  if(v.size() < tokens.size())
    v.resize(tokens.size());
  for (size_t i = 0; i < tokens.size(); i++)
    v[i] = strtod(tokens[i].c_str(), NULL);

  return *this;
}


Parse& Parse::operator>>(vector<int>& v)
{
  eatCommentsAndSpaces(*pFile_);
  string line;

  readUntilEndCharacter(line);
  vector<string> tokens;
  expand(line, tokens);

  if(v.size() < tokens.size())
    v.resize(tokens.size());

  for (size_t i = 0; i < tokens.size(); i++)
    v[i] = int( strtol(tokens[i].c_str(), NULL, 0) );

  return *this;
}



Parse& Parse::operator>>(vector<size_t>& v)
{
  eatCommentsAndSpaces(*pFile_);
  string line;

  readUntilEndCharacter(line);
  vector<string> tokens;
  expand(line, tokens);

  if(v.size() < tokens.size())
    v.resize(tokens.size());

  for (size_t i = 0; i < tokens.size(); i++)
    v[i] = size_t( strtol(tokens[i].c_str(), NULL, 0) );

  return *this;
}



Parse& Parse::operator>>(vector<string>& v)
{
  eatCommentsAndSpaces(*pFile_);
  string line;
  readUntilEndCharacter(line);
  expand(line, v);

  return *this;
}




Parse& Parse::operator>>(vector<vector<string> >& vv)
{
  vector<string> vline;
  
  eatCommentsAndSpaces(*pFile_);
  readTableOfData(vline);
  int cnt = vline.size();
  vv.resize(cnt);

  for (int i = 0; i < cnt; i++)
    expand(vline[i], vv[i]);

  return *this;
}


/*
Parse& Parse::operator>>(Matrix2D<double>& A)
{
  vector<string> vline;
  vector<string> v;

  eatCommentsAndSpaces(*pFile_);
  readTableOfDataMatrices(vline);
  int cnt = vline.size();

  for (int i = 0; i < cnt; i++)
  {
    expand(vline[i], v);
    for (unsigned int j = 0; j < v.size(); j++)
    {
      A.data.push_back(strtod(v[j].c_str(), NULL));
    }
  }
 
  A.ncols = vline.size();
  A.nrows = v.size();

  return *this;
}
*/



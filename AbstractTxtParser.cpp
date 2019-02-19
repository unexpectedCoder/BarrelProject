#include "AbstractTxtParser.h"

#include <stdlib.h>
#include <sstream>

using namespace std;


AbstractTXTParser::AbstractTXTParser()
{
  curr_str = "...eof...\n";
  digit = {{'0', 0},
           {'1', 1},
           {'2', 2},
           {'3', 3},
           {'4', 4},
           {'5', 5},
           {'6', 6},
           {'7', 7},
           {'8', 8},
           {'9', 9}};
}

AbstractTXTParser::AbstractTXTParser(const string &path, char mode)
{
  curr_str = "...eof...\n";
  open(path, mode);
  digit = {{'0', 0},
           {'1', 1},
           {'2', 2},
           {'3', 3},
           {'4', 4},
           {'5', 5},
           {'6', 6},
           {'7', 7},
           {'8', 8},
           {'9', 9}};
}

void AbstractTXTParser::open(const string &path, char mode)
{
  if (!isExisting(path))
    throw "Error: file is not found or existing!";
  if (!checkMode(mode))
    throw "Error: unknown file opening mode!";

  if (mode == 'w')
  {
    file.open(path, ios_base::out);
    return;
  }
  if (mode == 'r')
  {
    file.open(path, ios_base::in);
    return;
  }
  if (mode == 'a')
  {
    file.open(path, ios_base::app);
    return;
  }
}

bool AbstractTXTParser::isExisting(const string &path) const
{
  fstream f(path.c_str());
  if (f.is_open())
    return true;
  return false;
}

void AbstractTXTParser::close()
{
  if (isOpen())
    file.close();
}

bool AbstractTXTParser::isOpen() const
{
  if (file.is_open())
    return true;
  return false;
}

bool AbstractTXTParser::isEnd() const
{
  if (file.eof())
    return true;
  return false;
}

void AbstractTXTParser::createFileTXT(const string &path)
{
  ofstream f{ path.c_str() };
}

void AbstractTXTParser::write(const string &text)
{
  if (!isOpen())
    throw "Error: file is not opened!";
  if (!text.empty())
    file << text;
}

void AbstractTXTParser::write(double v)
{
  if (!isOpen())
    throw "Error: file is not opened!";
  file << v;
}

void AbstractTXTParser::write(char split, double v)
{
	if (!isOpen())
		throw "Error: file is not opened!";
	file << split << v;
}

void AbstractTXTParser::write(double v, char split)
{
  if (!isOpen())
    throw "Error: file is not opened!";
  file << v << split;
}

void AbstractTXTParser::write(double v, const string &text)
{
  if (!isOpen())
    throw "Error: file is not opened!";
  file << v << text;
}

void AbstractTXTParser::write(const string &text, double v, char split)
{
  if (!isOpen())
    throw "Error: file is not opened!";
  file << text << v << split;
}

void AbstractTXTParser::write(const string &pre_text, double v, const string &post_text)
{
  if (!isOpen())
    throw "Error: file is not opened!";
  file << pre_text << v << post_text;
}

const string& AbstractTXTParser::readStr()
{
  if (!isOpen())
    throw "Error: file is not opened!";

	char buf[MAX_CHAR_LINE];
  if (!file.eof())
    file >> buf;
	curr_str = buf;
  return curr_str;
}

double AbstractTXTParser::readDouble()
{
  if (!isOpen())
    throw "Error: file is not opened!";

  double v;
  file >> v;
  return v;
}

double** AbstractTXTParser::readMatrix(int &r, int &c)
{
  vector<double*> buf;
  char line[MAX_CHAR_LINE];
  int c_buf = 0;
  r = 0; c = 0;

  while (true)
  {
    if (c != c_buf)
      throw "Error: invalid matrix in file! Number of columns doesn't match in rows!";

    getLine(line, MAX_CHAR_LINE);
    if (isEnd() || strcmp(line, "") == 0) break;
    double *arr = readArray(line, c);
    c_buf = c;
    buf.push_back(arr);
  }

  r = buf.size();
  if (r < 1 || c < 1)
    throw "Error: invalid matrix size! Number of rows or columns < 1!";
  matr = new double*[r];
  for (int i = 0; i < r; i++)
    matr[i] = buf[i];
  return matr;
}

double* AbstractTXTParser::readArray(const char *line, int &size)
{
  if (isEnd() || !isOpen())
    throw "Error: end of file or file is not opened!";
  vector<double> buf;
  return strToArray(line, size);
}


void AbstractTXTParser::getLine(char *s, int size, char split)
{
  if (isEnd() || !isOpen())
    throw"Error: end of file or file is not opened!";
  file.getline(s, size, split);
}

void AbstractTXTParser::newLine()
{
  if (!isOpen())
    throw "Error: file is not opened!";
  file << '\n';
}


bool AbstractTXTParser::checkMode(char mode)
{
  if (mode == 'w' || mode == 'r' || mode == 'a')
    return true;
  return false;
}

double AbstractTXTParser::strToDouble(const string &str)
{
  const char* s = str.c_str();
  return strToDouble(s);
}

double AbstractTXTParser::strToDouble(const char *s)
{
  string val_str;
  string_num(s, val_str);

  double val = 0;
  int i;
  int pos_neg = 1;
  if (val_str[0] == '-')
    pos_neg = -1;

  bool is_integer = true;
  double ks = 1;
  int kb = 1;

  for ((val_str[0] == '-')? i = 1 : i = 0; i < val_str.size(); i++)
  {
    if (val_str[i] == '.')
    {
      is_integer = false;
      i++;
    }

    if (is_integer)
    {
      kb *= 10;
      integer(val, kb, digit[val_str[i]]);
    }
    else
    {
      ks *= 0.1;
      fraction(val, ks, digit[val_str[i]]);
    }
  }

  return pos_neg * val;
}

double* AbstractTXTParser::strToArray(const char *s, int &size)
{
  string val_str = "";
  vector<double> buf;
  while (true)
  {
    if (isalpha(*s) || isspace(*s) || !*s)
    {
      if (!val_str.empty())
      {
        buf.push_back(strToDouble(val_str));
        val_str.clear();
      }
      if (!*s)
        break;
      s++;
      continue;
    }
    val_str += *s;
    s++;
  }

  size = buf.size();
  arr = new double[size];
  for (int i = 0; i < size; i++)
    arr[i] = buf[i];
  return arr;
}

void AbstractTXTParser::string_num(const char *s, string &val_str)
{
  val_str = "";
  while (*s != '\0')
  {
    if (isalpha(*s) || isspace(*s)) break;
    if (!isdigit(*s) && *s != '-' && *s != '.')
    {
      s++;
      continue;
    }
    if ((val_str.empty() || val_str == "-") && *s == '.')
    {
      val_str = "0.";
      s++;
      continue;
    }

    val_str += *s;
    s++;
  }
}

void AbstractTXTParser::integer(double &v, int k, int a)
{
  v = k * v + a;
}

void AbstractTXTParser::fraction(double &v, double k, int a)
{
  v += k * a;
}

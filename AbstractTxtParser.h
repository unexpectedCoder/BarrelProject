#ifndef TXT_PARSER
#define TXT_PARSER

/*
	* The class implements the basic functions of writing
	* and reading text files.
	* Intended for subsequent inheritance.
*/

#include <fstream>
#include <cstring>
#include <stack>
#include <map>
#include <vector>

#define MAX_CHAR_LINE 128


class AbstractTXTParser
{
public:
  AbstractTXTParser();
  AbstractTXTParser(const std::string &path, char mode);
	~AbstractTXTParser() {}

  static void createFileTXT(const std::string &path);
  void open(const std::string &path, char mode);
  void close();

  bool isExisting(const std::string &path) const;
  bool isOpen() const;
  bool isEnd() const;

  void write(const std::string &text);
  void write(double v);
	void write(char split, double v);
  void write(double v, char split);
  void write(double v, const std::string &text);
  void write(const std::string &text, double v, char split = '\n');
  void write(const std::string &pre_text, double v, const std::string &post_text);

  template<typename T>
  void writeArray(const T* arr, int size, char split = '\t');
  template<typename T>
  void writeMatrix(T** matr, int rows, int columns);

  const std::string& readStr();
  double readDouble();
  double* readArray(const char *line, int &size);
  double** readMatrix(int &r, int &c);

  double strToDouble(const std::string &str);
  double strToDouble(const char *s);
  double* strToArray(const char *s, int &size);
  void getLine(char *s, int size = MAX_CHAR_LINE, char split = '\n');
  void newLine();

protected:
  bool checkMode(char mode);
  void string_num(const char *s, std::string &str);
  void integer(double &v, int k, int a);
  void fraction(double &v, double k, int a);

protected:
  std::fstream file;
  std::string curr_str;
  std::map<char, int> digit;
  double *arr, **matr;
  int rows;
};


template<typename T>
void AbstractTXTParser::writeArray(const T* arr, int n, char split)
{
  if (n < 1)
    throw "Error: invalid array size! Array size must be > 0!";

  for (int i = 0; i < n - 1; i++)
    file << arr[i] << split;
	if (split != '\n')
		file << arr[n - 1] << "\n\n";	// Double carriage transfer so that
																	// when reading an array is not perceived as a matrix.
	else
		file << arr[n - 1] << '\n';
}

template<typename T>
void AbstractTXTParser::writeMatrix(T** m, int r, int c)
{
  if (r < 1 || c < 1)
    throw "Error: invalid matrix size! Number of rows and columns must be > 0!";

  for (int i = 0; i < r; ++i)
  {
    for (int j = 0; j < c - 1; ++j)
      file << m[i][j] << '\t';
    file << m[i][c - 1] << '\n';
  }
  file << '\n';
}

#endif

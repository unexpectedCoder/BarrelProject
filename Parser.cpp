#include "Parser.h"

using namespace std;

bool Parser::open(const string &path)
{
	if (file.is_open())
		file.close();

	file.open(path.c_str());
	if (file.is_open())
		return true;
	return false;
}

bool Parser::close()
{
	file.close();

	if (!file.eof())
		return true;
	return false;
}

double Parser::readNext()
{
	double x;
	file >> x;
	return x;
}

string Parser::readStr()
{
	string s;
	file >> s;
	return s;
}

void Parser::write(const string &txt)
{
	file << txt;
}

void Parser::write(double x, char split)
{
	file << x << split;
}

bool Parser::isEnd() const
{
	return file.eof();
}

void Parser::createFile(const string &path)
{
	ofstream f{ path.c_str() };
}

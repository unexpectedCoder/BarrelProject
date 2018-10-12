#include "Parser.h"

using namespace std;

Parser::Parser(const string &path, char mode)
{
	switch (mode)
	{
	case 'w': file.open(path.c_str(), ios_base::out);
		break;
	case 'r': file.open(path.c_str(), ios_base::in);
		break;
	case 'a': file.open(path.c_str(), ios_base::app);
		break;
	default: throw "Error: incorrect opening mode!";
	}
}

bool Parser::open(const string &path, char mode)
{
	if (file.is_open())
		file.close();

	switch (mode)
	{
	case 'w': file.open(path.c_str(), ios_base::out);
		break;
	case 'r': file.open(path.c_str(), ios_base::in);
		break;
	case 'a': file.open(path.c_str(), ios_base::app);
		break;
	default: throw "Error: incorrect opening mode!";
	}

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

Barrel& Parser::readBarrel()
{
	file >> barr.Cq >> barr.CE >> barr.CE15 >> barr.eta_omega >> barr.omega_q >>
		barr.pm_kr >> barr.pm >> barr.hi1 >>  barr.ns;

	return barr;
}

void Parser::write(const string &txt)
{
	file << txt;
}

void Parser::write(double x, char split)
{
	file << x << split;
}

void Parser::writeBarrel(const Barrel &barr)
{
	file << barr.Cq << '\n' << barr.CE << '\n' << barr.CE15 << '\n' << barr.eta_omega << '\n' << barr.omega_q << '\n' <<
		barr.pm_kr << '\n' << barr.pm << '\n' << barr.hi1 << '\n' << barr.ns;
}

bool Parser::isEnd() const
{
	return file.eof();
}

void Parser::createFile(const string &path)
{
	ofstream f{ path.c_str() };
}

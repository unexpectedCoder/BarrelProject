#include "Parser.h"

#include <sstream>

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
		barr.pm_kr >> barr.pm >> barr.hi >> barr.ns;

	return barr;
}

const Powders& Parser::readXMLPowders(const string &path)
{
	if (!powders.empty()) powders.clear();

	Powder powder;
	pugi::xml_document doc;
	pugi::xml_parse_result res = doc.load_file(path.c_str());
	if (!res) throw "Error in <readPowders()>!";

	pugi::xml_node nd, ndpow;
	nd = doc.child("root").child("POWDERS");
	ndpow = nd.child("Powder");
	while (ndpow)
	{
		powder.name = ndpow.attribute("name").as_string();
		powder.f = ndpow.child("Data").attribute("f").as_double() * 1e6;
		powder.alpha = ndpow.child("Data").attribute("alpha").as_double() * 1e-3;
		powder.delta = ndpow.child("Data").attribute("delta").as_double() * 1e3;
		powder.Ik = ndpow.child("Data").attribute("Ik").as_double() * 1e6;
		powder.k = ndpow.child("Data").attribute("k").as_double();
		powder.kappa1 = ndpow.child("Data").attribute("kappa1").as_double();
		powder.kappa2 = ndpow.child("Data").attribute("kappa2").as_double();
		powder.k_f = ndpow.child("Data").attribute("kappa_f").as_double();
		powder.k_I = ndpow.child("Data").attribute("k_f").as_double();
		powder.lambda1 = ndpow.child("Data").attribute("lambda1").as_double();
		powder.lambda2 = ndpow.child("Data").attribute("lambda2").as_double();
		powder.T = ndpow.child("Data").attribute("T").as_double();
		powder.zk = ndpow.child("Data").attribute("zk").as_double();

		powders.push_back(powder);

		ndpow = ndpow.next_sibling("Powder");
	}

	return powders;
}

void Parser::write(const string &txt)
{
	file << txt;
}

void Parser::write(const string &txt, double x, char split)
{
	file << txt << x << split;
}

void Parser::write(double x)
{
	file << x;
}

void Parser::write(char split, double x)
{
	file << split << x;
}

void Parser::write(double x, char split)
{
	file << x << split;
}

void Parser::writeBarrel(const Barrel &barr)
{
	file << barr.Cq << '\n' << barr.CE << '\n' << barr.CE15 << '\n' << barr.eta_omega << '\n' <<
		barr.omega_q << '\n' << barr.pm_kr << '\n' << barr.pm << '\n' << barr.hi << '\n' << barr.ns;
}

bool Parser::isOpen() const
{
	if (file.is_open())
		return true;
	return false;
}

bool Parser::isEnd() const
{
	return file.eof();
}

void Parser::createFile(const string &path)
{
	ofstream f{ path.c_str() };
}

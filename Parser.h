#ifndef PARSER_H
#define PARSER_H

#include <fstream>
#include <string>
#include "pugixml.hpp"

#include "AbstractTxtParser.h"
#include "Types.h"

class Parser : public AbstractTXTParser
{
private:
	Barrel barr;
	Powders powders;

public:
	Parser() {}
	Parser(const std::string &path, char mode);
	~Parser();

	Barrel& readBarrel();
	const Powders& readXMLPowders(const std::string &path);
	void writeBarrel(const Barrel &barr);
};

#endif

#ifndef PARSER_H
#define PARSER_H

#include <fstream>
#include <string>
#include "pugixml.hpp"

#include "Types.h"

class Parser
{
private:
	std::fstream file;
	Barrel barr;
	Powders powders;

public:
	Parser() {}
	Parser(const std::string &path, char mode);
	~Parser() {
		file.close();
	}

	bool open(const std::string &path, char mode);
	bool close();

	double readNext();
	std::string readStr();
	Barrel& readBarrel();
	const Powders& readXMLPowders(const std::string &path);

	void write(const std::string &txt);
	void write(const std::string &txt, double x, char split = '\n');
	void write(double x);
	void write(char split, double x);
	void write(double x, char split);
	void writeBarrel(const Barrel &barr);

	bool isOpen() const;
	bool isEnd() const;

	static void createFile(const std::string &path);
};

#endif

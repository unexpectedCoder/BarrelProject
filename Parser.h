#ifndef PARSER_H
#define PARSER_H

#include <fstream>
#include <string>

#include "Types.h"

class Parser
{
private:
	std::fstream file;
	Barrel barr;

public:
	Parser(const std::string &path, char mode);
	~Parser() {
		file.close();
	}

	bool open(const std::string &path, char mode);
	bool close();

	double readNext();
	std::string readStr();
	Barrel& readBarrel();
	void write(const std::string &txt);
	void write(double x, char split);
	void writeBarrel(const Barrel &barr);

	bool isEnd() const;

	static void createFile(const std::string &path);
};

#endif

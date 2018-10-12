#ifndef PARSER_H
#define PARSER_H

#include <fstream>
#include <string>

class Parser
{
private:
	std::fstream file;

public:
	Parser(const std::string &path) {
		file.open(path.c_str());
	}
	~Parser() {
		file.close();
	}

	bool open(const std::string &path);
	bool close();

	double readNext();
	std::string readStr();
	void write(const std::string &txt);
	void write(double x, char split);

	bool isEnd() const;

	static void createFile(const std::string &path);
};

#endif

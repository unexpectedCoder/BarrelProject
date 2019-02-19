#ifndef ANALOGSSOLVER_H
#define ANALOGSSOLVER_H

#include "AbstractSolver.h"
#include "Parser.h"

#define ANALOGS_PATH "files/analogs.txt"
#define ANALOGS_RES_PATH "files/analogs_res.txt"


class AnalogsSolver : public AbstractSolver
{
private:
	Analogs analogs;
	std::string status;

public:
	AnalogsSolver() : status("successfully") {}
	~AnalogsSolver() { analogs.~vector(); }

	void printIntro();
	void printOutro();
	void solve();

private:
	void fillAnalogsData();
	void writeAnalogsFile();
	void writeSolveFile(Parser &par);
};

#endif

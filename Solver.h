#ifndef SOLVER_H
#define SOLVER_H

#include <string>

#include "Types.h"

class Solver
{
private:
	int num_of_as;
	Analogs analogs;

	Chuev linInterp(double x, const std::string &path_chuev = "chuev.txt");

public:
	Solver() : num_of_as(0) {}
	~Solver() {}

	void fillAnalogs(const std::string &path = "analogs.txt");
	Analogs& calcAnalogs(const std::string &path = "analogs.txt");

	double countCE15(double cq, double ce, double eta);
};

#endif

#ifndef SOLVER_H
#define SOLVER_H

#include <string>
#include <vector>

#include "Types.h"

class Solver
{
private:
	int num_of_as;
	Analogs analogs;
	Barrels barrs;
	std::vector<double> Delta;
	std::vector<double> eta_K;

	Chuev linInterp(double x, const std::string &path = "chuev.txt");
	void makeTableTxt(const Barrel &barr, double pm_nround, const std::string &path);
	void fillData(const std::string &txt, std::vector<double> &data);

public:
	Solver() : num_of_as(0) {}
	~Solver() {}

	void makeTest();
	void fillAnalogs(const std::string &path = "analogs.txt");
	Analogs& calcAnalogs(const std::string &path = "analogs.txt");
	void calcBarrelPressure(const std::string &path = "barrel_pressure.txt");
	Barrels& solveInvProblem();

	static double calcCE15(double cq, double ce, double eta);
};

#endif

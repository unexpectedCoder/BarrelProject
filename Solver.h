#ifndef SOLVER_H
#define SOLVER_H

#include <string>
#include <vector>

#include "Types.h"

#define NEW_BARR_K 1.2

#define ROWS 18
#define COLOUMS 11

#define ANALOGS_PATH "files/analogs.txt"
#define BARR_TABLE_PATH "files/barrel_table.txt"
#define TEST_PATH "files/test.txt"
#define P_CE15_PATH "files/p_CE15.txt"
#define BARR_SRC_PATH "files/barrel_src.txt"
#define B_DELTA_PATH "files/B(Delta).txt"
#define CHUEV_TABLE_33_PATH "files/chuev_table_3.3.txt"
#define Z_SLUH_PATH "files/ZSluh.txt"

class Solver
{
private:
	int num_of_as;
	Analogs analogs;
	Barrels barrs;
	std::vector<double> Delta;
	std::vector<double> eta_K;

	static int fcounter;

	void makeTableTxt(const Barrel &barr, double pm_nround, const std::string &path);
	void fillData(const std::string &txt, std::vector<double> &data);
	void writeBarrelToFile(const Barrel &barr);
	double belinearInterp(double hi_n, double Lambda_D, double **table);

public:
	Solver() : num_of_as(0) {}
	~Solver() {}

	void makeTest();
	void fillAnalogs(const std::string &path = ANALOGS_PATH);
	Analogs& calcAnalogs(const std::string &path = ANALOGS_PATH);
	void calcBarrelPressure(const std::string &path = BARR_TABLE_PATH);
	Barrels& solveInvProblem();

	static double calcCE15(double cq, double ce, double eta);
};

#endif

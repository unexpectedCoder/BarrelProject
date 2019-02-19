#ifndef ANALITICSOLVER_H
#define ANALITICSOLVER_H

#include "AbstractSolver.h"
#include "Parser.h"
#include "Types.h"

#define NEW_BARR_K 1.2	// The coefficient taking into account
												// the use of new construction materials
#define BARR_TABLE_PATH "files/barrel_table.txt"
#define BARR_SRC_PATH "files/barrel_src.txt"
#define BARR_LOG_PATH "files/barrel_log.txt"
#define CHUEV_TABLE_33_PATH "chuev_table_3.3.txt"
#define Z_SLUH_PATH "results/ZSluh.txt"
#define B_DELTA_PATH "files/B(Delta).txt"
#define RESULTS_PATH "results/results.txt"

#define LD_PATH "results/LD.txt"
#define IK_PATH "results/Ik.txt"
#define W0_PATH "results/W0.txt"
#define W_PATH "results/W.txt"
#define W_Q_PATH "results/w_q.txt"


class AnaliticSolver : public AbstractSolver
{
private:
	Barrels barrs;
	std::string status;
	std::vector<double> Delta;
	std::vector<double> eta_K;

public:
	AnaliticSolver() : status("successfully") {}
	~AnaliticSolver();

	void printIntro();
	void printOutro();
	void calcMaxPressure();
	void solve();

private:
	void makeTableTxt(const Barrel &barr, double pm_nround, const std::string &path);
	void fillBarrelData(const std::string &txt, std::vector<double> &data);
	void writeBarrelsToFile();
	void writeFileFromMaxPress(const Barrel &b, double pm_nround);
	void createNecessaryFiles();
};

#endif

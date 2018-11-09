#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include "Parser.h"
#include "Types.h"

#define NEW_BARR_K 1.2

#define ROWS 18
#define COLOUMS 11

#define ANALOGS_PATH "files/analogs.txt"
#define ANALOGS_RES_PATH "files/analogs_res.txt"
#define BARR_TABLE_PATH "files/barrel_table.txt"
#define TEST_PATH "files/test.txt"
#define P_CE15_PATH "files/p_CE15.txt"
#define BARR_SRC_PATH "files/barrel_src.txt"
#define B_DELTA_PATH "files/B(Delta).txt"
#define CHUEV_TABLE_33_PATH "chuev_table_3.3.txt"
#define Z_SLUH_PATH "results/ZSluh.txt"
#define RESULTS_PATH "results/results.txt"
#define BARR_LOG_PATH "files/barrel_log.txt"

#define POWDERS_PATH "powders.xml"
#define PM_PATH "direct_res/pm"

// -------- БАЗОВЫЙ КЛАСС -------- //
class Solver
{
public:
	virtual void printInfo() = 0;
	virtual void solve() = 0;
	virtual void printResults() = 0;

	virtual double calcCE15(double cq, double ce, double eta) {
		return 0.5 * 15.0 / cq * (ce - 3.0 * eta * cq + sqrt(pow(ce - 3.0 * eta * cq, 2.0) + 4.0 / 5.0 * ce * eta * pow(cq, 2.0)));
	}
	virtual double bilinearInterp(double hi_n, double Lambda_D, double **table);

	virtual ~Solver() {}
};
// ------------------------------- //

// Тестовый решатель
class TestSolver : public Solver
{
private:
	Barrel *barr;
	std::string status;

public:
	TestSolver() : status("successfully") {
		barr = new Barrel(StartData(0.122, 21.76, 690, 3e7, 1.05, 1.04));
	}

	void printInfo() {
		std::cout << "\n\t<Тестовый решатель (TestSolver)>\n";
	}
	void solve();
	void printResults() {
		std::cout << "\tСТАТУС ВЫПОЛНЕНИЯ: " << status << ".\n";
		std::cout << "\tРЕЗУЛЬТАТЫ: см. файл " << TEST_PATH << ".\n\n";
	}

	~TestSolver() {
		delete barr;
	}
};

// Решатель для аналогов
class AnalogsSolver : public Solver
{
private:
	Analogs analogs;
	std::string status;

	void fillAnalogsData();

public:
	AnalogsSolver() : status("successfully") {}

	void printInfo() {
		std::cout << "\n\t<Решатель для аналогов (AnalogsSolver)>\n";
	}
	void solve();
	void printResults() {
		std::cout << "\tСТАТУС ВЫПОЛНЕНИЯ: " << status << ".\n";
		std::cout << "\tРЕЗУЛЬТАТЫ: см. файл " << ANALOGS_PATH << ".\n\n";
	}

	~AnalogsSolver() {
		analogs.~vector();
	}
};

// Решатель для аналитического решения ОЗВБ
class AnaliticSolver : public Solver
{
private:
	Barrels barrs;
	std::string status;
	std::vector<double> Delta;
	std::vector<double> eta_K;
	static int file_count;

	void makeTableTxt(const Barrel &barr, double pm_nround, const std::string &path);
	void fillBarrelData(const std::string &txt, std::vector<double> &data);
	void writeBarrelsToFile();

public:
	AnaliticSolver() : status("successfully") {}
	~AnaliticSolver() {
		barrs.~vector();
		Delta.~vector();
		eta_K.~vector();
	}

	void printInfo() {
		std::cout << "\n\t<Аналитический решатель (AnaliticSolver)>\n";
	}
	void calcMaxPressure();
	void solve();
	void printResults() {
		std::cout << "\tСТАТУС ВЫПОЛНЕНИЯ: " << status << ".\n";
		std::cout << "\tРЕЗУЛЬТАТЫ: см. в папке results.";
	}
};

class DirectSolver : public Solver
{
private:
	std::string status;
	Powder pwd;
	Powders pwds;
	double *Delta, *w_q;
	double d, q, Vd, ns, S, K;
	double p0, pm;
	int size_d, size_wq;
	int key_V, key_S, key_Z;

	void fillData();
	int choosePowder();
	void searchPmax(double dt, const Params &params, Results &rs);
	void continueCalc(double dt, const Result &start, Results &rs);
	Results::const_iterator minDeltaPm(const Results &rs);
	void rksolve(double dt, const Params &params, Result &res);
	void writeResultsToFile(const std::string &path, const Results &rs);
	void findPmLine(double dt);

	void setPathPm(std::string &path, unsigned num = 0);
	void createFilePm(const std::string &path);
	void writeFilePm(const std::string &path, const Result &res);

	// Система ОДУ
	double dz(double p) {
		return p / pwd.Ik * key_Z;
	}
	double dpsi(double z, double p) {
		return (pwd.kappa1 * (1.0 + 2.0 * pwd.lambda1 * z) * key_S +
			pwd.kappa2 * (1.0 + 2.0 * pwd.lambda2 * (z - 1.0)) * (1 - key_S)) * dz(p);
	}
	double dL(double V) {
		return V;
	}
	double dV(double p, double w_q, double fi) {
		return p * S / (fi * q) * key_V;
	}
	double dW(double omega, double z, double p, double V) {
		return (1.0 - pwd.alpha * pwd.delta) / pwd.delta * omega * dpsi(z, p) + S * V;
	}
	double dp(double W, double omega, double _Delta, double p, double z, double L, double V, double F0) {
		return 1.0 / W * (pwd.f * omega * dpsi(z, p) -
			(pwd.k - 1.0) * Consts::sigma_T * Consts::nu_T * p * (F0 + pi * d * L) / pwd.R() -
			pwd.k * p * dW(omega, z, p, V));
	}
	//

public:
	DirectSolver(double _pm = 410e6, double _d = 0.1, double _q = 4.35, double _Vd = 1600, double _K = 1.03, double _p0 = 1e7, double _ns = 1) :
		status("successfully"), pm(_pm), d(_d), q(_q), Vd(_Vd), K(_K), p0(_p0), ns(_ns), key_V(0), key_S(1), key_Z(1)
	{
		pwds = Parser().readXMLPowders(POWDERS_PATH);
		S = 0.25 * pi * d * d * ns;
	}
	~DirectSolver() {
		delete[] Delta;
		delete[] w_q;
		pwds.~vector();
	}

	void printInfo() {
		std::cout << "\n\t<Решатель прямой задачи (DirectSolver)>\n";
	}

	void solve();
	int showPowders();

	void printResults() {
		std::cout << "\tСТАТУС ВЫПОЛНЕНИЯ: " << status << ".\n";
		std::cout << "\tРЕЗУЛЬТАТЫ: см. в папке direct_res.";
	}
};

#endif

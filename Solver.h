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

#define DIRSOL_TEST_PATH "direct_res/test.txt"
#define POWDERS_PATH "powders.xml"
#define PM_PATH "direct_res/pm"
#define I_DIAG_PATH "direct_res/IDiag"
#define LM_W_PATH "direct_res/lm(W)"
#define LM_DELTA_PATH "direct_res/lm(Delta)"
#define Z_PATH "direct_res/Z"
#define Z_MAX_PATH "direct_res/Zmax"
#define SOL_ONCE_PATH "direct_res/FullSol"

class Error
{
	std::string mess;

public:
	Error() : mess("No errors.") {}
	~Error() {
		mess.~basic_string();
	}

	const std::string& sendMess(const std::string &_mess) {
		return mess = _mess;
	}
};

// -------- БАЗОВЫЙ КЛАСС -------- //
class Solver
{
public:
	virtual void printIntro() = 0;
	virtual void solve() = 0;
	virtual void printOutro() = 0;

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

	void printIntro() {
		std::cout << "\n\t<Тестовый решатель (TestSolver)>\n";
	}
	void solve();
	void printOutro() {
		std::cout << "СТАТУС ВЫПОЛНЕНИЯ: " << status << ".\n";
		if (status != "failed")
			std::cout << "РЕЗУЛЬТАТЫ: см. файл " << TEST_PATH << ".\n\n";
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

	void printIntro() {
		std::cout << "\n\t<Решатель для аналогов (AnalogsSolver)>\n";
	}
	void solve();
	void printOutro() {
		std::cout << "СТАТУС ВЫПОЛНЕНИЯ: " << status << ".\n";
		if (status != "failed")
			std::cout << "РЕЗУЛЬТАТЫ: см. файл " << ANALOGS_PATH << ".\n\n";
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

	void printIntro() {
		std::cout << "\n\t<Аналитический решатель (AnaliticSolver)>\n";
	}
	void calcMaxPressure();
	void solve();
	void printOutro() {
		std::cout << "СТАТУС ВЫПОЛНЕНИЯ: " << status << ".\n";
		if (status != "failed")
			std::cout << "РЕЗУЛЬТАТЫ: см. в папке results.\n\n";
	}
};

class DirectSolver : public Solver
{
private:
	unsigned size_d, size_wq;
	Error err;
	std::string status;
	Powder pwd;
	Powders pwds;
	double *Delta, *w_q;
	double d, q, Vd, ns, S, K;
	double p0, pm;
	double l_d_max;
	Results results;
	double Ld;

	void fillDelta();
	void fill_wq();
	double setTimeStep();
	int choosePowder();
	void searchPmaxConds(double dt, double delta, Result &res);
	void calcToPmax(double dt, Result &res);
	void continueCalc(double dt, Result &res);
	void calcToPmax(double dt, Result &res, Results &rs);
	void continueCalcVd(double dt, Result &res, Results &rs);		// Обрезка по дул. скорости
	void continueCalcLd(double dt, Result &res, Results &rs);		// ... по длине ведущей части
	
	void set_l_d_max();
	void calcPmLine(double dt, unsigned indx);
	void calcIndicatDiag(double dt, unsigned indx);
	void rksolve(double dt, Result &res);

	void setPath(const std::string &base_path, std::string &res_path, int num = 0);
	void createFile(const std::string &path, const std::string &head = "", bool w_name = true);
	void writeFilePm(const std::string &path, const Result &res);
	void writeResultsToFile(const std::string &path, const Results &rs);
	void writeFileDiag(const std::string &path, const Result &res);

	double funcW0(double w);
	void writeLmFile(unsigned indx = 0);

	void fillCriterionData(CriterionParams &cr);
	void calcCriterionCoeffs(const CResults &rs, CriterionParams &cp);
	Criterions::iterator maxCriterion(const Criterions::iterator &start,
																		const Criterions::iterator &end);
	void fillCriterions(const CResults &rs, const CriterionParams &cp, Criterions &crs);
	void readResults(const std::string &path, CResults &rs);
	void getMaxCriterion(const std::string &path, Criterion &max_cr);
	void readCriterionFile(const std::string &path, Criterions &crs);
	void writeCriterionsFile(const std::string &path, const Criterions &crs);
	void writeMaxCriterionFile(const std::string &path, const Criterions &crs);

	void solveStandartTemp(double dt);
	void writeResultFile(const std::string &path, const Results &rs);

	// Система ОДУ
	double dz(double z, double p) {
		return p / pwd.Ik * ksi_end(z);
	}
	double dpsi(double z, double p) {
		return (pwd.kappa1 * (1.0 + 2.0 * pwd.lambda1 * z) * ksi_s(z) +
			pwd.kappa2 * (1.0 + 2.0 * pwd.lambda2 * (z - 1.0)) * (1 - ksi_s(z))) * dz(z, p);
	}
	double dL(double V) {
		return V;
	}
	double dV(double fi, double V, double p) {
		return p * S / (fi * q) * ksi_v(p, V);
	}
	double dW(double omega, double z, double p, double V) {
		return (1.0 - pwd.alpha * pwd.delta) / pwd.delta * omega * dpsi(z, p) + S * V;
	}
	double dp(double omega, double F0, double W, double p, double z, double L, double V) {
		return 1.0 / W * (pwd.f * omega * dpsi(z, p) -
			(pwd.k - 1.0) * Consts::sigma_T * Consts::nu_T * p * (F0 + pi * d * L) / pwd.R() -
			pwd.k * p * dW(omega, z, p, V));
	}
	//
	// Функции-переключатели
	int ksi_end(double z) {
		if (z > pwd.zk)
			return 0;
		return 1;
	}
	int ksi_v(double p, double V) {
		if (p > p0 || V > 0)
			return 1;
		return 0;
	}
	int ksi_s(double z) {
		if (z > 1)
			return 0;
		return 1;
	}
	//

public:
	DirectSolver(double _pm = 410e6, double _d = 0.1, double _q = 4.35, double _Vd = 1600, double _K = 1.03, double _p0 = 1e7, double _ns = 1) :
		status("successfully"), pm(_pm), d(_d), q(_q), Vd(_Vd), K(_K), p0(_p0), ns(_ns)
	{
		pwds = Parser().readXMLPowders(POWDERS_PATH);
		S = 0.25 * pi * pow(d, 2.0) * ns;
		Ld = 0;
	}
	~DirectSolver() {
		delete[] Delta;
		delete[] w_q;
		pwds.~vector();
	}

	void printIntro() {
		std::cout << "\n\t<Решатель прямой задачи (DirectSolver)>\n";
	}

	void makeTest(const TestParams &tp);
	void solve();
	void calcCriterions();
	void solveOnce(double dt, int start_temp_C);
	int showPowders();

	double getBarLen() const {
		return Ld;
	}

	void printOutro() {
		std::cout << "СТАТУС ВЫПОЛНЕНИЯ: " << status << ".\n";
		if (status != "failed")
			std::cout << "РЕЗУЛЬТАТЫ: см. в папке direct_res.\n\n";
	}
};

#endif

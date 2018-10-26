#ifndef SOLVER_H
#define SOLVER_H

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

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
#define CHUEV_TABLE_33_PATH "files/chuev_table_3.3.txt"
#define Z_SLUH_PATH "files/ZSluh.txt"
#define RESULTS_PATH "files/results.txt"

// -------- БАЗОВЫЙ КЛАСС -------- //
class Solver
{
public:
	virtual void printInfo() = 0;
	virtual int solve() = 0;
	virtual void printResults() = 0;

	virtual double calcCE15(double cq, double ce, double eta) {
		return 0.5 * 15.0 / cq * (ce - 3.0 * eta * cq + sqrt(pow(ce - 3.0 * eta * cq, 2.0) + 4.0 / 5.0 * ce * eta * pow(cq, 2.0)));
	}
	virtual double belinearInterp(double hi_n, double Lambda_D, double **table);

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
	int solve();
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
	int solve();
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
	~AnaliticSolver() {}

	void printInfo() {
		std::cout << "\n\t<Аналитический решатель (AnaliticSolver)>\n";
	}
	void calcMaxPressure();
	int solve();
	void printResults() {
		std::cout << "\tСТАТУС ВЫПОЛНЕНИЯ: " << status << ".\n";
		std::cout << "\tРЕЗУЛЬТАТЫ: см. файл barrels_xxxx.txt, " <<
			P_CE15_PATH << " и " << Z_SLUH_PATH << ".\n\n";
	}
};

#endif

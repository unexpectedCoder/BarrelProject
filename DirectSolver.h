#ifndef DIRECTSOLVER_H
#define DIRECTSOLVER_H

#include "AbstractSolver.h"
#include "Parser.h"
#include "Types.h"

#define DIRSOL_TEST_PATH "direct_res/test.txt"
#define POWDERS_PATH "powders.xml"
#define PM_PATH "direct_res/pm"
#define I_DIAG_PATH "direct_res/IDiag"
#define LM_W_PATH "direct_res/lm(W)"
#define LM_DELTA_PATH "direct_res/lm(Delta)"
#define Z_PATH "direct_res/Z"
#define Z_MAX_PATH "direct_res/Zmax"
#define SOL_ONCE_PATH "direct_res/FullSol"


class DirectSolver : public AbstractSolver
{
private:
	unsigned size_d, size_wq;
	std::string status;
	Powder pwd;
	Powders pwds;
	double *Delta, *w_q;
	double d, q, Vd, ns, S, K;
	double p0, pm;
	double l_d_max;
	Results results;
	double Ld;

public:
	DirectSolver(double _pm = 410e6, double _d = 0.1,
							 double _q = 4.35, double _Vd = 1600, double _K = 1.03,
							 double _p0 = 1e7, double _ns = 1);
	~DirectSolver();

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

private:
	void fillDelta();
	void fill_wq();
	double setTimeStep();
	int choosePowder();
	void searchPmaxConds(double dt, double delta, Result &res);
	void calcToPmax(double dt, Result &res);
	void continueCalc(double dt, Result &res);
	void calcToPmax(double dt, Result &res, Results &rs);
	void continueCalcVd(double dt, Result &res, Results &rs);		// Окончание по дул. скорости
	void continueCalcLd(double dt, Result &res, Results &rs);		// и по длине ведущей части

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

	// Switch-functions
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
};

#endif

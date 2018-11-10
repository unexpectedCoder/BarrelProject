#ifndef TYPES_H
#define TYPES_H

#include <math.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

const double eps = 1e-4;
const double pi = 3.141592654;

namespace Consts {
	// Общие константы для решения ОЗВБ
	const double g = 9.80665;
	const double fi1 = 1.02;
	const double k = 1.25;
	const double Nkr = 1.23;
	const double f = 1e6;
	const double alpha_k = 1e-3;
	const double delta = 1600;
	const double kapa = 1;
	const double lambda = 0;
	const double p_flash = 1e6;
	// Для критерия Слухоцкого
	const double C = 1e6;
	// Для прямого перебора
	const double sigma_T = 360;
	const double nu_T = 0.7;
}

struct Analog
{
	std::string name;
	double d;
	double q;
	double vd;

	double Cq, CE;
	double CE15;
	double eta_omega;
	double pm;
	double l_d;

	Analog(const std::string &_name, double _d, double _q, double _vd, double _l_d, double _pm) :
		name(_name), d(_d), q(_q), vd(_vd), l_d(_l_d), pm(_pm)
	{
		Cq = q / pow(d * 10, 3.0);
		CE = q * pow(vd, 2.0) / (2e3 * Consts::g * pow(d * 10.0, 3.0));
	}
};
typedef std::vector<Analog> Analogs;

struct Chuev
{
	double pm_kr;
	double eta_omega;
	double hi;

	Chuev(double CE)
	{
		double CE1, CE2;
		double p1, p2;
		double eta1, eta2;
		double hi1, hi2;

		std::fstream file("chuev.txt", std::ios_base::in);
		while (!file.eof())
		{
			file >> CE1;
			file >> p1;
			file >> eta1;
			file >> hi1;

			if (CE - CE1 < 100)
			{
				file >> CE2;
				file >> p2;
				file >> eta2;
				file >> hi2;

				break;
			}
		}

		pm_kr = p1 + (CE - CE1) / (CE2 - CE1) * (p2 - p1);
		eta_omega = eta1 + (CE - CE1) / (CE2 - CE1) * (eta2 - eta1);
		hi = hi1 + (CE - CE1) / (CE2 - CE1) * (hi2 - hi1);
	}
};

struct StartData
{
	double d;
	double q;
	double vd;

	double p0;
	double K;
	double ns;

	StartData(double _d = 0.1, double _q = 4.35, double _vd = 1600, double _p0 = 1e7, double _K = 1.03, double _ns = 1.0) :
		d(_d), q(_q), vd(_vd), p0(_p0), K(_K), ns(_ns) {}
};

struct Barrel
{
	// Исходные данные
	double d, q, vd;
	// Рассчитываемые величины
	double Cq;
	double CE;
	double CE15;
	double eta_omega;
	double omega_q;
	double fi;
	double pm_kr;
	double p0, pm;
	double hi, hi1;
	double ns;
	double p_mid;
	// Для ОЗВБ
	double Delta;
	double psi0;
	double sigma0;
	double z0;
	double K, K1;
	double B, B1;
	double gamma1;
	double alpha1;
	double beta1, beta2, beta_m;
	double beta_m_1, beta_m_2;
	double fi_m;
	double pi_m;
	double beta_k;
	double beta_k1, beta_k2;
	double fi_k;
	double Lambda_K;
	double r_K;
	// Переменные результатов
	double eta_K;
	double Lambda_D;
	double r_D, r_Dmin;
	double omega;
	double W0, W;
	double L0;
	double LD;
	double Ik;
	double Z_Sluh;

	Barrel(const StartData &data = StartData()) :
		d(data.d), q(data.q), vd(data.vd), p0(data.p0), K(data.K), ns(data.ns)
	{
		Cq = q / pow(d * 10.0, 3.0);
		CE = q * pow(vd, 2.0) / (2e3 * Consts::g * pow(d * 10.0, 3.0));
		r_Dmin = 1.0 / 6.0 * pow(vd, 2.0) / (Consts::f / (Consts::k - 1.0));
	}
	Barrel(const Barrel &barr, const StartData &data) :
		p0(data.p0), K(data.K)
	{
		d = barr.d;
		q = barr.q;
		vd = barr.vd;

		Cq = barr.Cq;
		CE = barr.CE;
		CE15 = barr.CE15;
		eta_omega = barr.eta_omega;
		omega = q * barr.omega;
		omega_q = barr.omega_q;
		pm_kr = barr.pm_kr;
		pm = barr.pm;
		hi = barr.hi;
		ns = barr.ns;

		r_Dmin = 1.0 / 6.0 * pow(vd, 2.0) / (Consts::f / (Consts::k - 1.0));
	}

	void calcHi1();
	void calcForTest(double _B);
	void calcB(double a, double b);
	void halfSearchB(double a, double b);
	void calcLambdaK();
	void calcForEtaK(double eta_K);
	double calcZSluh();
	double calcZSluh(double v1_vd);
};
typedef std::vector<Barrel> Barrels;

struct Matrix
{
	double **data;
	int x, y;

	Matrix() : x(0), y(0) {}
	Matrix(int _x, int _y) {
		if (_x > 0 && _y > 0)
		{
			x = _x;
			y = _y;

			data = new double*[x];
			for (int i = 0; i < x; i++)
				data[i] = new double[y];
		}
		else throw "Error (in struct \'Matrix\'): invalid matrix size!";
	}
	Matrix(double **mtrx, int _x, int _y) {
		if (_x > 0 && _y > 0)
		{
			x = _x;
			y = _y;

			data = new double*[x];
			for (int i = 0; i < x; i++)
			{
				data[i] = new double[y];
				for (int j = 0; j < y; j++)
					data[i][j] = mtrx[i][j];
			}
		}
		else throw "Error (in struct \'Matrix\'): invalid matrix size!";
	}
	Matrix(const Matrix &mtrx) {
		x = mtrx.x;
		y = mtrx.y;

		data = new double*[x];
		for (int i = 0; i < x; i++)
		{
			data[i] = new double[y];
			for (int j = 0; j < y; j++)
				data[i][j] = mtrx.data[i][j];
		}
	}

	~Matrix() {
		for (int i = 0; i < x; i++)
			delete[] data[i];
		delete[] data;
	}

	Matrix& T();
	Matrix& zeros();

	friend std::ostream& operator<<(std::ostream &os, const Matrix &mtrx);
};
inline std::ostream& operator<<(std::ostream &os, const Matrix &mtrx)
{
	for (int i = 0; i < mtrx.x; i++)
	{
		os << "[\t";
		for (int j = 0; j < mtrx.y; j++)
			os << mtrx.data[i][j] << '\t';
		os << "]\n";
	}

	return os;
}

struct Powder
{
	std::string name;
	double
		f,
		k,
		alpha,
		T,
		delta,
		Ik,
		zk,
		kappa1,
		lambda1,
		kappa2,
		lambda2,
		kappa_f,
		k_f;

	double R() {
		return f / T;
	}
	friend std::ostream& operator<<(std::ostream &os, const Powder &p);
};
typedef std::vector<Powder> Powders;

inline std::ostream& operator<<(std::ostream &os, const Powder &p)
{
	os << "\tПорох " << p.name << ":\n";
	os << "\t - f = " << p.f * 1e-6 << " (МДж/кг);\n";
	os << "\t - k = " << p.k << ";\n";
	os << "\t - alpha = " << p.alpha * 1e3 << " (дм^3/кг);\n";
	os << "\t - delta = " << p.delta * 1e-3 << " (кг/дм^3);\n";
	os << "\t - T = " << p.T << " (К);\n";
	os << "\t - Ik = " << p.Ik * 1e-6 << " (МПа*с);\n";
	os << "\t - zk = " << p.zk << ";\n";
	os << "\t - kappa1 = " << p.kappa1 << ";\n";
	os << "\t - lambda1 = " << p.lambda1 << ";\n";
	os << "\t - kappa2 = " << p.kappa2 << ";\n";
	os << "\t - lambda2 = " << p.lambda2 << ";\n";
	os << "\t - kappa_f = " << p.kappa_f << ";\n";
	os << "\t - k_f = " << p.k_f << ";\n";

	return os;
}

struct Result
{
	double t;
	double p, W, W0;
	double V, L;
	double psi, z;
	double Delta, w_q;

	Result() {
		t = 0.0;
		L = 0.0;
		p = Consts::p_flash;
		psi = 0.0;
		V = 0.0;
		z = 0.0;
	}

	Result& operator=(const Result &other);
	friend std::ostream& operator<<(std::ostream &os, const Result &res);
};
typedef std::vector<Result> Results;

inline std::ostream& operator<<(std::ostream &os, const Result &r)
{
	os << r.t << ' ' << r.Delta << ' ' << r.w_q << ' ' << r.p << ' ' <<
		r.V << ' ' << r.L << ' ' << ' ' << r.W0 << ' ' << r.W << ' ' << r.z << '\n';

	return os;
}

struct Params
{
	double
		Delta,
		w_q,
		W0,
		W,
		fi,
		F0;
};

#endif

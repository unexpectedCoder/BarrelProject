#ifndef TYPES_H
#define TYPES_H

#include <math.h>
#include <vector>
#include <string>
#include <iostream>

const double eps = 1e-4;

namespace Consts {
	// Общие константы для решения ОЗВБ
	const double pi = 3.141592654;
	const double g = 9.80655;
	const double fi1 = 1.02;
	const double K = 1.03;
	const double k = 1.25;
	const double Nkr = 1.23;
	const double f = 1e6;
	const double alpha_k = 1e-3;
	const double delta = 1600;
	const double kapa = 1;
	const double lambda = 0;
	const double p_flash = 1e6;
}

struct Analog
{
	std::string name;
	double d;
	double q;
	double vd;

	double CE15;
	double pm;
};
typedef std::vector<Analog> Analogs;

///////////////////////////////////////////////////////////////////////////////
struct Chuev
{
	double CE;
	double pm_kr;
	double eta_omega;

	Chuev(double _CE = 0, double _pm_kr = 0, double _eta_omega = 0) :
		CE(_CE), pm_kr(_pm_kr), eta_omega(_eta_omega) {}

	friend std::ostream& operator<<(std::ostream &os, const Chuev &chuev);
};

inline std::ostream& operator<<(std::ostream &os, const Chuev &chuev)
{
	os << "CE = " << chuev.CE << '\t'
		<< "p = " << chuev.pm_kr << '\t' << "eta_omega = " << chuev.eta_omega;
	return os;
}
///////////////////////////////////////////////////////////////////////////////

struct Barrel
{
	// Исходные данные
	double d;
	double q;
	double vd;
	// Рассчитываемые величины
	double Cq;
	double CE;
	double CE15;
	double eta_omega;
	double omega_q;
	double fi;
	double pm_kr;
	double p0, pm;
	double hi1;
	double ns;
	double p_mid;
	// Для ОЗВБ
	double Delta;
	double psi0;
	double sigma0;
	double z0;
	double K1;
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

	double eta_K;
	double Lambda_D;
	double r_D;
	double omega;
	double W0;
	double L0;
	double LD;
	double Ik;

	Barrel(double _d = 0.1, double _q = 4.35, double _vd = 1600, double _p0 = 1e7, double _ns = 1) :
		d(_d), q(_q), vd(_vd), p0(_p0), ns(_ns)
	{
		Cq = q / pow(d * 10, 3.0);
		CE = q * pow(vd, 2.0) / (2e3 * Consts::g * pow(d, 3.0)) * 1e-3;
	}

	void calcForTest(double _B)
	{
		B = _B;
		Delta = 500;
		p_mid = 0.5 * pm;
		eta_K = 0.5;

		hi1 = 1.0 - (Consts::k - 1.0) * p_mid / Consts::f * (1.0 - Consts::alpha_k * Consts::delta) / Consts::delta;

		psi0 = (1.0 / Delta - 1.0 / Consts::delta) / (Consts::f / p0 - (1.0 - Consts::alpha_k * Consts::delta) / Consts::delta);
		sigma0 = sqrt(1.0 + 4.0 * Consts::lambda / Consts::kapa * psi0);
		z0 = 2.0 * psi0 / (Consts::kapa * (1.0 + sigma0));
		K1 = Consts::kapa * sigma0;

		B1 = 0.5 * (Consts::k - 1) * B - Consts::kapa * Consts::lambda * hi1;
		gamma1 = psi0 * B1 / pow(K1, 2.0);
		alpha1 = 2 * Consts::kapa * Consts::lambda / B1;
		beta1 = 0.5 * hi1 + sqrt(gamma1 + 0.25 * pow(hi1, 2.0));
		beta2 = 0.5 * hi1 - sqrt(gamma1 + 0.25 * pow(hi1, 2.0));

		beta_m = (Consts::k * hi1 - 1) / (2.0 * Consts::k + alpha1);
		beta_m_1 = beta_m / beta1;
		beta_m_2 = beta_m / beta2;

		fi_m = pow(1.0 - beta_m_2, (1.0 + alpha1 * beta2) / (beta1 - beta2) - 1.0) /
			pow(1.0 - beta_m_1, (1.0 + alpha1 * beta1) / (beta1 - beta2) + 1.0);

		pi_m = (1 - beta_m_1) * (1 - beta_m_2) * pow(fi_m, -1.0 / (Consts::k - 1.0));

		calcLambdaK();
		calcOtherThings();
	}

	void calcB(double a, double b)
	{
		psi0 = (1 / Delta - 1 / Consts::delta) / (Consts::f / p0 - (1 - Consts::alpha_k * Consts::delta) / Consts::delta);
		sigma0 = sqrt(1 + 4.0 * Consts::lambda / Consts::kapa * psi0);
		z0 = 2 * psi0 / (Consts::kapa * (1 + sigma0));
		K1 = Consts::kapa * sigma0;

		halfSearchB(a, b);
	}

	void halfSearchB(double a, double b)
	{
		double pi_m_star = pm / p0;

		while (true)
		{
			B = 0.5 * (a + b);

			B1 = 0.5 * (Consts::k - 1) * B - Consts::kapa * Consts::lambda * hi1;
			gamma1 = psi0 * B1 / pow(K1, 2.0);
			alpha1 = 2 * Consts::kapa * Consts::lambda / B1;
			beta1 = 0.5 * hi1 + sqrt(gamma1 + 0.25 * pow(hi1, 2.0));
			beta2 = 0.5 * hi1 - sqrt(gamma1 + 0.25 * pow(hi1, 2.0));

			beta_m = (Consts::k * hi1 - 1) / (2.0 * Consts::k + alpha1);
			beta_m_1 = beta_m / beta1;
			beta_m_2 = beta_m / beta2;

			fi_m = pow(1.0 - beta_m_2, (1.0 + alpha1 * beta2) / (beta1 - beta2) - 1.0) /
				pow(1.0 - beta_m_1, (1.0 + alpha1 * beta1) / (beta1 - beta2) + 1.0);

			pi_m = (1 - beta_m_1) * (1 - beta_m_2) * pow(fi_m, -1.0 / (Consts::k - 1.0));

			if (fabs(pi_m - pi_m_star) < eps) break;

			if (pi_m > pi_m_star) a = B;
			else b = B;
		}
	}

	void calcLambdaK()
	{
		beta_k = B1 / K1 * (1.0 - z0);
		beta_k1 = beta_k / beta1;
		beta_k2 = beta_k / beta2;
		fi_k = pow(1.0 - beta_k2, (1 + alpha1 * beta2) / (beta1 - beta2) - 1.0) /
			pow(1.0 - beta_k1, (1 + alpha1 * beta1) / (beta1 - beta2) + 1.0);

		Lambda_K = Consts::f * Delta * psi0 / p0 * (pow(fi_k, 1.0 / (Consts::k - 1.0)) - 1.0) -
			(1.0 - Consts::alpha_k * Consts::delta) / Consts::delta * psi0 * Delta * beta_k / gamma1 *
			(1.0 + 0.5 * alpha1 * beta_k);
		r_K = 0.5 * (Consts::k - 1.0) * B * pow(1.0 - z0, 2.0);
	}

	void calcOtherThings()
	{
		Lambda_D = Lambda_K / eta_K;
		r_D = hi1 - (hi1 - r_K) * pow((1.0 - Consts::alpha_k * Delta + Lambda_K) / (1.0 - Consts::alpha_k * Delta + Lambda_D), Consts::k - 1.0);

		omega = q * Consts::K / (2 * Consts::f / (Consts::k - 1.0) * r_D / pow(vd, 2.0) - 1.0 / 3.0);
		W0 = omega / Delta;
		double S = 0.25 * Consts::pi * pow(d, 2.0) * ns;
		L0 = W0 / S;
		LD = Lambda_D * L0;
		fi = Consts::K + 1.0 / 3.0 * omega / q;
		Ik = sqrt(Consts::f * omega * fi * q * B) / S;
	}
};
typedef std::vector<Barrel> Barrels;

#endif

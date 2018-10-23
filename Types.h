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
	// ����� ��������� ��� ������� ����
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
	// ��� �������� ����������
	const double C = 1e6;
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

	Analog(const std::string &_name, double _d, double _q, double _vd) :
		name(_name), d(_d), q(_q), vd(_vd)
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

		std::fstream file("files/chuev.txt", std::ios_base::in);
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

struct OutputResult
{
	double eta_K;
	double Lambda_D;
	double r_D, r_Dmin;
	double omega;
	double W0;
	double L0;
	double LD;
	double Ik;
};
typedef std::vector<OutputResult> OutputResults;

struct Barrel
{
	// �������� ������
	double d, q, vd;
	// �������������� ��������
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
	// ��� ����
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

	OutputResult o_r;
	OutputResults o_rs;

	std::vector<double> Z_Sluh;

	/*std::vector<double> vEta_K;
	double Lambda_D;
	double r_D, r_Dmin;
	double omega;
	double W0;
	double L0;
	double LD;
	double Ik;*/

	Barrel(const StartData &data = StartData()) :
		d(data.d), q(data.q), vd(data.vd), p0(data.p0), K(data.K), ns(data.ns)
	{
		Cq = q / pow(d * 10.0, 3.0);
		CE = q * pow(vd, 2.0) / (2e3 * Consts::g * pow(d * 10.0, 3.0));
		o_r.r_Dmin = 1.0 / 6.0 * pow(vd, 2.0) / (Consts::f / (Consts::k - 1.0));
	}

	Barrel(const Barrel &barr, const StartData &data = StartData()) :
		d(data.d), q(data.q), vd(data.vd), p0(data.p0), K(data.K), ns(data.ns)
	{
		Cq = barr.Cq;
		CE = barr.CE;
		CE15 = barr.CE15;
		eta_omega = barr.eta_omega;
		o_r.omega = q * barr.omega_q;
		omega_q = barr.omega_q;
		pm_kr = barr.pm_kr;
		pm = barr.pm;
		hi = barr.hi;
		ns = barr.ns;
	}

	void calcForTest(double _B)
	{
		B = _B;
		Delta = 500;
		pm = 240e6;
		p_mid = 0.5 * pm;
		double eta_K = 0.5;

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
		calcForEtaK(eta_K);
	}

	int calcB(double a, double b)
	{
		p_mid = 0.5 * pm;
		hi1 = 1 - (Consts::k - 1) * p_mid / Consts::f * (1 - Consts::alpha_k * Consts::delta) / Consts::delta;

		psi0 = (1 / Delta - 1 / Consts::delta) / (Consts::f / p0 - (1 - Consts::alpha_k * Consts::delta) / Consts::delta);
		sigma0 = sqrt(1 + 4.0 * Consts::lambda / Consts::kapa * psi0);
		z0 = 2 * psi0 / (Consts::kapa * (1 + sigma0));
		K1 = Consts::kapa * sigma0;

		if (halfSearchB(a, b) == 1) return 1;
		return 0;
	}

	int halfSearchB(double a, double b)
	{
		double pi_m_star = pm / p0;

		int i = 0;
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

			if (fabs(pi_m - pi_m_star) < eps) break;					// ������� ������ �� �����
			
			if (pi_m > pi_m_star) a = B;
			else b = B;

			i++;
			if (i > 10000) return 1;
		}

		return 0;
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

	int calcForEtaK(double eta_K)
	{
		o_r.Lambda_D = Lambda_K / eta_K;
		o_r.r_D = hi1 - (hi1 - r_K) * pow((1.0 - Consts::alpha_k * Delta + Lambda_K) /
			(1.0 - Consts::alpha_k * Delta + o_r.Lambda_D), Consts::k - 1.0);

		if (o_r.r_D > o_r.r_Dmin && o_r.Lambda_D >= 1.0 && o_r.Lambda_D <= 10.0)
		{
			o_r.omega = q * K / (2 * Consts::f / (Consts::k - 1.0) * o_r.r_D / pow(vd, 2.0) - 1.0 / 3.0);
			o_r.W0 = o_r.omega / Delta;
			double S = 0.25 * pi * pow(d, 2.0) * ns;
			o_r.L0 = o_r.W0 / S;
			o_r.LD = o_r.Lambda_D * o_r.L0;
			fi = K + 1.0 / 3.0 * o_r.omega / q;
			o_r.Ik = sqrt(Consts::f * o_r.omega * fi * q * B) / S;

			// ��������� � ���������
			o_rs.push_back(o_r);

			return 1;
		}

		return 0;
	}

	void calcZSluh(double v1_vd) {
		Z_Sluh.push_back(Consts::C * sqrt(1.0 / o_r.Lambda_D + 1.0) / (pow(o_r.omega / q, 1.5) * pow(o_r.LD / d, 4.0) * v1_vd));
	}
};
typedef std::vector<Barrel> Barrels;

#endif

#ifndef TYPES_H
#define TYPES_H

#include <math.h>
#include <vector>
#include <string>

namespace Consts {
	// Общие константы для решения ОЗВБ
	const double g = 9.80655;
	const double K = 1.03;
	const double fi1 = 1.02;
	const double k = 1.25;
	const double Nkr = 1.23;
	const double f = 1e6;
	const double alpha_k = 1e-3;
	const double delta = 1600;
	const double kapa = 1;
	const double lambda = 0;

	const double p_flash = 1e6;
	const double p0 = 1e7;
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
	double K;
	double fi;
	double Nkr;
	double pm_kr;
	double pm;
	double p0;
	double l_d;
	double hi1;
	double ns;
	double p_mid;

	Barrel(double _d = 0.1, double _q = 4.35, double _vd = 1600, double _ns = 1) :
		d(_d), q(_q), vd(_vd), ns(_ns)
	{
		Cq = q / pow(d * 10, 3.0);
		CE = q * pow(vd, 2.0) / (2e3 * Consts::g * pow(d, 3.0)) * 1e-3;
	}
};

#endif

#ifndef TYPES_H
#define TYPES_H

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
}

struct Analog
{
	std::string name;
	double d;
	double q;
	double vd;

	double CE15;
	double p;
};
typedef std::vector<Analog> Analogs;

struct Chuev
{
	double CE;
	double p;
	double eta_omega;

	Chuev(double _CE = 0, double _p = 0, double _eta_omega = 0) :
		CE(_CE), p(_p), eta_omega(_eta_omega) {}

	friend std::ostream& operator<<(std::ostream &os, const Chuev &chuev);
};

inline std::ostream& operator<<(std::ostream &os, const Chuev &chuev)
{
	os << "CE = " << chuev.CE << '\t'
		<< "p = " << chuev.p << '\t' << "eta_omega = " << chuev.eta_omega;
	return os;
}
#endif

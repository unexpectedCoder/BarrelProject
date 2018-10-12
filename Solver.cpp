#include "Solver.h"
#include "Parser.h"

#include <iostream>
#include <math.h>

using namespace std;

void Solver::fillAnalogs(const string &path)
{
	Analog a;
	char ch;
	cout << "\tAdd analog? (+/-): ";
	cin >> ch;

	while (ch == '+')
	{
		cout << "\t - name: "; cin >> a.name;
		cout << "\t - d, mm: "; cin >> a.d;
		cout << "\t - q, kg: "; cin >> a.q;
		cout << "\t - Vd, mps: "; cin >> a.vd;

		analogs.push_back(a);
		num_of_as++;

		cout << "\tAdd analog? (+/-): "; cin >> ch;
	}

	// Запись в файл
	Parser::createFile(path);
	Parser p(path);
	for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); itr++)
	{
		p.write(itr->name + "\t");
		p.write(itr->d, '\t');
		p.write(itr->q, '\t');
		p.write(itr->vd, '\n');
	}
}

Analogs& Solver::calcAnalogs(const string &path)
{
	Parser p(path);
	Analog a;
	double CE;
	char ch;

	analogs.clear();
	while (true)
	{
		a.name = p.readStr();
		a.d = p.readNext();
		a.q = p.readNext();
		a.vd = p.readNext();

		if (p.isEnd()) break;

		double Cq = a.q / pow(a.d * 10, 3);
		CE = a.q * pow(a.vd, 2) / (2e3 * Consts::g * pow(a.d, 3)) * 1e-3;

		do
		{
			double eta = linInterp(CE).eta_omega;
			a.CE15 = calcCE15(Cq, CE, eta);
			double pm_kr = linInterp(a.CE15).pm_kr;
			double omega = a.q * pow(a.vd, 2.0) / (2.0 * eta * Consts::g * 1e3);
			double fi = Consts::K + 1.0 / 3.0 * omega / a.q;
			a.pm = pm_kr * fi * Consts::Nkr / (Consts::fi1 + 0.5 * omega / a.q) * 1.2;

			cout << "\tAnalog's p(CE_15) (in tech system) for " << a.name << ":\n" <<
				"\t\tCE15 = " << a.CE15 << ", pm = " << a.pm << endl;
			
			cout << "\tClarify the solution? (+/-): ";
			cin >> ch;
			if (ch == '+') CE = a.CE15;
		} while (ch == '+');
		cout << '\n';

		analogs.push_back(a);
	}

	Parser::createFile("p_CE15.txt");
	if (p.open("p_CE15.txt"))
		for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); itr++)
		{
			p.write(itr->CE15, '\t');
			p.write(itr->pm, '\n');
		}

	return analogs;
}

Barrel& Solver::calcBarrelPressure(const string &path)
{
	double CE = barr.CE;		// Буфер для хранения CE, чтобы можно было рекурсивно искать pm
	char ch;
	do
	{
		barr.eta_omega = linInterp(CE).eta_omega;
		barr.CE15 = calcCE15(barr.Cq, CE, barr.eta_omega);
		barr.pm_kr = linInterp(barr.CE15).pm_kr;
		barr.omega_q = pow(barr.vd, 2.0) / (2e3 * Consts::g * barr.eta_omega);
		barr.fi = Consts::K + 1.0 / 3.0 * barr.omega_q;
		barr.pm = barr.pm_kr * barr.fi * Consts::Nkr / (Consts::fi1 + 0.5 * barr.omega_q) * 1.2;

		cout << "\t CE15 = " << barr.CE15 << ", pm = " << barr.pm << endl;
		cout << "\tClarify the solution? (+/-): ";
		cin >> ch;
		if (ch == '+') CE = barr.CE15;
	} while (ch == '+');

	double pm_nround = barr.pm;				// Хранит неокругленное значение pm
	barr.pm /= Consts::g / 1e6;				// Перевод в Па

	cout << "\n\tFor your own sample pm = " << barr.pm * 1e-6 << " MPa\n";
	cout << "\t - please, enter rounded pm, MPa: pm = ";
	cin >> barr.pm;
	barr.pm *= 1e6;										// Перевод в Па
	// Расчет остальных характеристик
	barr.p_mid = 0.5 * barr.pm;
	barr.hi1 = 1 - (Consts::k - 1) * barr.p_mid / Consts::f *
		(1 - Consts::alpha_k * Consts::delta) / Consts::delta;

	// Запись в файл в виде таблички
	makeTableTxt(pm_nround, path);

	return barr;
}

Chuev Solver::linInterp(double CE, const string &path_chuev)
{
	Parser par(path_chuev.c_str());
	double CE1, CE2, p1, p2, eta1, eta2;

	while (!par.isEnd())
	{
		CE1 = par.readNext();
		p1 = par.readNext();
		eta1 = par.readNext();

		if (CE - CE1 < 100)
		{
			CE2 = par.readNext();
			p2 = par.readNext();
			eta2 = par.readNext();

			break;
		}
	}

	return Chuev(CE,
		p1 + (CE - CE1)/(CE2 - CE1) * (p2 - p1),
		eta1 + (CE - CE1) / (CE2 - CE1) * (eta2 - eta1));
}

double Solver::calcCE15(double cq, double ce, double eta)
{
	return 0.5 * 15.0 / cq * (ce - 3.0 * eta * cq + sqrt(pow(ce - 3.0 * eta * cq, 2.0) + 4.0 / 5.0 * ce * eta * pow(cq, 2.0)));
}

void Solver::makeTableTxt(double pm_nround, const string &path)
{
	Parser::createFile(path);
	Parser p(path);

	p.write("Cq: ");
	p.write(barr.Cq, '\t');
	p.write("w/q: ");
	p.write(barr.omega_q, '\t');
	p.write("pm_kr, atm: ");
	p.write(barr.pm_kr, '\t');
	p.write("p0, MPa: ");
	p.write(Consts::p0, '\n');

	p.write("CE: ");
	p.write(barr.CE, '\t');
	p.write("K: ");
	p.write(Consts::K, '\t');
	p.write("\tpm, atm: ");
	p.write(pm_nround, '\t');
	p.write("(l/d)max: ");
	p.write("???\n");

	p.write("CE15: ");
	p.write(barr.CE15, '\t');
	p.write("fi: ");
	p.write(barr.fi, '\t');
	p.write("pm, Mpa: ");
	p.write(pm_nround / Consts::g, '\t');
	p.write("hi: ");
	p.write(barr.hi1, '\n');

	p.write("ew: ");
	p.write(barr.eta_omega, '\t');
	p.write("Nkr: ");
	p.write(Consts::Nkr, '\t');
	p.write("pm_r, MPa: ");
	p.write(barr.pm * 1e-6, '\t');
	p.write("\tns: ");
	p.write(barr.ns, '\n');
}

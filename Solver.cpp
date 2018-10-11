#include "Solver.h"
#include "Parser.h"

#include <iostream>
#include <math.h>

using namespace std;

void Solver::fillAnalogs(const string &path)
{
	Analog a;
	char ch;
	cout << "\tAdd analog? (+/-): "; cin >> ch;

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
	Parser par(path);
	for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); itr++)
	{
		par.write(itr->name, '\t');
		par.write(itr->d, '\t');
		par.write(itr->q, '\t');
		par.write(itr->vd, '\n');
	}
}

Analogs& Solver::calcAnalogs(const string &path)
{
	Parser par(path);
	Analog a;
	double CE;
	char ch;

	analogs.clear();
	while (true)
	{
		a.name = par.readStr();
		a.d = par.readNext();
		a.q = par.readNext();
		a.vd = par.readNext();

		if (par.isEnd()) break;

		double Cq = a.q / pow(a.d * 10, 3);
		CE = a.q * pow(a.vd, 2) / (2e3 * Consts::g * pow(a.d, 3)) * 1e-3;

		do
		{
			double eta = linInterp(CE).eta_omega;
			a.CE15 = countCE15(Cq, CE, eta);
			double pm_kr = linInterp(a.CE15).p;
			double omega = a.q * pow(a.vd, 2.0) / (2.0 * eta * Consts::g * 1e3);
			double fi = Consts::K + 1.0 / 3.0 * omega / a.q;
			a.p = pm_kr * fi * Consts::Nkr / (Consts::fi1 + 0.5 * omega / a.q) * 1.2;

			cout << "\tAnalog's p(CE_15) (in tech system) for " << a.name << ":\n" <<
				"\t\tCE15 = " << a.CE15 << ", pm = " << a.p << endl;
			
			cout << "\tClarify the solution? (+/-): ";
			cin >> ch;
			if (ch == '+') CE = a.CE15;
		} while (ch == '+');
		cout << '\n';

		analogs.push_back(a);
	}

	par.close();
	Parser::createFile("p_CE15.txt");
	if (par.open("p_CE15.txt"))
		for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); itr++)
		{
			par.write(itr->CE15, '\t');
			par.write(itr->p, '\n');
		}

	return analogs;
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

double Solver::countCE15(double cq, double ce, double eta)
{
	return 0.5 * 15.0 / cq * (ce - 3.0 * eta * cq + sqrt(pow(ce - 3.0 * eta * cq, 2.0) + 4.0 / 5.0 * ce * eta * pow(cq, 2.0)));
}

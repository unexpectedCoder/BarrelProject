#include "Solver.h"
#include "Parser.h"

#include <cstdlib>
#include <iostream>
#include <math.h>

using namespace std;

int Solver::fcounter = 0;

void Solver::makeTest()
{
	Barrel barr(0.122, 21.76, 690, 3e7, 1.04);
	barr.pm = 240e6;
	barr.calcForTest(1.25);

	Parser::createFile("files/test.txt");
	Parser par("files/test.txt", 'w');

	par.write("pi_m* = ", barr.pm / barr.p0);
	par.write("p_mid = ", barr.p_mid * 1e-6);
	par.write("hi1 = ", barr.hi1);
	par.write("Delta = ", barr.Delta);
	par.write("psi0 = ", barr.psi0);
	par.write("sigma0 = ", barr.sigma0);
	par.write("z0 = ", barr.z0);
	par.write("K1 = ", barr.K1);
	par.write("B = ", barr.B);
	par.write("B1 = ", barr.B1);
	par.write("gamma1 = ", barr.gamma1);
	par.write("alpha1 = ", barr.alpha1);
	par.write("beta1 = ", barr.beta1);
	par.write("beta2 = ", barr.beta2);
	par.write("beta_m = ", barr.beta_m);
	par.write("beta_m1 = ", barr.beta_m_1);
	par.write("beta_m2 = ", barr.beta_m_2);
	par.write("fi_m = ", barr.fi_m);
	par.write("pi_m = ", barr.pi_m);

	par.write("beta_k = ", barr.beta_k);
	par.write("beta_k1 = ", barr.beta_k1);
	par.write("beta_k2 = ", barr.beta_k2);
	par.write("fi_k = ", barr.fi_k);
	par.write("Lambda_K = ", barr.Lambda_K);
	par.write("r_K = ", barr.r_K);

	par.write("Lambda_D = ", barr.Lambda_D);
	par.write("r_D = ", barr.r_D);
	par.write("omega = ", barr.omega);
	par.write("W0 = ", barr.W0);
	par.write("L0 = ", barr.L0);
	par.write("LD = ", barr.LD);
	par.write("Ik = ", barr.Ik * 1e-6);
}

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
	Parser p(path, 'w');
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
	Parser p(path, 'r');
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

	Parser::createFile("files/p_CE15.txt");
	if (p.open("files/p_CE15.txt", 'w'))
		for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); itr++)
		{
			p.write(itr->CE15, '\t');
			p.write(itr->pm, '\n');
		}

	return analogs;
}

void Solver::calcBarrelPressure(const string &path)
{
	cout << "\n\tSearching the max pressure pm for your own sample:\n";

	Barrel barr;
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

		cout << "\tNow it has\n";
		cout << "\tCE15 = " << barr.CE15 << ", pm = " << barr.pm << endl;
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
	makeTableTxt(barr, pm_nround, path);
	// Добавление в файл p_CE15.txt
	Parser p("files/p_CE15.txt", 'a');
	p.write(barr.CE15, '\t');
	p.write(barr.pm * 1e-6 * Consts::g, '\n');

	Parser::createFile("files/barrel_src.txt");
	p.open("files/barrel_src.txt", 'w');
	p.writeBarrel(barr);
}

Barrels& Solver::solveInvProblem()
{
	Parser par("files/barrel_src.txt", 'r');
	Barrel barr = par.readBarrel();

	fillData("Set the loading density (Delta), kg/m^3:", Delta);
	fillData("Set the dimesionless coordinate of the end of burning (eta_K):", eta_K);

	double a, b;
	cout << "\tSet the search range for B*:\n";
	cout << "\t - left: "; cin >> a;
	cout << "\t - right: "; cin >> b;

	Parser::createFile("files/B(Delta).txt");
	par.open("files/B(Delta).txt", 'w');

	for (vector<double>::iterator itr1 = Delta.begin(); itr1 != Delta.end(); itr1++)
	{
		barr.Delta = *itr1;
		barr.calcB(a, b);
		barr.calcLambdaK();

		for (vector<double>::iterator itr2 = eta_K.begin(); itr2 != eta_K.end(); itr2++)
		{
			barr.eta_K = *itr2;
			barr.calcForEtaK();

			if (barr.r_D > barr.r_Dmin && (barr.Lambda_D > 3 && barr.Lambda_D < 10))
			{
				barrs.push_back(barr);
				writeBarrelToFile(barr);
			}
		}

		par.write(barr.Delta, '\t');
		par.write(barr.B, '\n');
	}

	return barrs;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Вспомогательные private-функции ///////////////////////////////////
void Solver::makeTableTxt(const Barrel &barr, double pm_nround, const string &path)
{
	Parser::createFile(path);
	Parser p(path, 'w');

	p.write("Cq: ");
	p.write(barr.Cq, '\t');
	p.write("w/q: ");
	p.write(barr.omega_q, '\t');
	p.write("pm_kr, atm: ");
	p.write(barr.pm_kr, '\t');
	p.write("p0, MPa: ");
	p.write(barr.p0, '\n');

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

Chuev Solver::linInterp(double CE, const string &path)
{
	Parser par(path, 'r');
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
		p1 + (CE - CE1) / (CE2 - CE1) * (p2 - p1),
		eta1 + (CE - CE1) / (CE2 - CE1) * (eta2 - eta1));
}

double Solver::calcCE15(double cq, double ce, double eta)
{
	return 0.5 * 15.0 / cq * (ce - 3.0 * eta * cq + sqrt(pow(ce - 3.0 * eta * cq, 2.0) + 4.0 / 5.0 * ce * eta * pow(cq, 2.0)));
}

void Solver::fillData(const string &head_txt, vector<double> &data)
{
	if (!data.empty()) data.clear();

	double start, end, step;
	cout << '\t' << head_txt << '\n';
	cout << "\t - from: "; cin >> start;
	cout << "\t - to:  "; cin >> end;
	cout << "\t - step: "; cin >> step;

	while (start < end + 0.5 * step)
	{
		data.push_back(start);
		start += step;
	}
}

void Solver::writeBarrelToFile(const Barrel &barr)
{
	int lim = 99999;
	char buf[6];
	string path = "barrels/barrel_";

	fcounter++;
	if (fcounter > lim)
		throw "Runtime error in <writeBarrelsToFile()> func!";
	_itoa_s(fcounter, buf, 10);

	path += buf;
	path += ".txt";

	Parser::createFile(path);
	Parser par(path, 'w');

	par.write("Delta = ", barr.Delta);
	par.write("eta_K = ", barr.eta_K);
	par.write("\n");

	par.write("pi_m* = ", barr.pm / barr.p0);
	par.write("p_mid = ", barr.p_mid * 1e-6);
	par.write("hi1 = ", barr.hi1);
	par.write("psi0 = ", barr.psi0);
	par.write("sigma0 = ", barr.sigma0);
	par.write("z0 = ", barr.z0);
	par.write("K1 = ", barr.K1);
	par.write("B = ", barr.B);
	par.write("B1 = ", barr.B1);
	par.write("gamma1 = ", barr.gamma1);
	par.write("alpha1 = ", barr.alpha1);
	par.write("beta1 = ", barr.beta1);
	par.write("beta2 = ", barr.beta2);
	par.write("beta_m = ", barr.beta_m);
	par.write("beta_m1 = ", barr.beta_m_1);
	par.write("beta_m2 = ", barr.beta_m_2);
	par.write("fi_m = ", barr.fi_m);
	par.write("pi_m = ", barr.pi_m);

	par.write("beta_k = ", barr.beta_k);
	par.write("beta_k1 = ", barr.beta_k1);
	par.write("beta_k2 = ", barr.beta_k2);
	par.write("fi_k = ", barr.fi_k);
	par.write("Lambda_K = ", barr.Lambda_K);
	par.write("r_K = ", barr.r_K);

	par.write("Lambda_D = ", barr.Lambda_D);
	par.write("r_D = ", barr.r_D);
	par.write("omega = ", barr.omega);
	par.write("W0 = ", barr.W0);
	par.write("L0 = ", barr.L0);
	par.write("LD = ", barr.LD);
	par.write("Ik = ", barr.Ik * 1e-6);
}

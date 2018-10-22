#include "Solver.h"
#include "Parser.h"

#include <cstdlib>
#include <iostream>
#include <math.h>

using namespace std;

int Solver::fcounter = 0;

void Solver::makeTest()
{
	Barrel barr(StartData(0.122, 21.76, 690, 3e7, 1.05, 1.04));
	barr.calcForTest(1.25);

	Parser::createFile(TEST_PATH);
	Parser par(TEST_PATH, 'w');

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
	char ch;
	cout << "\tДобавить аналог? (+/-): ";
	cin >> ch;

	while (ch == '+')
	{
		string name;
		double d, q, vd;

		cout << "\tАналог №" << num_of_as << ":\n";
		cout << "\t - название: "; cin >> name;
		cout << "\t - d, мм: "; cin >> d;
		cout << "\t - q, кг: "; cin >> q;
		cout << "\t - Vd, м/с: "; cin >> vd;

		analogs.push_back(Analog(name, d, q, vd));
		num_of_as++;

		cout << "\tДобавить аналог? (+/-): "; cin >> ch;
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
	char ch;

	double K;
	cout << "\tКоэф-т учета второстепенных работ для аналогов:\n\t K = ";
	cin >> K;

	analogs.clear();
	while (true)
	{
		string name = p.readStr();
		double d = p.readNext() * 1e-3;
		double q = p.readNext();
		double vd = p.readNext();

		if (p.isEnd()) break;								// Выход из цикла. При таком расположении
																				// не произойдет работы с переменными, в которые записывается последняя (пустая)
																				// строка, из-за чего значения будут неадекватными.

		Analog a(name, d, q, vd);
		a.CE15 = a.CE;
		do
		{
			a.eta_omega = Chuev(a.CE15).eta_omega;
			a.CE15 = calcCE15(a.Cq, a.CE, a.eta_omega);
			double pm_kr = Chuev(a.CE15).pm_kr;
			double omega_q = a.CE15 / (a.Cq * a.eta_omega);
			double fi = K + 1.0 / 3.0 * omega_q;
			a.pm = pm_kr * fi * Consts::Nkr / (Consts::fi1 + 0.5 * omega_q) * NEW_BARR_K;

			cout << "\tp(CE_15) (в тех. системе) для " << a.name << ":\n" <<
				"\t\tCE15 = " << a.CE15 << ", pm = " << a.pm << endl;
			
			cout << "\tУточнить решение? (+/-): ";
			cin >> ch;
		} while (ch == '+');
		cout << '\n';

		analogs.push_back(a);
	}

	Parser::createFile(P_CE15_PATH);
	if (p.open(P_CE15_PATH, 'w'))
		for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); itr++)
		{
			p.write(itr->CE15, '\t');
			p.write(itr->pm, '\n');
		}

	return analogs;
}

void Solver::calcBarrelPressure(const string &path)
{
	cout << "\n\tПоиск макс. давления для собственного образца:\n";

	Barrel barr;
	char ch;

	barr.CE15 = barr.CE;
	do
	{
		barr.eta_omega = Chuev(barr.CE15).eta_omega;
		barr.CE15 = calcCE15(barr.Cq, barr.CE, barr.eta_omega);
		barr.pm_kr = Chuev(barr.CE15).pm_kr;
		barr.omega_q = barr.CE15 / (barr.Cq * barr.eta_omega);
		barr.fi = barr.K + 1.0 / 3.0 * barr.omega_q;
		barr.pm = barr.pm_kr * barr.fi * Consts::Nkr / (Consts::fi1 + 0.5 * barr.omega_q) * NEW_BARR_K;

		cout << "\tСейчас\n";
		cout << "\t CE15 = " << barr.CE15 << ", pm = " << barr.pm << endl;
		cout << "\tУточнить решение? (+/-): ";
		cin >> ch;
	} while (ch == '+');
	double pm_nround = barr.pm;				// Хранит неокругленное значение pm
	barr.pm /= Consts::g / 1e6;				// Перевод в Па

	cout << "\n\tДля вашего образца: pm = " << barr.pm * 1e-6 << " МПа\n";
	cout << "\t - введите округленное значение pm, МПа: pm = ";
	cin >> barr.pm;
	barr.pm *= 1e6;										// Перевод в Па
	// Расчет остальных характеристик
	barr.p_mid = 0.5 * barr.pm;
	barr.hi1 = 1 - (Consts::k - 1) * barr.p_mid / Consts::f *
		(1 - Consts::alpha_k * Consts::delta) / Consts::delta;

	// Запись в файл в виде таблички
	makeTableTxt(barr, pm_nround, path);
	// Добавление в файл p_CE15.txt
	Parser p(P_CE15_PATH, 'a');
	p.write(barr.CE15, '\t');
	p.write(barr.pm * 1e-6 * Consts::g, '\n');

	Parser::createFile(BARR_SRC_PATH);
	p.open(BARR_SRC_PATH, 'w');
	p.writeBarrel(barr);
}

Barrels& Solver::solveInvProblem()
{
	Parser par(BARR_SRC_PATH, 'r');
	Barrel barr = par.readBarrel();

	fillData("Задайте плотность заряжания (Delta), кг/м^3:", Delta);
	fillData("Задайте безразмерную координату конца горения (eta_K):", eta_K);

	double a, b;
	cout << "\tЗадайте границы поиска параметра Дроздова B*:\n";
	cout << "\t - слева: "; cin >> a;
	cout << "\t - справа: "; cin >> b;

	// Чтение из табл. 3.3 (Чуев)
	par.open(CHUEV_TABLE_33_PATH, 'r');
	double **table = new double*[ROWS];
	for (int i = 0; i < ROWS; i++)
	{
		*(table + i) = new double[COLOUMS];
		for (int j = 0; j < COLOUMS; j++)
			table[i][j] = par.readNext();
	}
	// Для контрольного вывода
	/*cout << "\nTable:\n";
	for (int i = 0; i < ROWS; i++)
	{
		for (int j = 0; j < COLOUMS; j++)
			cout << table[i][j] << '\t';
		cout << '\n';
	}*/

	Parser::createFile(B_DELTA_PATH);
	par.open(B_DELTA_PATH, 'w');
	Parser::createFile(Z_SLUH_PATH);
	Parser parsluh(Z_SLUH_PATH, 'w');

	parsluh.write(0.0);
	for (vector<double>::iterator itr = eta_K.begin(); itr != eta_K.end(); itr++)
		parsluh.write('\t', *itr);

	for (vector<double>::iterator itr1 = Delta.begin(); itr1 != Delta.end(); itr1++)
	{
		barr.Delta = *itr1;
		// Поиск к-та Дроздова B* и относительной длины процесса горения Lambda_K
		barr.calcB(a, b);
		barr.calcLambdaK();
		// Действия с 22 по 34
		parsluh.write('\n', *itr1);
		for (vector<double>::iterator itr2 = eta_K.begin(); itr2 != eta_K.end(); itr2++)
		{
			barr.eta_K = *itr2;
			barr.calcForEtaK();

			if (barr.r_D > barr.r_Dmin && (barr.Lambda_D > 3 && barr.Lambda_D < 10))
			{
				double hi_n = 1.0 / (1.0 / Chuev(barr.CE15).hi + 0.75 * barr.d / barr.L0);
				barr.calcZSluh(belinearInterp(hi_n, barr.Lambda_D, table));

				barrs.push_back(barr);
				writeBarrelToFile(barr);

				parsluh.write('\t', barr.Z_Sluh);
			}
			else
				parsluh.write('\t', 0.0);
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
	p.write(barr.K, '\t');
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

double Solver::calcCE15(double cq, double ce, double eta) {
	return 0.5 * 15.0 / cq * (ce - 3.0 * eta * cq + sqrt(pow(ce - 3.0 * eta * cq, 2.0) + 4.0 / 5.0 * ce * eta * pow(cq, 2.0)));
}

void Solver::fillData(const string &head_txt, vector<double> &data)
{
	if (!data.empty()) data.clear();

	double start, end, step;
	cout << '\t' << head_txt << '\n';
	cout << "\t - от: "; cin >> start;
	cout << "\t - до: "; cin >> end;
	cout << "\t - шаг: "; cin >> step;

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

double Solver::belinearInterp(double x, double y, double **table)
{
	// Массив значений по оси ROWS
	double *arrx = new double[ROWS];
	for (int i = 0; i < ROWS; i++)
		arrx[i] = table[i][0];
	// Массив значений по оси COLOUMS
	double *arry = new double[COLOUMS];
	for (int j = 0; j < COLOUMS; j++)
		arry[j] = table[0][j];

	// Поиск индексов эл-ов, между которыми находится заданная точка...
	// ...по оси ROWS
	int ix, iy;
	for (int i = 1; i < ROWS - 1; i++)
		if (x >= arrx[i] && x <= arrx[i + 1])
		{
			ix = i;
			break;
		}
	// ...по оси COLOUMS
	for (int j = 1; j < COLOUMS; j++)
		if (y - arry[j] <= 1.0)
		{
			iy = j;
			break;
		}

	// Расчетная часть
	double fQ11 = table[ix][iy];
	double fQ21 = table[ix][iy + 1];
	double fQ12 = table[ix + 1][iy];
	double fQ22 = table[ix + 1][iy + 1];

	double fR1 = 1.0 / (arry[iy + 1] - arry[iy]) * ((arry[iy + 1] - y) * fQ11 + (y - arry[iy]) * fQ21);
	double fR2 = 1.0 / (arry[iy + 1] - arry[iy]) * ((arry[iy + 1] - y) * fQ12 + (y - arry[iy]) * fQ22);

	return 1.0 / (arrx[ix + 1] - arrx[ix]) * ((arrx[ix + 1] - x) * fR1 + (x - arrx[ix]) * fR2);
}

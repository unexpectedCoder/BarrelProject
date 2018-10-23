#include "Solver.h"
#include "Parser.h"

#include <cstdlib>
#include <math.h>

using namespace std;

int AnaliticSolver::file_count = 0;

double Solver::belinearInterp(double x, double y, double **table)
{
	int ix, iy;

	/*
		* ������� table[0][0] �� ���������������, �.�. ����� ������ ��� �������� ������� (��. ���� chuev_table_33.txt.).
		* ��� �������� � ����, ��� � ������ i ���������� � 1, � �� � 0.
	*/
	// ����� �������� ��-��, ����� �������� ��������� �������� �����...
	// ...�� ��� COLOUMS (X)
	for (int i = 1; i < COLOUMS; i++)
		if (x - table[0][i] <= 1.0)
		{
			ix = i;
			break;
		}
	// ...�� ��� ROWS (Y)
	for (int i = 1; i < ROWS - 1; i++)
		if (y >= table[i][0] && y <= table[i + 1][0])
		{
			iy = i;
			break;
		}

	// ��������� �����
	double fQ11 = table[iy][ix];					// ���������������� �������� - ��������� ������������
	double fQ21 = table[iy][ix + 1];			// �������� ������ � ������:
	double fQ12 = table[iy + 1][ix];			// 2�� ������ (�������) ������������� ��� X;
	double fQ22 = table[iy + 1][ix + 1];	// 1�� ������ (������) - ��� Y.

	double fR1 = 1.0 / (table[0][ix + 1] - table[0][ix]) * ((table[0][ix + 1] - x) * fQ11 + (x - table[0][ix]) * fQ21);
	double fR2 = 1.0 / (table[0][ix + 1] - table[0][ix]) * ((table[0][ix + 1] - x) * fQ12 + (x - table[0][ix]) * fQ22);

	return 1.0 / (table[iy + 1][0] - table[iy][0]) * ((table[iy + 1][0] - y) * fR1 + (y - table[iy][0]) * fR2);
}

int TestSolver::solve()
{
	cout << "\t<������� ��������� ��������>\n";

	barr->calcForTest(1.25);

	Parser::createFile(TEST_PATH);
	Parser par(TEST_PATH, 'w');

	par.write("pi_m* = ", barr->pm / barr->p0);
	par.write("p_mid = ", barr->p_mid * 1e-6);
	par.write("hi1 = ", barr->hi1);
	par.write("Delta = ", barr->Delta);
	par.write("psi0 = ", barr->psi0);
	par.write("sigma0 = ", barr->sigma0);
	par.write("z0 = ", barr->z0);
	par.write("K1 = ", barr->K1);
	par.write("B = ", barr->B);
	par.write("B1 = ", barr->B1);
	par.write("gamma1 = ", barr->gamma1);
	par.write("alpha1 = ", barr->alpha1);
	par.write("beta1 = ", barr->beta1);
	par.write("beta2 = ", barr->beta2);
	par.write("beta_m = ", barr->beta_m);
	par.write("beta_m1 = ", barr->beta_m_1);
	par.write("beta_m2 = ", barr->beta_m_2);
	par.write("fi_m = ", barr->fi_m);
	par.write("pi_m = ", barr->pi_m);

	par.write("beta_k = ", barr->beta_k);
	par.write("beta_k1 = ", barr->beta_k1);
	par.write("beta_k2 = ", barr->beta_k2);
	par.write("fi_k = ", barr->fi_k);
	par.write("Lambda_K = ", barr->Lambda_K);
	par.write("r_K = ", barr->r_K);

	par.write("Lambda_D = ", barr->o_r.Lambda_D);
	par.write("r_D = ", barr->o_r.r_D);
	par.write("omega = ", barr->o_r.omega);
	par.write("W0 = ", barr->o_r.W0);
	par.write("L0 = ", barr->o_r.L0);
	par.write("LD = ", barr->o_r.LD);
	par.write("Ik = ", barr->o_r.Ik * 1e-6);

	return 0;
}

int AnalogsSolver::solve()
{
	cout << "\t<������� ������� ��������>\n";

	char ch;
	cout << "\t������ ������ �� ��������? (+/-):";
	cin >> ch;
	if (ch == '+')
		fillAnalogsData();

	Parser par(ANALOGS_PATH, 'r');
	
	double K;
	cout << "\t - ����-� ����� �������������� ����� ��� ��������: K = ";
	cin >> K;

	analogs.clear();
	while (true)
	{
		string name = par.readStr();
		double d = par.readNext() * 1e-3;
		double q = par.readNext();
		double vd = par.readNext();

		if (par.isEnd()) break;							// ����� �� �����. ��� ����� ������������
																				// �� ���������� ������ � �����������, � ������� ������������ ��������� (������)
																				// ������, ��-�� ���� �������� ����� �������������.

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

			cout << "\tp(CE_15) (� ���. �������) ��� " << a.name << ":\n" <<
				"\t\tCE15 = " << a.CE15 << ", pm = " << a.pm << endl;
			
			cout << "\t�������� �������? (+/-): ";
			cin >> ch;
		} while (ch == '+');
		cout << '\n';

		analogs.push_back(a);
	}

	Parser::createFile(P_CE15_PATH);
	par.open(P_CE15_PATH, 'w');
	if (!par.isOpen())
	{
		status = "faild";
		return 1;
	}

	for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); itr++)
		{
			par.write(itr->CE15, '\t');
			par.write(itr->pm, '\n');
		}

	printResults();
	return 0;
}

void AnalogsSolver::fillAnalogsData()
{
	cout << "\t<������� ���������� ������ �� ��������>\n";

	char ch;
	cout << "\t\t�������� ������? (+/-): ";
	cin >> ch;

	int num = 0;
	while (ch == '+')
	{
		string name;
		double d, q, vd;

		cout << "\t\t������ �" << ++num << ":\n";
		cout << "\t\t - ��������: "; cin >> name;
		cout << "\t\t - d, ��: "; cin >> d;
		cout << "\t\t - q, ��: "; cin >> q;
		cout << "\t\t - Vd, �/�: "; cin >> vd;

		analogs.push_back(Analog(name, d, q, vd));

		cout << "\t\t�������� ������? (+/-): "; cin >> ch;
	}

	// ������ � ����
	Parser::createFile(ANALOGS_PATH);
	Parser p(ANALOGS_PATH, 'w');
	for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); itr++)
	{
		p.write(itr->name + "\t");
		p.write(itr->d, '\t');
		p.write(itr->q, '\t');
		p.write(itr->vd, '\n');
	}
}

void AnaliticSolver::calcMaxPressure()
{
	cout << "\t<������� ������ ������ ����. ��������>\n";

	Barrel b;
	char ch;

	b.CE15 = b.CE;
	do
	{
		b.eta_omega = Chuev(b.CE15).eta_omega;
		b.CE15 = calcCE15(b.Cq, b.CE, b.eta_omega);
		b.pm_kr = Chuev(b.CE15).pm_kr;
		b.omega_q = b.CE15 / (b.Cq * b.eta_omega);
		b.fi = b.K + 1.0 / 3.0 * b.omega_q;
		b.pm = b.pm_kr * b.fi * Consts::Nkr / (Consts::fi1 + 0.5 * b.omega_q) * NEW_BARR_K;

		cout << "\t\t������ " << "CE15 = " << b.CE15 << ", pm = " << b.pm << endl;
		cout << "\t\t�������� �������? (+/-): ";
		cin >> ch;
	} while (ch == '+');

	double pm_nround = b.pm;		// ������ ������������� �������� pm
	b.pm *= Consts::g * 1e4;		// ������� � �� �� ���. ���.

	cout << "\n\t\t��� ������ �������: pm = " << b.pm * 1e-6 << " ���\n";
	cout << "\t\t - ������� ����������� �������� pm, ���: pm = ";
	cin >> b.pm;
	b.pm *= 1e6;								// ������� � �� �� ���

	// ������ ��������� �������������
	b.p_mid = 0.5 * b.pm;
	b.hi = Chuev(b.CE15).hi;

	// ������ � ���� � ���� ��������
	makeTableTxt(b, pm_nround, BARR_TABLE_PATH);
	// ���������� � ���� p_CE15.txt
	Parser p(P_CE15_PATH, 'a');
	p.write(b.CE15, '\t');
	p.write(b.pm / Consts::g * 1e-4, '\n');
	// ������ ������ � ������ � ���� ��� ����������� ����������
	Parser::createFile(BARR_SRC_PATH);
	p.open(BARR_SRC_PATH, 'w');
	p.writeBarrel(b);
}

int AnaliticSolver::solve()
{
	cout << "\t<������� �������������� �������>\n";

	// ������ ������ � ������ �� �����-���������
	Parser par(BARR_SRC_PATH, 'r');
	barr = new Barrel(par.readBarrel());

	// ������� �������� ��������� ��������� � eta_K
	fillBarrelData("������� ��������� ��������� (Delta), ��/�^3:", Delta);
	fillBarrelData("������� ������������ ���������� ����� ������� (eta_K):", eta_K);

	// ������ �� ����. 3.3 (����)
	par.open(CHUEV_TABLE_33_PATH, 'r');
	double **table = new double*[ROWS];
	for (int i = 0; i < ROWS; i++)
	{
		*(table + i) = new double[COLOUMS];
		for (int j = 0; j < COLOUMS; j++)
			table[i][j] = par.readNext();
	}
	// ��� ������������ ������
	/*cout << "\nTable:\n";
	for (int i = 0; i < ROWS; i++)
	{
		for (int j = 0; j < COLOUMS; j++)
			cout << table[i][j] << '\t';
		cout << '\n';
	}*/

	// �������� � ���������� ������
	Parser::createFile(B_DELTA_PATH);
	par.open(B_DELTA_PATH, 'w');
	Parser::createFile(Z_SLUH_PATH);
	Parser parsluh(Z_SLUH_PATH, 'w');
	// ���������� ������� ������ �����
	parsluh.write(0.0);
	for (vector<double>::iterator itr = eta_K.begin(); itr != eta_K.end(); itr++)
		parsluh.write('\t', *itr);

	// ��� ������ ��������� �������� B*
	double a, b;
	cout << "\t\t������� ������� ������ ��������� �������� B*:\n";
	cout << "\t\t - �����: "; cin >> a;
	cout << "\t\t - ������: "; cin >> b;
	if (a > b)
	{
		double buf = a;
		a = b;
		b = buf;
	}

	// ��������� �����
	for (vector<double>::iterator itr1 = Delta.begin(); itr1 != Delta.end(); itr1++)
	{
		barr->o_rs.clear();

		barr->Delta = *itr1;
		// ����� �-�� �������� B* � ������������� ����� �������� ������� Lambda_K
		barr->calcB(a, b);
		barr->calcLambdaK();

		// �������� � 22 �� 34
		parsluh.write('\n', *itr1);
		bool flag = true;
		for (vector<double>::iterator itr2 = eta_K.begin(); itr2 != eta_K.end(); itr2++)
		{
			barr->o_r.eta_K = *itr2;

			if (barr->calcForEtaK(*itr2) == 1)
			{
				double hi_n = 1.0 / (1.0 / Chuev(barr->CE15).hi + 0.75 * barr->d / barr->o_r.L0);
				barr->calcZSluh(belinearInterp(barr->o_r.Lambda_D, hi_n, table));

				parsluh.write('\t', *(barr->Z_Sluh.end() - 1));
			}
			else
			{
				parsluh.write('\t', 0.0);
				flag = false;
			}
		}

		if (flag)
			barrs.push_back(*barr);

		// ������ � �����
		writeBarrelToFile(*barr);
		
		par.write(barr->Delta, '\t');
		par.write(barr->B, '\n');
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// ��������������� private-������� ///////////////////////////////////
void AnaliticSolver::makeTableTxt(const Barrel &barr, double pm_nround, const string &path)
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
	p.write(48, '\n');

	p.write("CE15: ");
	p.write(barr.CE15, '\t');
	p.write("fi: ");
	p.write(barr.fi, '\t');
	p.write("pm, MPa: ");
	p.write(pm_nround * Consts::g * 1e-2, '\t');
	p.write("hi: ");
	p.write(barr.hi, '\n');

	p.write("ew: ");
	p.write(barr.eta_omega, '\t');
	p.write("Nkr: ");
	p.write(Consts::Nkr, '\t');
	p.write("pm_r, MPa: ");
	p.write(barr.pm * 1e-6, '\t');
	p.write("\tns: ");
	p.write(barr.ns, '\n');
}

void AnaliticSolver::fillBarrelData(const string &head_txt, vector<double> &data)
{
	if (!data.empty()) data.clear();

	double start, end, step;
	cout << "\t\t" << head_txt << '\n';
	cout << "\t\t - ��: "; cin >> start;
	cout << "\t\t - ��: "; cin >> end;
	cout << "\t\t - ���: "; cin >> step;

	while (start < end + 0.5 * step)
	{
		data.push_back(start);
		start += step;
	}
}

void AnaliticSolver::writeBarrelToFile(const Barrel &barr)
{
	int lim = 99999;
	char buf[6];
	string path = "barrels/barrel_";

	file_count++;
	if (file_count > lim)
		throw "Runtime error in <writeBarrelsToFile()> func!";
	_itoa_s(file_count, buf, 10);

	path += buf;
	path += ".txt";

	Parser::createFile(path);
	Parser par(path, 'w');

	par.write("Delta = ", barr.Delta);
	par.write("\nB* = ", barr.B);
	par.write("\nLambda_K = ", barr.Lambda_K);
	par.write("\n");
	
	for (OutputResults::const_iterator itr = barr.o_rs.begin(); itr != barr.o_rs.end(); itr++)
	{
		par.write("eta_K = ", itr->eta_K, '\t');
		par.write("Lambda_D = ", itr->Lambda_D, '\t');
		par.write("omega = ", itr->omega, '\t');
		par.write("W0 = ", itr->W0, '\t');
		par.write("LD = ", itr->LD, '\t');
		par.write("Ik = ", itr->Ik * 1e-6);
	}

	// ��� ���������� ������
	/*par.write("pi_m* = ", barr.pm / barr.p0);
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
	par.write("Ik = ", barr.Ik * 1e-6);*/
}

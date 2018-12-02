#include "Solver.h"

#include <algorithm>
#include <cstdlib>
#include <math.h>

using namespace std;

double Solver::bilinearInterp(double x, double y, double **table)
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

void TestSolver::solve()
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

	par.write("Lambda_D = ", barr->Lambda_D);
	par.write("r_D = ", barr->r_D);
	par.write("omega = ", barr->omega);
	par.write("W0 = ", barr->W0);
	par.write("L0 = ", barr->L0);
	par.write("LD = ", barr->LD);
	par.write("Ik = ", barr->Ik * 1e-6);
}

void AnalogsSolver::solve()
{
	cout << "\t<������� ������� ��������>\n";

	char ch;
	cout << "\t������ ������ �� ��������? (+/-):";
	cin >> ch;
	if (ch == '+')
		fillAnalogsData();

	Parser par(ANALOGS_PATH, 'r');

	double K;
	cout << "\t - ����-� ����� �������������� �����: K = ";
	cin >> K;

	analogs.clear();
	while (true)
	{
		string name = par.readStr();
		double d = par.readNext() * 1e-3;
		double q = par.readNext();
		double vd = par.readNext();
		double l_d = par.readNext();
		double pm = par.readNext();

		if (par.isEnd()) break;							// ����� �� �����. ��� ����� ������������
																				// �� ���������� ������ � �����������, � ������� ������������ ��������� (������)
																				// ������, ��-�� ���� �������� ����� �������������.

		Analog a(name, d, q, vd, l_d, pm);
		a.CE15 = a.CE;
		a.eta_omega = Chuev(a.CE15).eta_omega;
		do
		{
			a.CE15 = calcCE15(a.Cq, a.CE, a.eta_omega);
			a.eta_omega = Chuev(a.CE15).eta_omega;

			cout << "\t������ " << a.name << ":\n";
			cout << "\t\t - �������� �������? (+/-): ";
			cin >> ch;
		} while (ch == '+');

		analogs.push_back(a);
		cout << '\n';
	}

	Parser::createFile(ANALOGS_RES_PATH);
	Parser parres(ANALOGS_RES_PATH, 'w');
	Parser::createFile(P_CE15_PATH);
	par.open(P_CE15_PATH, 'w');
	if (!par.isOpen() || !parres.isOpen())
	{
		status = "faild";
		return;
	}

	// ������ �������� �������� � ����
	for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); itr++)
	{
		parres.write(itr->name + "\t");
		parres.write("d = ", itr->d, '\t');
		parres.write("q = ", itr->q, '\t');
		parres.write("Vd = ", itr->vd, '\t');
		parres.write("Cq = ", itr->Cq, '\t');
		parres.write("CE = ", itr->CE, '\t');
		parres.write("CE15 = ", itr->CE15, '\t');
		parres.write("pm = ", itr->pm, '\t');
		parres.write("ld = ", itr->l_d, '\n');

		par.write(itr->CE15, '\t');
		par.write(itr->pm, '\n');
	}

	printOutro();
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
		double d, q, vd, l_d, pm;

		cout << "\t\t������ �" << ++num << ":\n";
		cout << "\t\t - ��������: "; cin >> name;
		cout << "\t\t - d, ��: "; cin >> d;
		cout << "\t\t - q, ��: "; cin >> q;
		cout << "\t\t - Vd, �/�: "; cin >> vd;
		cout << "\t\t - ����� ������� ����� � ��������: "; cin >> l_d;
		cout << "\t\t - ������������ ��������, ��/��^2: "; cin >> pm;

		analogs.push_back(Analog(name, d, q, vd, l_d, pm));

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
		p.write(itr->vd, '\t');
		p.write(itr->l_d, '\t');
		p.write(itr->pm, '\n');
	}
}

void AnaliticSolver::calcMaxPressure()
{
	cout << "\t<������� ������ ������ ����. ��������>\n";

	Barrel b;
	char ch;

	b.CE15 = b.CE;
	b.eta_omega = Chuev(b.CE15).eta_omega;
	do
	{
		b.CE15 = calcCE15(b.Cq, b.CE, b.eta_omega);
		b.eta_omega = Chuev(b.CE15).eta_omega;				// �������� eta_omega ��� ������ CE15
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
	b.calcHi1();
	b.hi = Chuev(b.CE15).hi;

	Parser::createFile(BARR_LOG_PATH);
	Parser par(BARR_LOG_PATH, 'w');
	par.write("p_mid = ", b.p_mid, '\n');
	par.write("hi1 = ", b.hi1, '\n');
	// ������ � ���� � ���� ��������
	makeTableTxt(b, pm_nround, BARR_TABLE_PATH);
	// ���������� � ���� p_CE15.txt
	par.open(P_CE15_PATH, 'a');
	par.write(b.CE15, '\t');
	par.write(b.pm / Consts::g * 1e-4, '\n');
	// ������ ������ � ������ � ���� ��� ����������� ����������
	Parser::createFile(BARR_SRC_PATH);
	par.open(BARR_SRC_PATH, 'w');
	par.writeBarrel(b);
}

void AnaliticSolver::solve()
{
	cout << "\t<������� �������������� �������>\n";

	// ������ ������ � ������ �� �����-���������
	Parser par(BARR_SRC_PATH, 'r');
	Barrel barr(par.readBarrel());

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

	Parser::createFile(B_DELTA_PATH);
	par.open(B_DELTA_PATH, 'w');

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

	bool with_interp = false;
	char ch;
	cout << "\t\t������������ ���������� ������������ ��� ������� ��������� ����������? (+/-): ";
	cin >> ch;
	if (ch == '+') with_interp = true;

	Parser::createFile(Z_SLUH_PATH);
	Parser parsluh(Z_SLUH_PATH, 'w');
	Parser::createFile("results/LD.txt");
	Parser parLD("results/LD.txt", 'w');
	Parser::createFile("results/Ik.txt");
	Parser parIk("results/Ik.txt", 'w');
	Parser::createFile("results/W0.txt");
	Parser parW0("results/W0.txt", 'w');
	Parser::createFile("results/W.txt");
	Parser parW("results/W.txt", 'w');
	Parser::createFile("results/w_q.txt");
	Parser parWQ("results/w_q.txt", 'w');
	// ���������� ������� ������ �����
	for (vector<double>::iterator itr = eta_K.begin(); itr != eta_K.end(); itr++)
	{
		parsluh.write('\t', *itr);
		parLD.write('\t', *itr);
		parW0.write('\t', *itr);
		parW.write('\t', *itr);
		parIk.write('\t', *itr);
		parWQ.write('\t', *itr);
	}

	// ��������� �����
	for (vector<double>::iterator itr1 = Delta.begin(); itr1 != Delta.end(); itr1++)
	{
		barr.Delta = *itr1;
		// ����� �-�� �������� B* � ������������� ����� �������� ������� Lambda_K
		barr.calcB(a, b);
		barr.calcLambdaK();

		// �������� � 22 �� 34
		parsluh.write('\n', *itr1);
		parLD.write('\n', *itr1);
		parW.write('\n', *itr1);
		parW0.write('\n', *itr1);
		parIk.write('\n', *itr1);
		parWQ.write('\n', *itr1);
		for (vector<double>::iterator itr2 = eta_K.begin(); itr2 != eta_K.end(); itr2++)
		{
			barr.eta_K = *itr2;
			barr.calcForEtaK(*itr2);

			if (barr.r_D >= barr.r_Dmin) {
				double Z;
				if (with_interp)
				{
					double hi_n = 1.0 / (1.0 / Chuev(barr.CE15).hi + 0.75 * barr.d / barr.L0);
					Z = barr.calcZSluh(bilinearInterp(hi_n, barr.Lambda_D, table));
				}
				else Z = barr.calcZSluh();

				barrs.push_back(barr);

				parsluh.write('\t', Z);
				parLD.write('\t', barr.LD / barr.d);
				parW0.write('\t', barr.W0 * 1e3);			// � ��^3
				parW.write('\t', (barr.W + barr.W0) * 1e3);			// � ��^3
				parIk.write('\t', barr.Ik * 1e-6);							// � ���*�
				parWQ.write('\t', barr.omega_q);
			}
			else
			{
				parsluh.write('\t', 0.0);
				parLD.write('\t', 0.0);
				parW0.write('\t', 0.0);
				parW.write('\t', 0.0);
				parIk.write('\t', 0.0);
				parWQ.write('\t', 0.0);
			}
		}

		par.write(barr.Delta, '\t');
		par.write(barr.B, '\n');
	}

	// ������ � ����
	writeBarrelsToFile();
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

void AnaliticSolver::writeBarrelsToFile()
{
	Parser::createFile(RESULTS_PATH);
	Parser par(RESULTS_PATH, 'w');

	const int nvals = 9;
	Matrix mtrx(barrs.size(), nvals);
	Matrix m_sluh(21, 17 * 21);
	m_sluh.zeros();

	int i = 0;
	for (Barrels::const_iterator itr = barrs.begin(); itr != barrs.end(); itr++)
	{
		int j = 0;
		mtrx.data[i][j++] = itr->Delta;
		mtrx.data[i][j++] = itr->eta_K;

		mtrx.data[i][j++] = itr->Lambda_D;
		mtrx.data[i][j++] = itr->L0 / (itr->hi * itr->d);
		mtrx.data[i][j++] = itr->LD / (itr->d);
		mtrx.data[i][j++] = itr->omega_q;
		mtrx.data[i][j++] = (itr->W + itr->W0) * 1e3;
		mtrx.data[i][j++] = itr->Ik * 1e-6;
		mtrx.data[i][j] = itr->Z_Sluh;

		i++;
	}

	for (int i = 0; i < mtrx.x; i++)
	{
		for (int j = 0; j < mtrx.y - 1; j++)
			par.write(mtrx.data[i][j], '\t');
		par.write(mtrx.data[i][mtrx.y - 1], '\n');
	}
}

void DirectSolver::makeTest(const TestParams &tp)
{
	Parser::createFile(DIRSOL_TEST_PATH);
	Parser par(DIRSOL_TEST_PATH, 'w');

	double dt = tp.dt;
	pwd = tp.pwd;
	cout << pwd << endl;

	Result res;
	res.byDefault();
	res.Delta = tp.Delta;
	res.w_q = tp.w_q;						// ��������� �����������
	res.update(d, q, S, K, pwd.delta);

	calcToPmax(dt, res);
	par.write("p_max = ", res.p_max * 1e-6);

	continueCalc(dt, res);
	par.write("Vd = ", res.V);
	par.write("Ld = ", res.L);
	par.write("W0 = ", res.W0 * 1e3);
	par.write("W_ch = ", res.W_ch * 1e3);

	cout << "\t���������� ��. � " << DIRSOL_TEST_PATH << endl;
}

void DirectSolver::solve()
{
	cout << "\t<������� ������� ��������>\n";
	set_l_d_max();
	fillDelta();
	double dt = setTimeStep();

	showPowders();
	while (true)
	{
		results.clear();

		int indx = choosePowder();
		if (indx < 1 || indx > pwds.size())	// ������� ������ �� ����� -
			break;														// ���������� ���������� ������

		calcPmLine(dt, indx);
		fill_wq();
		calcIndicatDiag(dt, indx);
	}
	writeLmFile();
}

void DirectSolver::set_l_d_max()
{
	cout << "\t����������� �� ����� ������� �����, � ��������: ";
	cin >> l_d_max;
	if (l_d_max < 1.0 || l_d_max > 60.0)
	{
		status = "failed";
		throw err.sendMess("Error: max(l / d) must be > 1 and < 60!");
	}
}

void DirectSolver::fillDelta()
{
	cout << "\t������� ������������ ��������� ���������:\n";
	double from, to, step;
	cout << "\t - ��\t"; cin >> from;
	cout << "\t - ��\t"; cin >> to;
	do {
		cout << "\t - ���\t"; cin >> step;
	} while (step < eps);
	if (from > to)
	{
		double buf = from;
		from = to;
		to = buf;
	}

	size_d = unsigned((to - from) / step) + 1;
	Delta = new double[size_d];
	for (unsigned i = 0; i < size_d; i++)
		Delta[i] = from + i * step;
}

double DirectSolver::setTimeStep()
{
	double dt;
	cout << "\t��� �� �������, ���: ";
	cin >> dt;
	if (dt < 1.0 || dt > 5e4)
	{
		status = "failed";
		throw err.sendMess("Error: time step must be > 1.0 or < 5e+4 mcsec!");
	}

	return dt * 1e-6;
}

int DirectSolver::showPowders()
{
	int i = 0;
	cout << "\t������ �������:\n";
	for (Powders::const_iterator itr = pwds.begin(); itr != pwds.end(); itr++)
		cout << ++i << ") " << itr->name << endl;
	return i;
}

int DirectSolver::choosePowder()
{
	unsigned i;
	cout << "\t - ����� ������ (�� ������): ";
	cin >> i;
	if (i == 0 || i > pwds.size())
		return -1;
	i--;

	pwd = pwds[i];
	cout << "\t�������� �����: " << pwd.name << '\n';

	return i + 1;
}

void DirectSolver::calcPmLine(double dt, unsigned indx)
{
	// �������������� ������ ���� ��� ������ �����������
	string path;
	setPath(PM_PATH, path, indx);
	createFile(path, "Delta\tw/q\tW0\tW\tp");

	// ������� ������� ���
	cout << "\n\t omega / q:\n";
	for (unsigned i = 0; i < size_d; i++)
	{
		Result res;
		searchPmaxConds(dt, Delta[i], res);
		continueCalc(dt, res);

		cout << "\t\t " << res.w_q << '\n';
		writeFilePm(path, res);
	}
	cout << "\t���������� ��. � " << path << ".\n";
}

void DirectSolver::searchPmaxConds(double dt, double delta, Result &res)
{
	res.Delta = delta;
	res.w_q = 1.5;						// ��������� �����������

	double buf_w_q = 0.0;
	while (true)
	{
		res.byDefault();
		
		calcToPmax(dt, res);
		if (fabs(res.w_q - buf_w_q) < 1e-5)
			break;

		buf_w_q = res.w_q;
		res.w_q *= sqrt(pm / res.p);
		res.update(d, q, S, K, pwd.delta);
	}
}

void DirectSolver::continueCalc(double dt, Result &res)
{
	double buf_V = res.V;
	while (res.V < Vd)
	{
		rksolve(dt, res);
		if (res.V - buf_V < 1e-5 && res.V > 0.0)	// ���� ��� ������-���� ������
			break;																	// ���������� ������� �������� �������� Vd
		buf_V = res.V;
	}
	res.W_ch = res.W0 + S * res.L;
}

void DirectSolver::writeFilePm(const string &path, const Result &res)
{
	Parser par(path, 'a');

	par.write("\n", res.Delta, '\t');
	par.write(res.w_q, '\t');
	par.write(res.W0 * 1e3, '\t');
	par.write(res.W_ch * 1e3, '\t');
	par.write(res.p_max * 1e-6);
}

void DirectSolver::fill_wq()
{
	delete[] w_q;

	cout << "\t������� ������������ ������������� ����� ������:\n";
	double from, to;
	int num;
	cout << "\t - ��\t\t"; cin >> from;
	cout << "\t - ��\t\t"; cin >> to;
	cout << "\t - ���������\t"; cin >> num;

	if (from > to)
	{
		double buf = from;
		from = to;
		to = buf;
	}
	if (num < 1)
	{
		status = "failed";
		throw err.sendMess("Error: number of partitions must be > 0!");
	}

	size_wq = num + 1;
	double step = (to - from) / num;
	w_q = new double[size_wq];
	for (unsigned i = 0; i < size_wq; i++)
		w_q[i] = from + i * step;
}

void DirectSolver::calcIndicatDiag(double dt, unsigned indx)
{
	string path;
	setPath(I_DIAG_PATH, path, indx);
	createFile(path, "t\tDelta\tw/q\tV\tL/d\tW0\tW_ch\tpsi\tz\tp_max");

	for (unsigned i = 0; i < size_d; i++)
		for (unsigned j = 0; j < size_wq; j++)
		{
			Result res;
			// ��������� �������
			res.byDefault();
			res.Delta = Delta[i];
			res.w_q = w_q[j];
			res.update(d, q, S, K, pwd.delta);

			calcToPmax(dt, res);
			continueCalc(dt, res);
			results.push_back(res);

			writeFileDiag(path, res);
		}
	cout << "\t���������� ��. � " << path << ".\n";
}

void DirectSolver::calcToPmax(double dt, Result &res)
{
	Result buf;			// ����� ������� ���������� max(p)
	buf.p = 0;
	while (buf.p < res.p)
	{
		buf = res;
		rksolve(dt, res);
	}
	res = buf;
	res.p_max = res.p;
}

void DirectSolver::rksolve(double dt, Result &res)
{
	double fz[4], fpsi[4], fL[4], fV[4], fW[4], fp[4];

	fz[0] = dt * dz(res.z, res.p);
	fpsi[0] = dt * dpsi(res.z, res.p);
	fW[0] = dt * dW(q * res.w_q, res.z, res.p, res.V);
	fp[0] = dt * dp(q * res.w_q, res.F0, res.W, res.p, res.z, res.L, res.V);
	fV[0] = dt * dV(res.fi, res.V, res.p);
	fL[0] = dt * dL(res.V);

	fz[1] = dt * dz(res.z + 0.5 * fz[0], res.p + 0.5 * fp[0]);
	fpsi[1] = dt * dpsi(res.z + 0.5 * fz[0], res.p + 0.5 * fp[0]);
	fW[1] = dt * dW(q * res.w_q, res.z + 0.5 * fz[0], res.p + 0.5 * fp[0], res.V + 0.5 * fV[0]);
	fp[1] = dt * dp(q * res.w_q, res.F0,
		res.W + 0.5 * fW[0], res.p + 0.5 * fp[0], res.z + 0.5 * fz[0], res.L + 0.5 * fL[0], res.V + 0.5 * fV[0]);
	fV[1] = dt * dV(res.fi, res.V + 0.5 * fV[0], res.p + 0.5 * fp[0]);
	fL[1] = dt * dL(res.V + 0.5 * fV[0]);

	fz[2] = dt * dz(res.z + 0.5 * fz[1], res.p + 0.5 * fp[1]);
	fpsi[2] = dt * dpsi(res.z + 0.5 * fz[1], res.p + 0.5 * fp[1]);
	fW[2] = dt * dW(q * res.w_q, res.z + 0.5 * fz[1], res.p + 0.5 * fp[1], res.V + 0.5 * fV[1]);
	fp[2] = dt * dp(q * res.w_q, res.F0,
		res.W + 0.5 * fW[1], res.p + 0.5 * fp[1], res.z + 0.5 * fz[1], res.L + 0.5 * fL[1], res.V + 0.5 * fV[1]);
	fV[2] = dt * dV(res.fi, res.V + 0.5 * fV[1], res.p + 0.5 * fp[1]);
	fL[2] = dt * dL(res.V + 0.5 * fV[1]);

	fz[3] = dt * dz(res.z + fz[2], res.p + fp[2]);
	fpsi[3] = dt * dpsi(res.z + fz[2], res.p + fp[2]);
	fW[3] = dt * dW(q * res.w_q, res.z + fz[2], res.p + fp[2], res.V + fV[2]);
	fp[3] = dt * dp(q * res.w_q, res.F0,
		res.W + fW[2], res.p + fp[2], res.z + fz[2], res.L + fL[2], res.V + fV[2]);
	fV[3] = dt * dV(res.fi, res.V + fV[2], res.p + fp[2]);
	fL[3] = dt * dL(res.V + fV[2]);

	res.t += dt;
	res.z += (fz[0] + 2.0 * fz[1] + 2.0 * fz[2] + fz[3]) / 6.0;
	res.psi += (fpsi[0] + 2.0 * fpsi[1] + 2.0 * fpsi[2] + fpsi[3]) / 6.0;
	res.L += (fL[0] + 2.0 * fL[1] + 2.0 * fL[2] + fL[3]) / 6.0;
	res.V += (fV[0] + 2.0 * fV[1] + 2.0 * fV[2] + fV[3]) / 6.0;
	res.W += (fW[0] + 2.0 * fW[1] + 2.0 * fW[2] + fW[3]) / 6.0;
	res.p += (fp[0] + 2.0 * fp[1] + 2.0 * fp[2] + fp[3]) / 6.0;
}

void DirectSolver::writeFileDiag(const string &path, const Result &res)
{
	Parser par(path, 'a');

	par.write("\n", res.t, '\t');
	par.write(res.Delta, '\t');
	par.write(res.w_q, '\t');
	par.write(res.V, '\t');
	par.write(res.L / d, '\t');
	par.write(res.W0 * 1e3, '\t');
	par.write(res.W_ch * 1e3, '\t');
	par.write(res.psi, '\t');
	par.write(res.z, '\t');
	par.write(res.p_max * 1e-6);
}

void DirectSolver::writeResultsToFile(const std::string &path, const Results &rs)
{
	Parser::createFile(path);
	Parser par(path, 'w');

	par.write(pwd.name + "\n");
	par.write("t\tp\tV\tL\tpsi\tz\n");
	for (Results::const_iterator itr = rs.begin(); itr != rs.end(); itr++)
	{
		par.write(itr->t, '\t');
		par.write(itr->p, '\t');
		par.write(itr->V, '\t');
		par.write(itr->L / d, '\t');
		par.write(itr->psi, '\t');
		par.write(itr->z, '\n');
	}
}

void DirectSolver::writeLmFile(unsigned indx)
{
	string path;
	setPath(LM_PATH, path, indx);
	createFile(path, "W0\tW");

	double w[2];
	w[0] = 41e-3;
	w[1] = 46e-3;

	Parser par(path, 'a');
	par.write(funcW0(w[0]) * 1e3, '\t');
	par.write(w[0] * 1e3, '\n');
	par.write(funcW0(w[1]) * 1e3, '\t');
	par.write(w[1] * 1e3);
}

double DirectSolver::funcW0(double w)
{
	return w - l_d_max * d * S;
}

void DirectSolver::fillCriterionData(CriterionParams &cp)
{
	cp.d = d;
	cp.pm_star = pm;
	cp.l_d_max = l_d_max;

	cout << "\t - (l / d)_ref: ";
	cin >> cp.l_d_ref;
	cout << "\t - (W0 / d^3)_ref: ";
	cin >> cp.W0_d_ref;
	cout << "\t - (w / q)_ref: ";
	cin >> cp.w_q_ref;

	cout << "\n\t - alpha1: ";
	cin >> cp.alpha[0];
	cout << "\t - alpha2: ";
	cin >> cp.alpha[1];
	cout << "\t - alpha3: ";
	cin >> cp.alpha[2];
	cout << "\t - alpha4: ";
	cin >> cp.alpha[3];

	cout << "\n\t - ����-� �������� �-��� ����� ������� ����� � �������� (a): ";
	cin >> cp.a;
	cout << "\t - ����-� �������� �-��� �������� (b): ";
	cin >> cp.b;
}

void DirectSolver::calcCriterions()
{
	cout << "\n\t������ �������� �����������:\n";
	set_l_d_max();
	// ���������� ������ ��� ������� ��������
	CriterionParams cp;
	fillCriterionData(cp);

	showPowders();
	Criterions max_crs;
	while (true)
	{
		int indx = choosePowder();
		if (indx < 1 || indx > pwds.size())	// ������� ������ �� ����� -
			break;														// ���������� ���������� ������
		cp.pwd_name = pwd.name;

		string path;
		setPath(I_DIAG_PATH, path, indx);
		CResults rs;
		readResults(path, rs);

		calcCriterionCoeffs(rs, cp);
		Criterions crs;
		fillCriterions(rs, cp, crs);

		max_crs.push_back(*maxCriterion(crs.begin(), crs.end()));

		// ������ �������� � ����
		setPath(Z_PATH, path, indx);
		createFile(path);
		writeCriterionsFile(path, crs);
	}

	string path;
	setPath(Z_MAX_PATH, path);
	createFile(path, "", false);
	writeMaxCriterionFile(path, max_crs);
}

void DirectSolver::readResults(const string &path, CResults &rs)
{
	Parser par(path, 'r');
	for (int i = 0; i < 11; i++)		// ������ ���������
		par.readStr();

	while (true)
	{
		CResult res;
		res.t = par.readNext();
		res.Delta = par.readNext();
		res.w_q = par.readNext();
		res.V = par.readNext();
		res.L = par.readNext() * d;
		res.W0 = par.readNext() * 1e-3;
		res.W_ch = par.readNext() * 1e-3;
		res.psi = par.readNext();
		res.z = par.readNext();
		res.p_max = par.readNext() * 1e6;
		rs.push_back(res);

		if (par.isEnd())
			break;
	}
}

void DirectSolver::calcCriterionCoeffs(const CResults &rs, CriterionParams &cp)
{
	// ���������� � double ������� ��� �������� ������ �����������
	vector<double> v_ld;
	vector<double> v_W0;
	vector<double> v_wq;
	vector<double> v_pm;

	if (!rs.empty())
	{
		for (CResults::const_iterator itr = rs.begin(); itr != rs.end(); itr++)
		{
			v_ld.push_back(itr->L);
			v_W0.push_back(itr->W0);
			v_wq.push_back(itr->w_q);
			v_pm.push_back(itr->p_max);
		}
		cp.k[0] = 1;
		cp.k[1] = 1;
		cp.k[2] = 1;
		cp.k[3] = 1;
	}
}

Criterions::iterator DirectSolver::maxCriterion(const Criterions::iterator &start,
																								const Criterions::iterator &end)
{
	Criterions::iterator ans = start;
	double max = start->Z;
	for (Criterions::iterator itr = start; itr < end; itr++)
		if (max < itr->Z)
		{
			max = itr->Z;
			ans = itr;
		}
	return ans;
}

void DirectSolver::writeCriterionsFile(const std::string &path, const Criterions &crs)
{
	// ���������� ������ � ���� � ���� �������

	unsigned step = 0;
	for (Criterions::const_iterator i = crs.begin(); i != crs.end() - 1; i++)
	{
		step++;
		if (fabs((i + 1)->Delta - i->Delta) > eps)
			break;
	}

	Parser par(path, 'a');
	par.write('\n', 0);
	for (unsigned i = 0; i < crs.size(); i += step)
		par.write('\t', crs[i].Delta);

	for (unsigned i = 0; i < step; i++)
	{
		par.write('\n', crs[i].w_q);
		for (unsigned j = i; j < crs.size(); j += step)
			par.write('\t', crs[j].Z);
	}
}

void DirectSolver::writeMaxCriterionFile(const std::string &path, const Criterions &crs)
{
	Parser par(path, 'a');
	for (Criterions::const_iterator i = crs.begin(); i != crs.end(); i++)
	{
		par.write("\n" + i->pwd_name + "\t");
		par.write(i->Delta, '\t');
		par.write(i->w_q, '\t');
		par.write(i->Z);
	}
}

void DirectSolver::solveOnce()
{
	cout << "\n\t<������� ������ ������>\n";

	double delta, wq;
	char ch;
	cout << "\t������ ��� max(Z)? (+/-): ";
	cin >> ch;

	double dt = setTimeStep();

	if (ch == '+')
	{
		string path;
		setPath(Z_MAX_PATH, path);
		Criterion max_cr;
		getMaxCriterion(path, max_cr);
		cout << "\t������������ ��������: " << max_cr << endl;

		delta = max_cr.Delta;
		wq = max_cr.w_q;
	}
	else
	{
		cout << "\t - ��������� ���������, ��/�^3:\t";
		cin >> delta;
		cout << "\t - ������������� ����� ������:\t";
		cin >> wq;
	}

	showPowders();
	int indx = choosePowder();
	Result res;
	res.byDefault();
	res.Delta = delta;
	res.w_q = wq;
	res.update(d, q, S, K, pwd.delta);
	Results rs;
	calcToPmax(dt, res, rs);
	continueCalc(dt, res, rs);

	string path;
	setPath(SOL_ONCE_PATH, path, indx);
	createFile(path, "t\tp\tV\tL\tpsi\tz");
	writeResultFile(path, rs);
}

void DirectSolver::calcToPmax(double dt, Result &res, Results &rs)
{
	Result buf;			// ����� ������� ���������� max(p)
	buf.p = 0;
	while (buf.p < res.p)
	{
		rs.push_back(res);
		buf = res;
		rksolve(dt, res);
	}
	res = buf;
	res.p_max = res.p;
}

void DirectSolver::continueCalc(double dt, Result &res, Results &rs)
{
	rs.push_back(res);
	double buf_V = res.V;
	while (res.V < Vd)
	{
		rksolve(dt, res);
		rs.push_back(res);
		if (res.V - buf_V < 1e-5 && res.V > 0.0)	// ���� ��� ������-���� ������
			break;																	// ���������� ������� �������� �������� Vd
		buf_V = res.V;
	}
	res.W_ch = res.W0 + S * res.L;
}

void DirectSolver::getMaxCriterion(const string &path, Criterion &max_cr)
{
	Parser par(path, 'r');
	Criterions crs;
	while (!par.isEnd())
	{
		Criterion cr;
		cr.pwd_name = par.readStr();
		cr.Delta = par.readNext();
		cr.w_q = par.readNext();
		cr.Z = par.readNext();
		crs.push_back(cr);
	}
	max_cr = *maxCriterion(crs.begin(), crs.end());
}

void DirectSolver::writeResultFile(const string &path, const Results &rs)
{
	Parser par(path, 'a');
	for (Results::const_iterator itr = rs.begin(); itr != rs.end(); itr++)
	{
		par.write("\n", itr->t, '\t');
		par.write(itr->p * 1e-6, '\t');
		par.write(itr->V, '\t');
		par.write(itr->L, '\t');
		par.write(itr->psi, '\t');
		par.write(itr->z);
	}
}

void DirectSolver::fillCriterions(const CResults &rs, const CriterionParams &cp, Criterions &crs)
{
	char ch;
	cout << "�������� ����������/���������� ��������? (+/-): ";
	cin >> ch;

	if (ch == '+')
		for (CResults::const_iterator itr = rs.begin(); itr != rs.end(); itr++)
		{
			Criterion cr;
			cr.Delta = itr->Delta;
			cr.w_q = itr->w_q;

			double L0 = itr->W0 / S;
			double hi = 1.666;
			cr.calcCriterionSluh(L0, itr->L, hi, d, itr->p_max / pm, 1.5);
			crs.push_back(cr);
		}
	else
		for (CResults::const_iterator itr = rs.begin(); itr != rs.end(); itr++)
		{
			Criterion cr;
			cr.calcCriterion(*itr, cp);
			crs.push_back(cr);
		}
}

void DirectSolver::setPath(const string &base_path, string &res_path, unsigned num)
{
	res_path = base_path;
	if (num != 0)
	{
		char buf[3];
		_itoa_s(num, buf, 10);

		res_path = base_path;
		res_path += buf;
		res_path += ".txt";
		return;
	}
	res_path += ".txt";
}

void DirectSolver::createFile(const string &path, const string &head, bool w_name)
{
	Parser::createFile(path);
	if (!head.empty() || w_name)
	{
		Parser par(path, 'w');
		if (w_name)
			par.write(pwd.name);
		if (!head.empty())
			par.write("\n" + head);
	}
}

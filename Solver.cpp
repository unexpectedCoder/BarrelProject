#include "Solver.h"

#include <cstdlib>
#include <math.h>

using namespace std;

int AnaliticSolver::file_count = 0;

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

	printResults();
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
	// ��� ������������ ������
	/*cout << "\nTable:\n";
	for (int i = 0; i < ROWS; i++)
	{
	for (int j = 0; j < COLOUMS; j++)
	cout << table[i][j] << '\t';
	cout << '\n';
	}*/

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

void DirectSolver::solve()
{
	cout << "\t<������� ������� ��������>\n";
	fillData();

	// ������� ���� �� �������
	double dt;
	do
	{
		cout << "\t��� �� �������, ���: "; cin >> dt;
	} while (dt < 1.0 && dt > 1e4);
	dt *= 1e-6;

	char ch;
	////////////////////////////// ����� pm ///////////////////////////////
	cout << "\t����� ������ ��� ����� pm �� ������������ ���������? (+/-): ";
	cin >> ch;
	if (ch == '+')
	{
		showPowders();
		while (true)
		{
			int indx = choosePowder();
			if (indx == -1) break;

			// �������� ������ ����������� ��� ��������� �������
			char buf[3];
			_itoa_s(indx, buf, 10);
			string path = "direct_res/pm_"; path += buf; path += ".txt";
			Parser::createFile(path);
			Parser par(path, 'w');

			////////////////////////// ������� ������� ��� ////////////////////////////
			par.write(pwd.name + "\n");
			par.write("Delta\tw/q\tW0\tW\tp\n");
			for (int i = 0; i < size_d; i++)
			{
				Results rs;
				for (int j = 0; j < size_wq; j++)
					searchPmax(dt, Delta[i], w_q[j], rs);

				Result res = *minDeltaPm(rs);				// ������ ������������ � ������ ������� (��. ����)
				// ������ ����������� � ����
				par.write(res.Delta, '\t');
				par.write(res.w_q, '\t');
				par.write(res.W0 * 1e3, '\t');
				par.write(res.W * 1e3, '\t');
				par.write(res.p * 1e-6, '\n');

				rs.~vector();
			}

			cout << "\t���������� ��. � " << path << ".\n";
		}
	}
	
	////////////////////////////// ������ ������� ///////////////////////////////
	cout << "\t<��������� ������� ����>\n";
	showPowders();
	while (true)
	{
		int indx = choosePowder();
		if (indx == -1) break;

		char buf[3];
		_itoa_s(indx, buf, 10);
		string path = "direct_res/IDiag_"; path += buf; path += ".txt";
		Parser::createFile(path);
		Parser par(path, 'w');
		par.write(pwd.name + "\n");
		par.write("W0\tW\n");
		
		////////////////////////// ������� ������� ��� ////////////////////////////
		for (int i = 0; i < size_d; i++)
			for (int j = 0; j < size_wq; j++)
			{
				Results rs;

				// ������� �������� �� ���������� Vd
				sys(dt, Delta[i], w_q[j], rs);
				par.write((rs.end() - 1)->W0 * 1e3, '\t');
				par.write((rs.end() - 1)->W * 1e3, '\n');

				// �������� ������ ����������� ��� ��������� �������
				/*
				char buf1[3], buf2[3], buf3[3];
				_itoa_s(indx, buf1, 10);
				_itoa_s(i + 1, buf2, 10);
				_itoa_s(j + 1, buf3, 10);
				path = "direct_res/res_"; path += buf1; path += "_"; path += buf2; path += "_"; path += buf3; path += ".txt";
				// ������ � ����
				writeResultsToFile(path, rs);
				*/

				rs.~vector();
			}

		cout << "\t���������� ��. � " << path << ".\n";
	}
}

void DirectSolver::fillData()
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

	size_d = (to - from) / step + 1;
	Delta = new double[size_d];
	for (int i = 0; i < size_d; i++)
		Delta[i] = from + i * step;

	cout << "\t������� ������ ������������� ����� ������:\n";
	cout << "\t - ��\t\t"; cin >> from;
	cout << "\t - ��\t\t"; cin >> to;
	do {
		cout << "\t - ���������\t"; cin >> size_wq;
	} while (size_wq <= 0);
	if (from > to)
	{
		double buf = from;
		from = to;
		to = buf;
	}

	step = (to - from) / size_wq;
	size_wq++;
	w_q = new double[size_wq];
	for (int i = 0; i < size_wq + 1; i++)
		w_q[i] = from + i * step;
}

int DirectSolver::choosePowder()
{
	int i;
	cout << "\t - ����� ������ (�� ������): ";  cin >> i;
	i--;
	if (i < 0 && i >= pwds.size())
		return -1;

	pwd = pwds[i];
	cout << "\t�������� �����: " << pwd.name << '\n';

	return i + 1;
}

void DirectSolver::searchPmax(double dt, double delta, double wq, Results &rs)
{
	key_V = 0; key_S = 1; key_Z = 1;

	Result res;
	res.Delta = delta;
	res.w_q = wq;
	res.W0 = wq / delta;
	// ���. �������
	res.L = 0.0;
	res.p = Consts::p_flash;
	res.psi = 0.0;
	res.t = 0.0;
	res.V = 0.0;
	res.z = 0.0;
	res.W = q * wq / delta - q * wq / pwd.delta;

	// ��������������� �������
	double buf_p = 0.0;
	double fi = K + 1.0 / 1.3 * wq;
	double W0 = q * wq / delta;
	double F0 = 4.0 * W0 / d + 2.0 * S;

	// ������� �������
	double fz[4], fpsi[4], fL[4], fV[4], fW[4], fp[4];
	while (res.p > buf_p)
	{
		buf_p = res.p;

		fz[0] = dz(res.p);
		fpsi[0] = dpsi(res.z, res.p);
		fL[0] = dL(res.V);
		fV[0] = dV(res.p, wq, fi);
		fW[0] = dW(q * wq, res.z, res.p, res.V);
		fp[0] = dp(res.W, q * wq, delta, res.p, res.z, res.L, res.V, F0);

		fz[1] = dz(res.p + 0.5 * fp[0] * dt);
		fpsi[1] = dpsi(res.z + 0.5 * fz[0] * dt, res.p + 0.5 * fp[0] * dt);
		fL[1] = dL(res.V + 0.5 * fV[0] * dt);
		fV[1] = dV(res.p + 0.5 * fp[0] * dt, wq, fi);
		fW[1] = dW(q * wq, res.z + 0.5 * fz[0] * dt, res.p + 0.5 * fp[0] * dt, res.V + 0.5 * fV[0] * dt);
		fp[1] = dp(res.W + 0.5 * fW[0] * dt, q * wq, delta, res.p + 0.5 * fp[0] * dt,
			res.z + 0.5 * fz[0] * dt, res.L + 0.5 * fL[0] * dt, res.V + 0.5 * fV[0] * dt, F0);

		fz[2] = dz(res.p + 0.5 * fp[1] * dt);
		fpsi[2] = dpsi(res.z + 0.5 * fz[1] * dt, res.p + 0.5 * fp[1] * dt);
		fL[2] = dL(res.V + 0.5 * fV[1] * dt);
		fV[2] = dV(res.p + 0.5 * fp[1] * dt, wq, fi);
		fW[2] = dW(q * wq, res.z + 0.5 * fz[1] * dt, res.p + 0.5 * fp[1] * dt, res.V + 0.5 * fV[1] * dt);
		fp[2] = dp(res.W + 0.5 * fW[1] * dt, q * wq, delta, res.p + 0.5 * fp[1] * dt,
			res.z + 0.5 * fz[1] * dt, res.L + 0.5 * fL[1] * dt, res.V + 0.5 * fV[1] * dt, F0);

		fz[3] = dz(res.p + fp[2] * dt);
		fpsi[3] = dpsi(res.z + fz[2] * dt, res.p + fp[2] * dt);
		fL[3] = dL(res.V + fV[2] * dt);
		fV[3] = dV(res.p + fp[2] * dt, wq, fi);
		fW[3] = dW(q * wq, res.z + fz[2] * dt, res.p + fp[2] * dt, res.V + fV[2] * dt);
		fp[3] = dp(res.W + fW[2] * dt, q * wq, delta, res.p + fp[2] * dt,
			res.z + fz[2] * dt, res.L + fL[2] * dt, res.V + fV[2] * dt, F0);

		res.t += dt;
		res.z += (fz[0] + 2.0 * fz[1] + 2.0 * fz[2] + fz[3]) * dt / 6.0;
		res.psi += (fpsi[0] + 2.0 * fpsi[1] + 2.0 * fpsi[2] + fpsi[3]) * dt / 6.0;
		res.L += (fL[0] + 2.0 * fL[1] + 2.0 * fL[2] + fL[3]) * dt / 6.0;
		res.V += (fV[0] + 2.0 * fV[1] + 2.0 * fV[2] + fV[3]) * dt / 6.0;
		res.W += (fW[0] + 2.0 * fW[1] + 2.0 * fW[2] + fW[3]) * dt / 6.0;
		res.p += (fp[0] + 2.0 * fp[1] + 2.0 * fp[2] + fp[3]) * dt / 6.0;

		if (res.p > Consts::p_flash) key_V = 1;
		if (res.z > 1.0) key_S = 0;
		if (res.z > pwd.zk) key_Z = 0;
	}

	rs.push_back(res);
}

void DirectSolver::sys(double dt, double delta, double wq, Results &rs)
{
	key_V = 0; key_S = 1; key_Z = 1;

	Result res;
	res.Delta = delta;
	res.w_q = wq;
	res.W0 = wq / delta;
	// ���. �������
	res.L = 0.0;
	res.p = Consts::p_flash;
	res.psi = 0.0;
	res.t = 0.0;
	res.V = 0.0;
	res.z = 0.0;
	res.W = q * wq / delta - q * wq / pwd.delta;

	// ��������������� �������
	double buf_p = 0.0;
	double fi = K + 1.0 / 1.3 * wq;
	double W0 = q * wq / delta;
	double F0 = 4.0 * W0 / d + 2.0 * S;

	rs.push_back(res);
	double fz[4], fpsi[4], fL[4], fV[4], fW[4], fp[4];
	while (res.V < Vd)
	{
		fz[0] = dz(res.p);
		fpsi[0] = dpsi(res.z, res.p);
		fL[0] = dL(res.V);
		fV[0] = dV(res.p, wq, fi);
		fW[0] = dW(q * wq, res.z, res.p, res.V);
		fp[0] = dp(res.W, q * wq, delta, res.p, res.z, res.L, res.V, F0);

		fz[1] = dz(res.p + 0.5 * fp[0] * dt);
		fpsi[1] = dpsi(res.z + 0.5 * fz[0] * dt, res.p + 0.5 * fp[0] * dt);
		fL[1] = dL(res.V + 0.5 * fV[0] * dt);
		fV[1] = dV(res.p + 0.5 * fp[0] * dt, wq, fi);
		fW[1] = dW(q * wq, res.z + 0.5 * fz[0] * dt, res.p + 0.5 * fp[0] * dt, res.V + 0.5 * fV[0] * dt);
		fp[1] = dp(res.W + 0.5 * fW[0] * dt, q * wq, delta, res.p + 0.5 * fp[0] * dt,
			res.z + 0.5 * fz[0] * dt, res.L + 0.5 * fL[0] * dt, res.V + 0.5 * fV[0] * dt, F0);

		fz[2] = dz(res.p + 0.5 * fp[1] * dt);
		fpsi[2] = dpsi(res.z + 0.5 * fz[1] * dt, res.p + 0.5 * fp[1] * dt);
		fL[2] = dL(res.V + 0.5 * fV[1] * dt);
		fV[2] = dV(res.p + 0.5 * fp[1] * dt, wq, fi);
		fW[2] = dW(q * wq, res.z + 0.5 * fz[1] * dt, res.p + 0.5 * fp[1] * dt, res.V + 0.5 * fV[1] * dt);
		fp[2] = dp(res.W + 0.5 * fW[1] * dt, q * wq, delta, res.p + 0.5 * fp[1] * dt,
			res.z + 0.5 * fz[1] * dt, res.L + 0.5 * fL[1] * dt, res.V + 0.5 * fV[1] * dt, F0);

		fz[3] = dz(res.p + fp[2] * dt);
		fpsi[3] = dpsi(res.z + fz[2] * dt, res.p + fp[2] * dt);
		fL[3] = dL(res.V + fV[2] * dt);
		fV[3] = dV(res.p + fp[2] * dt, wq, fi);
		fW[3] = dW(q * wq, res.z + fz[2] * dt, res.p + fp[2] * dt, res.V + fV[2] * dt);
		fp[3] = dp(res.W + fW[2] * dt, q * wq, delta, res.p + fp[2] * dt,
			res.z + fz[2] * dt, res.L + fL[2] * dt, res.V + fV[2] * dt, F0);

		res.t += dt;
		res.z += (fz[0] + 2.0 * fz[1] + 2.0 * fz[2] + fz[3]) * dt / 6.0;
		res.psi += (fpsi[0] + 2.0 * fpsi[1] + 2.0 * fpsi[2] + fpsi[3]) * dt / 6.0;
		res.L += (fL[0] + 2.0 * fL[1] + 2.0 * fL[2] + fL[3]) * dt / 6.0;
		res.V += (fV[0] + 2.0 * fV[1] + 2.0 * fV[2] + fV[3]) * dt / 6.0;
		res.W += (fW[0] + 2.0 * fW[1] + 2.0 * fW[2] + fW[3]) * dt / 6.0;
		res.p += (fp[0] + 2.0 * fp[1] + 2.0 * fp[2] + fp[3]) * dt / 6.0;

		if (res.p > Consts::p_flash) key_V = 1;
		if (res.z > 1.0) key_S = 0;
		if (res.z > pwd.zk) key_Z = 0;

		rs.push_back(res);
	}
}

Results::const_iterator DirectSolver::minDeltaPm(const Results &rs)
{
	Results::const_iterator itr = rs.begin();
	double min = fabs(rs.begin()->p - pm);

	for (Results::const_iterator i = rs.begin(); i != rs.end(); i++)
		if (fabs(i->p - pm) < min && i->p < pm)
		{
			min = fabs(i->p - pm);
			itr = i;
		}

	return itr;
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

int DirectSolver::showPowders()
{
	int i = 0;
	cout << "\t������ �������:\n";
	for (Powders::const_iterator itr = pwds.begin(); itr != pwds.end(); itr++)
		cout << ++i << ") " << *itr << endl;

	return i;
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

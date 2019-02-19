#include "AnaliticSolver.h"

using namespace std;


AnaliticSolver::~AnaliticSolver()
{
	barrs.~vector();
	Delta.~vector();
	eta_K.~vector();
}

void AnaliticSolver::printIntro()
{
	cout << "\n\t<Аналитический решатель (AnaliticSolver)>\n";
}

void AnaliticSolver::printOutro()
{
	cout << "STATUS: " << status << ".\n";
	if (status != "failed")
		cout << "RESULTS: look in \'results\'.\n\n";
}

void AnaliticSolver::calcMaxPressure()
{
	cout << "\t<Функция поиска уровня макс. давления>\n";

	Barrel b;
	char ch;
	b.CE15 = b.CE;
	b.eta_omega = Chuev(b.CE15).eta_omega;
	do
	{
		b.CE15 = calcCE15(b.Cq, b.CE, b.eta_omega);
		b.eta_omega = Chuev(b.CE15).eta_omega;				// eta_omega recalculation for new CE15
		b.pm_kr = Chuev(b.CE15).pm_kr;
		b.omega_q = b.CE15 / (b.Cq * b.eta_omega);
		b.fi = b.K + 1.0 / 3.0 * b.omega_q;
		b.pm = b.pm_kr * b.fi * Consts::Nkr / (Consts::fi1 + 0.5 * b.omega_q) * NEW_BARR_K;

		cout << "\t\tСейчас " << "CE15 = " << b.CE15 << ", pm = " << b.pm << endl;
		cout << "\t\tУточнить решение? (+/-): ";
		cin >> ch;
	}
	while (ch == '+');

	double pm_nround = b.pm;	// Stores not rounded up pm
	b.pm *= Consts::g * 1e4;	// Transfer to Pa from atm

	cout << "\n\t\tДля вашего образца: pm = " << b.pm * 1e-6 << " МПа\n";
	cout << "\t\t - введите округленное значение pm, МПа: pm = ";
	cin >> b.pm;
	b.pm *= 1e6;							// Transfer to Pa from MPa

	b.p_mid = 0.5 * b.pm;
	b.calcHi1();
	b.hi = Chuev(b.CE15).hi;

	writeFileFromMaxPress(b, pm_nround);
}

void AnaliticSolver::writeFileFromMaxPress(const Barrel &b, double pm_nround)
{
	Parser::createFileTXT(BARR_LOG_PATH);
	Parser par(BARR_LOG_PATH, 'w');
	par.write("p_mid = ", b.p_mid, '\n');
	par.write("hi1 = ", b.hi1, '\n');
	// Write to the file as a table
	makeTableTxt(b, pm_nround, BARR_TABLE_PATH);
	// Add to the p_CE15.txt file
	par.open(P_CE15_PATH, 'a');
	par.write(b.CE15, '\t');
	par.write(b.pm / Consts::g * 1e-4, '\n');
	// Writing data about the barrel to a file for further reading
	Parser::createFileTXT(BARR_SRC_PATH);
	par.open(BARR_SRC_PATH, 'w');
	par.writeBarrel(b);
}

void AnaliticSolver::solve()
{
	cout << "\t<Функция аналитического решения>\n";

	Parser par(BARR_SRC_PATH, 'r');
	Barrel barr(par.readBarrel());
	fillBarrelData("Задайте плотность заряжания (Delta), кг/м^3:", Delta);
	fillBarrelData("Задайте безразмерную координату конца горения (eta_K):", eta_K);
	// Reading from table #3.3 (Chuev)
	par.open(CHUEV_TABLE_33_PATH, 'r');
	double **table = new double*[ROWS];
	for (int i = 0; i < ROWS; ++i)
	{
		table[i] = new double[COLOUMS];
		for (int j = 0; j < COLOUMS; ++j)
			table[i][j] = par.readDouble();
	}

	Parser::createFileTXT(B_DELTA_PATH);
	par.open(B_DELTA_PATH, 'w');
	// To search for the Drozdov parameter B*
	double a, b;
	cout << "\t\tЗадайте границы поиска параметра Дроздова B*:\n";
	cout << "\t\t - слева: "; cin >> a;
	cout << "\t\t - справа: "; cin >> b;
	if (a > b)
	{
		double buf = a;
		a = b;
		b = buf;
	}

	bool with_interp = false;
	char ch;
	cout << "\t\tИспользовать билинейную интерполяцию при расчете параметра Слухоцкого? (+/-): ";
	cin >> ch;
	if (ch == '+')
		with_interp = true;

	createNecessaryFiles();
	Parser parsluh(Z_SLUH_PATH, 'w');
	Parser parLD(LD_PATH, 'w');
	Parser parIk(IK_PATH, 'w');
	Parser parW0(W0_PATH, 'w');
	Parser parW(W_PATH, 'w');
	Parser parWQ(W_Q_PATH, 'w');
	// Writing head of file
	for (vector<double>::iterator itr = eta_K.begin(); itr != eta_K.end(); itr++)
	{
		parsluh.write('\t', *itr);
		parLD.write('\t', *itr);
		parW0.write('\t', *itr);
		parW.write('\t', *itr);
		parIk.write('\t', *itr);
		parWQ.write('\t', *itr);
	}
	// Solution
	for (vector<double>::iterator itr1 = Delta.begin(); itr1 != Delta.end(); itr1++)
	{
		barr.Delta = *itr1;
		// Search for the Drozdov coefficient B * and the relative length of the combustion process Lambda_K
		barr.calcB(a, b);
		barr.calcLambdaK();
		// Actions 22 to 34 (look in lectures)
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
				parW0.write('\t', barr.W0 * 1e3);						// dm^3
				parW.write('\t', (barr.W + barr.W0) * 1e3);	// dm^3
				parIk.write('\t', barr.Ik * 1e-6);					// MPa*s
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
	writeBarrelsToFile();
}

void AnaliticSolver::fillBarrelData(const string &head_txt, vector<double> &data)
{
	if (!data.empty()) data.clear();

	double start, end, step;
	cout << "\t\t" << head_txt << '\n';
	cout << "\t\t - от: "; cin >> start;
	cout << "\t\t - до: "; cin >> end;
	cout << "\t\t - шаг: "; cin >> step;

	while (start < end + 0.5 * step)
	{
		data.push_back(start);
		start += step;
	}
}

void AnaliticSolver::makeTableTxt(const Barrel &barr, double pm_nround, const string &path)
{
	Parser::createFileTXT(path);
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
	Parser::createFileTXT(RESULTS_PATH);
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
		++i;
	}

	for (int i = 0; i < mtrx.x; ++i)
	{
		for (int j = 0; j < mtrx.y - 1; ++j)
			par.write(mtrx.data[i][j], '\t');
		par.write(mtrx.data[i][mtrx.y - 1], '\n');
	}
}

void AnaliticSolver::createNecessaryFiles()
{
	Parser::createFileTXT(Z_SLUH_PATH);
	Parser::createFileTXT(LD_PATH);
	Parser::createFileTXT(IK_PATH);
	Parser::createFileTXT(W0_PATH);
	Parser::createFileTXT(W_PATH);
	Parser::createFileTXT(W_Q_PATH);
}

#include "DirectSolver.h"

using namespace std;


DirectSolver::DirectSolver(double _pm, double _d, double _q, double _Vd, double _K, double _p0, double _ns) :
	status("successfully"), pm(_pm), d(_d), q(_q), Vd(_Vd), K(_K), p0(_p0), ns(_ns)
{
	pwds = Parser().readXMLPowders(POWDERS_PATH);
	S = 0.25 * pi * pow(d, 2.0) * ns;
	Ld = 0;
}

DirectSolver::~DirectSolver()
{
	delete[] Delta;
	delete[] w_q;
	pwds.~vector();
}

void DirectSolver::makeTest(const TestParams &tp)
{
	Parser::createFileTXT(DIRSOL_TEST_PATH);
	Parser par(DIRSOL_TEST_PATH, 'w');

	double dt = tp.dt;
	pwd = tp.pwd;
	cout << pwd << endl;

	Result res;
	res.byDefault();
	res.Delta = tp.Delta;
	res.w_q = tp.w_q;				// Initial approximation
	res.update(d, q, S, K, pwd.delta);

	calcToPmax(dt, res);
	par.write("p_max = ", res.p_max * 1e-6);

	continueCalc(dt, res);
	par.write("Vd = ", res.V);
	par.write("Ld = ", res.L);
	par.write("W0 = ", res.W0 * 1e3);
	par.write("W_ch = ", res.W_ch * 1e3);

	cout << "\tРезультаты см. в " << DIRSOL_TEST_PATH << endl;
}

void DirectSolver::solve()
{
	cout << "\t<Функция прямого решателя>\n";
	set_l_d_max();
	fillDelta();
	double dt = setTimeStep();

	showPowders();
	while (true)
	{
		results.clear();
		int indx = choosePowder();
		if (indx < 1 || indx > pwds.size())	// Condition out of the loop -
			break;														// no powder selected

		calcPmLine(dt, indx);
		fill_wq();
		calcIndicatDiag(dt, indx);
	}
	writeLmFile();
}

void DirectSolver::set_l_d_max()
{
	cout << "\tОграничение по длине ведущей части, в калибрах: ";
	cin >> l_d_max;
	if (l_d_max < 1.0 || l_d_max > 60.0)
	{
		status = "failed";
		throw "Error: max(l / d) must be > 1 and < 60!";
	}
}

void DirectSolver::fillDelta()
{
	cout << "\tПределы варьирования плотности заряжания:\n";
	double from, to, step;
	cout << "\t - от\t"; cin >> from;
	cout << "\t - до\t"; cin >> to;
	do
	{
		cout << "\t - шаг\t";
		cin >> step;
	}
	while (step < eps);
	if (from > to)
	{
		double buf = from;
		from = to;
		to = buf;
	}

	size_d = unsigned((to - from) / step) + 1;
	Delta = new double[size_d];
	for (unsigned i = 0; i < size_d; ++i)
		Delta[i] = from + i * step;
}

double DirectSolver::setTimeStep()
{
	double dt;
	cout << "\tШаг по времени, мкс: ";
	cin >> dt;
	if (dt < 1.0 || dt > 5e4)
	{
		status = "failed";
		throw "Error: time step must be > 1.0 or < 5e+4 mcsec!";
	}
	return dt * 1e-6;
}

int DirectSolver::showPowders()
{
	int i = 0;
	cout << "\tСписок порохов:\n";
	for (Powders::const_iterator itr = pwds.begin(); itr != pwds.end(); ++itr)
		cout << ++i << ") " << itr->name << endl;
	return i;
}

int DirectSolver::choosePowder()
{
	unsigned i;
	cout << "\t - номер пороха (из списка): ";
	cin >> i;
	if (i == 0 || i > pwds.size())
		return -1;
	i--;

	pwd = pwds[i];
	cout << "\tВыбраный порох: " << pwd.name << '\n';

	return i + 1;
}

void DirectSolver::calcPmLine(double dt, unsigned indx)
{
	// Initializes a special way to record results.
	string path;
	setPath(PM_PATH, path, indx);
	createFile(path, "Delta\tw/q\tW0\tW\tp");
	// ODE system solution
	cout << "\n\t omega / q:\n";
	for (unsigned i = 0; i < size_d; ++i)
	{
		Result res;
		searchPmaxConds(dt, Delta[i], res);
		continueCalc(dt, res);

		cout << "\t\t " << res.w_q << '\n';
		writeFilePm(path, res);
	}
	cout << "\tРезультаты см. в " << path << ".\n";
}

void DirectSolver::searchPmaxConds(double dt, double delta, Result &res)
{
	res.Delta = delta;
	res.w_q = 1.5;				// Initial approximation

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
		if (res.V - buf_V < 1e-5 && res.V > 0.0)	// If for any powder it is impossible
			break;																	// to achieve a given speed Vd
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
	par.write(res.p_max * 1e-6, '\n');
}

void DirectSolver::fill_wq()
{
	delete[] w_q;

	cout << "\tПределы варьирования относительной массы заряда:\n";
	double from, to;
	int num;
	cout << "\t - от\t\t"; cin >> from;
	cout << "\t - до\t\t"; cin >> to;
	cout << "\t - разбиений\t"; cin >> num;
	if (from > to)
	{
		double buf = from;
		from = to;
		to = buf;
	}
	if (num < 1)
	{
		status = "failed";
		throw "Error: number of partitions must be > 0!";
	}

	size_wq = num + 1;
	double step = (to - from) / num;
	w_q = new double[size_wq];
	for (unsigned i = 0; i < size_wq; ++i)
		w_q[i] = from + i * step;
}

void DirectSolver::calcIndicatDiag(double dt, unsigned indx)
{
	string path;
	setPath(I_DIAG_PATH, path, indx);
	createFile(path, "t\tDelta\tw/q\tV\tL/d\tW0\tW_ch\tpsi\tz\tp_max");

	for (unsigned i = 0; i < size_d; ++i)
		for (unsigned j = 0; j < size_wq; ++j)
		{
			// Initial approximation
			Result res;
			res.byDefault();
			res.Delta = Delta[i];
			res.w_q = w_q[j];
			res.update(d, q, S, K, pwd.delta);

			calcToPmax(dt, res);
			continueCalc(dt, res);
			results.push_back(res);

			writeFileDiag(path, res);
		}
	cout << "\tРезультаты см. в " << path << ".\n";
}

void DirectSolver::calcToPmax(double dt, Result &res)
{
	Result buf;			// To catch the achievement max(p)
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
	Parser::createFileTXT(path);
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
	string path_W, path_D;
	setPath(LM_W_PATH, path_W, indx);
	setPath(LM_DELTA_PATH, path_D, indx);
	createFile(path_W, "W0\tW");
	createFile(path_D, "Delta\tw / q");

	double W[2];
	W[0] = 41e-3;
	W[1] = 46e-3;

	Parser par(path_W, 'a');
	par.write("\n", funcW0(W[0]) * 1e3, '\t');
	par.write(W[0] * 1e3, '\n');
	par.write(funcW0(W[1]) * 1e3, '\t');
	par.write(W[1] * 1e3);
	par.close();

	double delta[2];
	delta[0] = 760;
	delta[1] = 790;
	double wq[2];
	wq[0] = funcW0(W[0]) * delta[0] / q;
	wq[1] = funcW0(W[1]) * delta[1] / q;

	par.open(path_D, 'a');
	par.write("\n", delta[0], '\t');
	par.write(wq[0], '\n');
	par.write(delta[1], '\t');
	par.write(wq[1]);
}

double DirectSolver::funcW0(double W)
{
	return W - l_d_max * d * S;
}

void DirectSolver::calcCriterions()
{
	cout << "\n\tРасчет критерия оптимизации:\n";
	set_l_d_max();
	// Fill in the data to calculate the criterion
	CriterionParams cp;
	fillCriterionData(cp);

	showPowders();
	Criterions max_crs;
	while (true)
	{
		int indx = choosePowder();
		if (indx < 1 || indx > pwds.size())	// Condition out of the loop is
			break;														// no powder selected

		cp.pwd_name = pwd.name;
		string path;
		setPath(I_DIAG_PATH, path, indx);
		CResults rs;
		readResults(path, rs);
		calcCriterionCoeffs(rs, cp);
		Criterions crs;
		fillCriterions(rs, cp, crs);
		max_crs.push_back(*maxCriterion(crs.begin(), crs.end()));

		// Write criteria to file
		setPath(Z_PATH, path, indx);
		createFile(path);
		writeCriterionsFile(path, crs);
	}
	string path;
	setPath(Z_MAX_PATH, path);
	createFile(path, "", false);
	writeMaxCriterionFile(path, max_crs);
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

	cout << "\n\t - коэф-т штрафной ф-ции длины ведущей части в калибрах (a): ";
	cin >> cp.a;
	cout << "\t - коэф-т штрафной ф-ции давления (b): ";
	cin >> cp.b;
}

void DirectSolver::readResults(const string &path, CResults &rs)
{
	// Read title
	Parser par(path, 'r');
	for (int i = 0; i < 11; ++i)
		par.readStr();

	while (true)
	{
		CResult res;
		res.t = par.readDouble();
		res.Delta = par.readDouble();
		res.w_q = par.readDouble();
		res.V = par.readDouble();
		res.L = par.readDouble() * d;
		res.W0 = par.readDouble() * 1e-3;
		res.W_ch = par.readDouble() * 1e-3;
		res.psi = par.readDouble();
		res.z = par.readDouble();
		res.p_max = par.readDouble() * 1e6;
		rs.push_back(res);

		if (par.isEnd())
			break;
	}
}

void DirectSolver::calcCriterionCoeffs(const CResults &rs, CriterionParams &cp)
{
	// Overwriting in double vector for easy search of extrema
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
	// The function writes data to the file as a matrix.
	unsigned step = 0;
	for (Criterions::const_iterator i = crs.begin(); i != crs.end() - 1; ++i)
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
	for (Criterions::const_iterator i = crs.begin(); i != crs.end(); ++i)
	{
		par.write("\n" + i->pwd_name + "\t");
		par.write(i->Delta, '\t');
		par.write(i->w_q, '\t');
		par.write(i->Z);
	}
}

void DirectSolver::solveOnce(double dt, int start_temp_C)
{
	cout << "\n\t<Решение прямой задачи>\n";
	if (start_temp_C == T_N)
	{
		solveStandartTemp(dt);
		return;
	}
	if (Ld == 0)
	{
		cout << "Warning: требуется решение при стандартной температуре для определения длины ведущей части!\n";
		return;
	}

	double delta, wq;
	string path;
	setPath(Z_MAX_PATH, path);
	Criterion max_cr;
	getMaxCriterion(path, max_cr);
	cout << "\tМаксимальный критерий: " << max_cr << endl;
	delta = max_cr.Delta;
	wq = max_cr.w_q;

	showPowders();
	choosePowder();
	pwd.abnormalTemperature(start_temp_C);

	Result res;
	res.byDefault();
	res.Delta = delta;
	res.w_q = wq;
	res.update(d, q, S, K, pwd.delta);
	Results rs;
	calcToPmax(dt, res, rs);
	continueCalcLd(dt, res, rs);

	setPath(SOL_ONCE_PATH, path, start_temp_C);
	createFile(path, "t\tp\tV\tL\tpsi\tz");
	writeResultFile(path, rs);
}

void DirectSolver::solveStandartTemp(double dt)
{
	double delta, wq;
	string path;
	setPath(Z_MAX_PATH, path);
	Criterion max_cr;
	getMaxCriterion(path, max_cr);
	cout << "\tМаксимальный критерий: " << max_cr << endl;
	delta = max_cr.Delta;
	wq = max_cr.w_q;

	showPowders();
	choosePowder();
	// pwd.abnormalTemperature(start_temp_C);

	Result res;
	res.byDefault();
	res.Delta = delta;
	res.w_q = wq;
	res.update(d, q, S, K, pwd.delta);
	Results rs;
	calcToPmax(dt, res, rs);
	continueCalcVd(dt, res, rs);

	Ld = res.L;

	setPath(SOL_ONCE_PATH, path, T_N);
	createFile(path, "t\tp\tV\tL\tpsi\tz");
	writeResultFile(path, rs);
}

void DirectSolver::calcToPmax(double dt, Result &res, Results &rs)
{
	Result buf;			// To catch the achievement max(p)
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

void DirectSolver::continueCalcVd(double dt, Result &res, Results &rs)
{
	rs.push_back(res);
	double buf_V = res.V;
	while (res.V < Vd)
	{
		rksolve(dt, res);
		rs.push_back(res);
		if (res.V - buf_V < 1e-5 && res.V > 0.0)	// If for any powder it is impossible
			break;																	// to achieve a given speed Vd
		buf_V = res.V;
	}
	res.W_ch = res.W0 + S * res.L;
}

void DirectSolver::continueCalcLd(double dt, Result &res, Results &rs)
{
	rs.push_back(res);
	while (res.L < Ld)
	{
		rksolve(dt, res);
		rs.push_back(res);
	}
	rs.pop_back();
	res.W_ch = res.W0 + S * res.L;
}

void DirectSolver::getMaxCriterion(const string &path, Criterion &max_cr)
{
	Parser par(path, 'r');
	Criterions crs;
	readCriterionFile(path, crs);
	max_cr = *maxCriterion(crs.begin(), crs.end());
}

void DirectSolver::readCriterionFile(const string &path, Criterions &crs)
{
	Parser par(path, 'r');
	while (!par.isEnd())
	{
		Criterion cr;
		cr.pwd_name = par.readStr();
		cr.Delta = par.readDouble();
		cr.w_q = par.readDouble();
		cr.Z = par.readDouble();
		crs.push_back(cr);
	}
}

void DirectSolver::fillCriterions(const CResults &rs, const CriterionParams &cp, Criterions &crs)
{
	for (CResults::const_iterator itr = rs.begin(); itr != rs.end(); ++itr)
	{
		Criterion cr;
		cr.calcCriterion(*itr, cp);
		crs.push_back(cr);
	}
}

void DirectSolver::writeResultFile(const string &path, const Results &rs)
{
	Parser par(path, 'a');
	for (Results::const_iterator itr = rs.begin(); itr != rs.end(); ++itr)
	{
		par.write("\n", itr->t, '\t');
		par.write(itr->p * 1e-6, '\t');
		par.write(itr->V, '\t');
		par.write(itr->L, '\t');
		par.write(itr->psi, '\t');
		par.write(itr->z);
	}
}

void DirectSolver::setPath(const string &base_path, string &res_path, int num)
{
	res_path = base_path;
	if (num != 0)
	{
		char buf[3];
		_itoa_s(fabs(num), buf, 10);

		res_path = base_path;
		if (num < 0)
			res_path += "_";
		res_path += buf;
		res_path += ".txt";
		return;
	}
	res_path += ".txt";
}

void DirectSolver::createFile(const string &path, const string &head, bool w_name)
{
	Parser::createFileTXT(path);
	if (!head.empty() || w_name)
	{
		Parser par(path, 'w');
		if (w_name)
			par.write(pwd.name);
		if (!head.empty())
			par.write("\n" + head);
	}
}


#include "AnalogsSolver.h"

using namespace std;


void AnalogsSolver::printIntro()
{
	cout << "\n\t<AnalogsSolver>\n";
}

void AnalogsSolver::printOutro() {
	cout << "STATUS: " << status << ".\n";
	if (status != "failed")
		cout << "RESULTS: look in the file " << ANALOGS_PATH << ".\n\n";
}

void AnalogsSolver::solve()
{
	cout << "\t<Функция расчета аналогов>\n";

	char ch;
	cout << "\tВнести данные об аналогах? (+/-):";
	cin >> ch;
	if (ch == '+')
		fillAnalogsData();

	Parser par(ANALOGS_PATH, 'r');
	double K;
	cout << "\t - коэф-т учета второстепенных работ: K = ";
	cin >> K;
	analogs.clear();
	while (true)
	{
		string name = par.readStr();
		double d = par.readDouble() * 1e-3;
		double q = par.readDouble();
		double vd = par.readDouble();
		double l_d = par.readDouble();
		double pm = par.readDouble();

		if (par.isEnd())	// Exit loop. With this arrangement,
			break;					// there will be no work with variables in which the last (empty) line is written,
											// due to which the values will be inadequate.

		Analog a(name, d, q, vd, l_d, pm);
		a.CE15 = a.CE;
		a.eta_omega = Chuev(a.CE15).eta_omega;
		do
		{
			a.CE15 = calcCE15(a.Cq, a.CE, a.eta_omega);
			a.eta_omega = Chuev(a.CE15).eta_omega;

			cout << "\tОрудие " << a.name << ":\n";
			cout << "\t\t - уточнить решение? (+/-): ";
			cin >> ch;
		}
		while (ch == '+');

		analogs.push_back(a);
		cout << '\n';
	}

	writeSolveFile(par);
	printOutro();
}

void AnalogsSolver::fillAnalogsData()
{
	cout << "\t<Функция заполнения данных об аналогах>\n";

	char ch;
	cout << "\t\tДобавить аналог? (+/-): ";
	cin >> ch;
	int num = 0;
	while (ch == '+')
	{
		string name;
		double d, q, vd, l_d, pm;
		cout << "\t\tАналог №" << ++num << ":\n";
		cout << "\t\t - название: "; cin >> name;
		cout << "\t\t - d, мм: "; cin >> d;
		cout << "\t\t - q, кг: "; cin >> q;
		cout << "\t\t - Vd, м/с: "; cin >> vd;
		cout << "\t\t - длина ведущей части в калибрах: "; cin >> l_d;
		cout << "\t\t - максимальное давление, кг/см^2: "; cin >> pm;
		analogs.push_back(Analog(name, d, q, vd, l_d, pm));

		cout << "\t\tДобавить аналог? (+/-): "; cin >> ch;
	}
	writeAnalogsFile();
}

void AnalogsSolver::writeAnalogsFile()
{
	Parser::createFileTXT(ANALOGS_PATH);
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

void AnalogsSolver::writeSolveFile(Parser &par)
{
	Parser::createFileTXT(ANALOGS_RES_PATH);
	Parser parres(ANALOGS_RES_PATH, 'w');
	Parser::createFileTXT(P_CE15_PATH);
	par.open(P_CE15_PATH, 'w');
	if (!par.isOpen() || !parres.isOpen())
	{
		status = "faild";
		return;
	}

	for (Analogs::iterator itr = analogs.begin(); itr != analogs.end(); ++itr)
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
}

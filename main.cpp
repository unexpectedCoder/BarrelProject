#include <iostream>
#include <Windows.h>

#include "Solver.h"
#include "Types.h"

using namespace std;

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	cout << "\t\t___*** ВНУТРЕННЯЯ БАЛЛИСТИКА ***___\n\n\tМГТУ\n\t2018\n\n";

	char choice;

	cout << "Произвести тестовое решение? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		TestSolver test;
		test.printIntro();
		test.solve();
		test.printOutro();
	}

	// Расчет аналогов
	cout << "Рассчитать аналоги? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		AnalogsSolver asol;
		asol.printIntro();
		asol.solve();
		asol.printOutro();
	}

	cout << "Решить ОЗВБ аналитически? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		// Расчет pm собственного образца
		AnaliticSolver analitic;
		analitic.printIntro();

		cout << "\tРассчитать макс. давление собственного образца? (+/-): ";
		cin >> choice;
		if (choice == '+')
			analitic.calcMaxPressure();

		analitic.solve();
		analitic.printOutro();
	}

	DirectSolver dirsol;
	cout << "\n\nРешить прямую задачу перебором? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		dirsol.printIntro();

		cout << "\tПроизвести тестовое решение? (+/-): ";
		cin >> choice;
		if (choice == '+')
		{
			DirectSolver test(0, 0.122, 21.76, 690, 1.05, 30e6, 1.04);
			Powder pwd;
			pwd.f = 1.004e6;
			pwd.k = 1.226;
			pwd.alpha = 1.013e-3;
			pwd.T = 2890;
			pwd.delta = 1600;
			pwd.Ik = 0.81e6;
			pwd.zk = 1.522;
			pwd.kappa1 = 0.725;
			pwd.lambda1 = 0.183;
			pwd.kappa2 = 0.546;
			pwd.lambda2 = -0.957;
			pwd.kappa_f = 0.0003;
			pwd.k_f = 0.0016;

			test.makeTest(TestParams(1e-5, 500, 0.18, pwd));
		}

		try {
			dirsol.solve();
		}
		catch (const string &ex) {
			cout << ex << endl;
		}
		dirsol.printOutro();
	}

	cout << "Рассчитать критерий оптимизации? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		dirsol.calcCriterions();
		dirsol.printOutro();
	}

	cout << "Решить прямую задачу? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		dirsol.solveOnce();
		dirsol.printOutro();
	}

	cout << endl << endl;
	system("pause");
	return 0;
}

#include <iostream>

#include "Solver.h"

using namespace std;

int main()
{
	setlocale(LC_ALL, "RUS");

	cout << "\t\t___*** ВНУТРЕННЯЯ БАЛЛИСТИКА ***___\n\n\tМГТУ\n\t2018\n\n";

	char choice;

	cout << "Произвести тестовое решение? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		TestSolver *test = new TestSolver();
		test->printInfo();
		test->solve();
		test->printResults();
	}

	// Расчет аналогов
	cout << "Рассчитать аналоги? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		AnalogsSolver *asol = new AnalogsSolver;
		asol->printInfo();
		asol->solve();
		asol->printResults();
	}

	cout << "Решить ОЗВБ аналитически? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		// Расчет pm собственного образца
		AnaliticSolver *analitic = new AnaliticSolver;
		analitic->printInfo();

		cout << "\tРассчитать макс. давление собственного образца? (+/-): ";
		cin >> choice;
		if (choice == '+')
			analitic->calcMaxPressure();

		analitic->solve();
		analitic->printResults();
	}

	cout << endl;
	system("pause");
	return 0;
}

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
		AnaliticSolver *analitic = new AnaliticSolver();
		analitic->printInfo();

		cout << "\tРассчитать макс. давление собственного образца? (+/-): ";
		cin >> choice;
		if (choice == '+')
			analitic->calcMaxPressure();

		analitic->solve();
		analitic->printResults();
	}

	cout << "\n\nРешить прямую задачу перебором? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		DirectSolver *dirsol = new DirectSolver();
		dirsol->printInfo();
		dirsol->solve();
		dirsol->printResults();
	}

	cout << endl << endl;
	system("pause");
	return 0;
}

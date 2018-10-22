#include <iostream>

#include "Solver.h"

using namespace std;

int main()
{
	setlocale(LC_ALL, "RUS");

	cout << "\t\t___*** ВНУТРЕННЯЯ БАЛЛИСТИКА ***___\n\n\tМГТУ\n\t2018\n\n";

	char choice;

	Solver test;
	cout << "Произвести тестовое решение? (+/-): ";
	cin >> choice;
	if (choice == '+')
		test.makeTest();

	Solver sol;

	// Заполнение аналогов
	cout << "Заполнить данные об аналогах? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.fillAnalogs();

	// Расчет аналогов
	cout << "Рассчитать макс. давление для аналогов? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.calcAnalogs();

	// Расчет pm собственного образца
	cout << "Рассчитать макс. давление собственного образца? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.calcBarrelPressure();

	// Аналитическое решение ОЗВБ
	cout << "Решить ОЗВБ (аналитически)? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.solveInvProblem();

	cout << endl;
	system("pause");
	return 0;
}

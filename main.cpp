#include <iostream>

#include "Solver.h"

using namespace std;

int main()
{
	setlocale(LC_ALL, "RUS");

	cout << "\t\t___*** ���������� ���������� ***___\n\n\t����\n\t2018\n\n";

	char choice;

	cout << "���������� �������� �������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		TestSolver *test = new TestSolver();
		test->printInfo();
		test->solve();
		test->printResults();
	}

	// ������ ��������
	cout << "���������� �������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		AnalogsSolver *asol = new AnalogsSolver;
		asol->printInfo();
		asol->solve();
		asol->printResults();
	}

	cout << "������ ���� ������������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		// ������ pm ������������ �������
		AnaliticSolver *analitic = new AnaliticSolver;
		analitic->printInfo();

		cout << "\t���������� ����. �������� ������������ �������? (+/-): ";
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

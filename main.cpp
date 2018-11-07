#include <iostream>
#include <Windows.h>

#include "Solver.h"
#include "Types.h"

using namespace std;

int main()
{
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

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
		AnaliticSolver *analitic = new AnaliticSolver();
		analitic->printInfo();

		cout << "\t���������� ����. �������� ������������ �������? (+/-): ";
		cin >> choice;
		if (choice == '+')
			analitic->calcMaxPressure();

		analitic->solve();
		analitic->printResults();
	}

	cout << "\n\n������ ������ ������ ���������? (+/-): ";
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

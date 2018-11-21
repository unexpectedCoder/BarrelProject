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
		TestSolver test;
		test.printInfo();
		test.solve();
		test.printResults();
	}

	// ������ ��������
	cout << "���������� �������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		AnalogsSolver asol;
		asol.printInfo();
		asol.solve();
		asol.printResults();
	}

	cout << "������ ���� ������������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		// ������ pm ������������ �������
		AnaliticSolver analitic;
		analitic.printInfo();

		cout << "\t���������� ����. �������� ������������ �������? (+/-): ";
		cin >> choice;
		if (choice == '+')
			analitic.calcMaxPressure();

		analitic.solve();
		analitic.printResults();
	}

	cout << "\n\n������ ������ ������ ���������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		DirectSolver dirsol;
		dirsol.printInfo();

		cout << "\t���������� �������� �������? (+/-): ";
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
		dirsol.printResults();
	}

	cout << endl << endl;
	system("pause");
	return 0;
}

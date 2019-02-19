#include <iostream>
#include <Windows.h>

#include "AnaliticSolver.h"
#include "AnalogsSolver.h"
#include "DirectSolver.h"
#include "TestSolver.h"
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
		test.printIntro();
		test.solve();
		test.printOutro();
	}

	cout << "���������� �������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		AnalogsSolver asol;
		asol.printIntro();
		asol.solve();
		asol.printOutro();
	}

	cout << "������ ���� ������������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		// Calculate maximum of pressure own sample
		AnaliticSolver analitic;
		analitic.printIntro();

		cout << "\t���������� ����. �������� ������������ �������? (+/-): ";
		cin >> choice;
		if (choice == '+')
			analitic.calcMaxPressure();

		analitic.solve();
		analitic.printOutro();
	}

	DirectSolver dirsol;
	cout << "\n\n������ ������ ������ ���������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		dirsol.printIntro();

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
			pwd.k_f = 0.0003;
			pwd.k_I = 0.0016;

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

	cout << "���������� �������� �����������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		dirsol.calcCriterions();
		dirsol.printOutro();
	}

	cout << "������ ������ ������? (+/-): ";
	cin >> choice;
	if (choice == '+')
	{
		double dt = 0;
		cout << "\t��� �� �������, ���: ";
		while (true)
		{
			cin >> dt;
			dt *= 1e-6;
			if (dt > 1e-8 || dt < 1e-4)
				break;
			cout << "\t��? ����� ��-�����.\n";
		}

		dirsol.solveOnce(dt, 15);
		cout << "\t����� ������ ��� +15C, �: " << dirsol.getBarLen() << endl;

		dirsol.solveOnce(dt, -50);
		dirsol.solveOnce(dt, 50);

		dirsol.printOutro();
	}

	cout << endl << endl;
	system("pause");
	return 0;
}

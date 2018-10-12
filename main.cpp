#include <iostream>

#include "Solver.h"

using namespace std;

int main()
{
	cout << "\t\t___*** INTERNAL BALLISTICS ***___\n\n\tBMSTU\n\t2018\n\n";

	Solver sol;
	Chuev chuev;
	char choice;

	// ���������� ��������
	cout << "Fill analogs data? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.fillAnalogs();

	// ������ ��������
	cout << "Calculate analogs? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.calcAnalogs();

	// ������ pm ������������ �������
	cout << "Calculate pm of your own sample? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.calcBarrelPressure();

	// ������������� ������� ����
	cout << "Solve the inverse problem of internal ballistics? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.solveInvProblem();
		

	cout << endl;
	system("pause");
	return 0;
}

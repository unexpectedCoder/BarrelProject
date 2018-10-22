#include <iostream>

#include "Solver.h"

using namespace std;

int main()
{
	setlocale(LC_ALL, "RUS");

	cout << "\t\t___*** ���������� ���������� ***___\n\n\t����\n\t2018\n\n";

	char choice;

	Solver test;
	cout << "���������� �������� �������? (+/-): ";
	cin >> choice;
	if (choice == '+')
		test.makeTest();

	Solver sol;

	// ���������� ��������
	cout << "��������� ������ �� ��������? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.fillAnalogs();

	// ������ ��������
	cout << "���������� ����. �������� ��� ��������? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.calcAnalogs();

	// ������ pm ������������ �������
	cout << "���������� ����. �������� ������������ �������? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.calcBarrelPressure();

	// ������������� ������� ����
	cout << "������ ���� (������������)? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.solveInvProblem();

	cout << endl;
	system("pause");
	return 0;
}

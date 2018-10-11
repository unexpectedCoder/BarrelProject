#include <iostream>

#include "Solver.h"

using namespace std;

int main()
{
	Solver sol;
	Chuev chuev;
	char choice;

	cout << "Fill analogs data? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.fillAnalogs();

	cout << "Calculate analogs? (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.calcAnalogs();

	cout << "Calculate your own sample& (+/-): ";
	cin >> choice;
	if (choice == '+')
		sol.calcBarrelPressure();

	cout << endl;
	system("pause");
	return 0;
}

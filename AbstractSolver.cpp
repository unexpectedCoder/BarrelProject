#include "AbstractSolver.h"


double AbstractSolver::bilinearInterp(double x, double y, double **table)
{
	int ix, iy;

	/*
	* ������� table[0][0] �� ���������������, �.�. ����� ������ ��� �������� ������� (��. ���� chuev_table_33.txt.).
	* ��� �������� � ����, ��� � ������ i ���������� � 1, � �� � 0.
	*/
	// ����� �������� ��-��, ����� �������� ��������� �������� �����...
	// ...�� ��� COLOUMS (X)
	for (int i = 1; i < COLOUMS; i++)
		if (x - table[0][i] <= 1.0)
		{
			ix = i;
			break;
		}
	// ...�� ��� ROWS (Y)
	for (int i = 1; i < ROWS - 1; i++)
		if (y >= table[i][0] && y <= table[i + 1][0])
		{
			iy = i;
			break;
		}

	// ��������� �����
	double fQ11 = table[iy][ix];					// ���������������� �������� - ��������� ������������
	double fQ21 = table[iy][ix + 1];			// �������� ������ � ������:
	double fQ12 = table[iy + 1][ix];			// 2�� ������ (�������) ������������� ��� X;
	double fQ22 = table[iy + 1][ix + 1];	// 1�� ������ (������) - ��� Y.

	double fR1 = 1.0 / (table[0][ix + 1] - table[0][ix]) * ((table[0][ix + 1] - x) * fQ11 + (x - table[0][ix]) * fQ21);
	double fR2 = 1.0 / (table[0][ix + 1] - table[0][ix]) * ((table[0][ix + 1] - x) * fQ12 + (x - table[0][ix]) * fQ22);

	return 1.0 / (table[iy + 1][0] - table[iy][0]) * ((table[iy + 1][0] - y) * fR1 + (y - table[iy][0]) * fR2);
}

#include "AbstractSolver.h"


double AbstractSolver::bilinearInterp(double x, double y, double **table)
{
	int ix, iy;

	/*
	* Ёлемент table[0][0] не рассматриваетс€, т.к. нужен только дл€ создани€ матрицы (см. файл chuev_table_33.txt.).
	* Ёто приводит к тому, что в циклах i начинаетс€ с 1, а не с 0.
	*/
	// ѕоиск индексов эл-ов, между которыми находитс€ заданна€ точка...
	// ...по оси COLOUMS (X)
	for (int i = 1; i < COLOUMS; i++)
		if (x - table[0][i] <= 1.0)
		{
			ix = i;
			break;
		}
	// ...по оси ROWS (Y)
	for (int i = 1; i < ROWS - 1; i++)
		if (y >= table[i][0] && y <= table[i + 1][0])
		{
			iy = i;
			break;
		}

	// –асчетна€ часть
	double fQ11 = table[iy][ix];					// Ќетрадиционность индексов - следствие особенностей
	double fQ21 = table[iy][ix + 1];			// хранени€ матриц в пам€ти:
	double fQ12 = table[iy + 1][ix];			// 2ой индекс (столбцы) соответствует оси X;
	double fQ22 = table[iy + 1][ix + 1];	// 1ый индекс (строки) - оси Y.

	double fR1 = 1.0 / (table[0][ix + 1] - table[0][ix]) * ((table[0][ix + 1] - x) * fQ11 + (x - table[0][ix]) * fQ21);
	double fR2 = 1.0 / (table[0][ix + 1] - table[0][ix]) * ((table[0][ix + 1] - x) * fQ12 + (x - table[0][ix]) * fQ22);

	return 1.0 / (table[iy + 1][0] - table[iy][0]) * ((table[iy + 1][0] - y) * fR1 + (y - table[iy][0]) * fR2);
}

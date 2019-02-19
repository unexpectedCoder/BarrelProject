#ifndef ABSTRACTSOLVER_H
#define ABSTRACTSOLVER_H

#include <iostream>
#include <stack>
#include <cstring>
#include <math.h>

// For interpolation of Chuev
#define ROWS 18
#define COLOUMS 11

#define P_CE15_PATH "files/p_CE15.txt"


class AbstractSolver
{
public:
	AbstractSolver() {}
	virtual ~AbstractSolver() {}

	virtual void printIntro() = 0;
	virtual void solve() = 0;
	virtual void printOutro() = 0;

	virtual double calcCE15(double cq, double ce, double eta) {
		return 0.5 * 15.0 / cq * (ce - 3.0 * eta * cq + sqrt(pow(ce - 3.0 * eta * cq, 2.0) + 4.0 / 5.0 * ce * eta * pow(cq, 2.0)));
	}
	virtual double bilinearInterp(double hi_n, double Lambda_D, double **table);
};

#endif

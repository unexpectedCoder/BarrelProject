#ifndef TESTSOLVER_H
#define TESTSOLVER_H

/*
	* The class is designed to verify the calculations.
	* Results are compared with the “reference” program.
*/

#include "AbstractSolver.h"
#include "Types.h"

#define TEST_PATH "files/test.txt"


class TestSolver : public AbstractSolver
{
private:
	Barrel *barr;
	std::string status;

public:
	TestSolver();
	~TestSolver();

	void printIntro();
	void printOutro();
	void solve();
};

#endif

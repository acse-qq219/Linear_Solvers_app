#include <cmath>
#include <iostream>
#include <ctime>
#include "Matrix.h"
#include "Solver.h"
#include "CSRMatrix.h"
#include "Interface.h"
#include "Matrix.cpp"
#include "Solver.cpp"
#include "CSRMatrix.cpp"
#include "Interface.cpp"
#include "Test.h"
using namespace std;

int main() {

	auto* nice = new Interface<double>;

	delete nice;
}
#include "Test.h"
#include <math.h>
#include <ctime>
#include "Matrix.h"
#include "Matrix.cpp"
#include "Solver.h"
#include "Solver.cpp"
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <iostream>
#include <cassert>
using namespace std;
template<typename T>
bool IsAlmostEqual(T* x, T* y, int length, int ulp);
Test::Test()
{};
Test::~Test()
{};
void Test::main_test()
{
	cout << "       TEST BEGAIN" << endl;
	this->Solver_test();
	this->CSR_test();
	cout << endl << endl << "Thank you for your patient, You can start using" << endl << endl;

}
void Test::Jacob_test()
{
	int rows = 3, cols = 3;

	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	auto* dense_mat = new Matrix<double>(rows, cols, values);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* jaco = new Solver<double>();
	jaco->solverJacobi(100, 1e-5, 3, *dense_mat, b_vec, x_vec, counter);
	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  JACOB TEST SUCCESS!     ***" << endl;;

	}
	else
	{
		cout << "***  JACOB TEST FAIL!     *** " << endl;

	}


	delete jaco;
	delete[] x_vec;

};
void Test::GausSeid_test()
{

	int rows = 3, cols = 3;

	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	auto* dense_mat = new Matrix<double>(rows, cols, values);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* GS = new Solver<double>();
	GS->solverGausSeid(100, 1e-5, 3, *dense_mat, b_vec, x_vec, counter);
	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  GAUSSEID TEST SUCCESS!  ***" << endl;;

	}
	else
	{
		cout << "***  GAUSSEID TEST FAIL!  ***" << endl;

	}


	delete GS;
	delete[] x_vec;




};
void Test::GausElim_test()
{
	int rows = 3, cols = 3;

	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	auto* dense_mat = new Matrix<double>(rows, cols, values);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* GE = new Solver<double>();
	GE->solverGausElim(*dense_mat, b_vec, x_vec);
	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  GAUSELIM TEST SUCCESS!  ***" << endl;;

	}
	else
	{
		cout << "***  GAUSELIM TEST SUCCESS!  *** " << endl;

	}


	delete GE;
	delete[] x_vec;



};
void Test::LuDecom_test()
{
	int rows = 3, cols = 3;

	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	auto* dense_mat = new Matrix<double>(rows, cols, values);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* LU = new Solver<double>();
	LU->solverLuDecom(*dense_mat, b_vec, x_vec);
	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  LU TEST SUCCESS!        ***" << endl;;

	}
	else
	{
		cout << "***  LU TEST FAIL!           *** " << endl;

	}


	delete LU;
	delete[] x_vec;



};

void Test::SOR_test()
{
	int rows = 3, cols = 3;

	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	auto* dense_mat = new Matrix<double>(rows, cols, values);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* SOR = new Solver<double>();
	SOR->solverSor(200, 1e-5, 2, *dense_mat, b_vec, x_vec, counter);
	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  SOR TEST SUCCESS!       ***" << endl;;

	}
	else
	{
		cout << "***  SOR TEST FAIL!          *** " << endl;

	}


	delete SOR;
	delete[] x_vec;



}

void Test::Solver_test()
{
	cout << endl << endl << "        SOLVER TEST" << endl << endl;
	this->Jacob_test();
	this->GausSeid_test();
	this->GausElim_test();

	this->LuDecom_test();
	this->SOR_test();
}
void Test::CSR_test()
{
	cout << endl << endl << "        CSR TEST" << endl << endl;
	this->CSR_GET_DIAG_test();
	this->CSR_GET_REM_test();
	this->CSR_MATMUL_test();
	this->CSR_jacob_test();
	this->CSR_GausSeid_test();
	this->CSR_GausElim_test();
	this->CSR_LUDecom_test();
	this->CSR_SOR_test();

};
void Test::CSR_jacob_test()
{
	int rows = 3, cols = 3;
	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	int row_position[4] = { 0, 3, 6,9 };
	int col_index[9] = { 0,1,2,0,1,2,0,1,2 };
	auto* test_CSR = new CSRMatrix<double>(rows, cols, 9, values, row_position, col_index);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* jaco = new Solver<double>();
	jaco->solverJacobi(100, 1e-5, 3, *test_CSR, b_vec, x_vec, counter);
	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  CSR JACOB TEST SUCCESS!     ***" << endl;;

	}
	else
	{
		cout << "***  CSR JACOB TEST FAIL!        ***" << endl;

	}


	delete jaco;
	delete[] x_vec;



};
void Test::CSR_GausSeid_test()
{
	int rows = 3, cols = 3;
	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	int row_position[4] = { 0, 3, 6,9 };
	int col_index[9] = { 0,1,2,0,1,2,0,1,2 };
	auto* test_CSR = new CSRMatrix<double>(rows, cols, 9, values, row_position, col_index);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* GS = new Solver<double>();
	GS->solverGausSeid(100, 1e-5, 3, *test_CSR, b_vec, x_vec, counter);
	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  CSR GS TEST SUCCESS!        ***" << endl;;

	}
	else
	{
		cout << "***  CSR GS TEST FAIL!        *** " << endl;

	}


	delete GS;
	delete[] x_vec;

};
void Test::CSR_GausElim_test()
{
	int rows = 3, cols = 3;
	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	int row_position[4] = { 0, 3, 6,9 };
	int col_index[9] = { 0,1,2,0,1,2,0,1,2 };
	auto* test_CSR = new CSRMatrix<double>(rows, cols, 9, values, row_position, col_index);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* GE = new Solver<double>();
	GE->solverGausElim(*test_CSR, b_vec, x_vec);
	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  CSR GE TEST SUCCESS!        ***" << endl;;

	}
	else
	{
		cout << "***  CSR GE TEST FAIL!        *** " << endl;

	}


	delete GE;
	delete[] x_vec;


};
void Test::CSR_LUDecom_test()
{
	int rows = 3, cols = 3;
	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	int row_position[4] = { 0, 3, 6,9 };
	int col_index[9] = { 0,1,2,0,1,2,0,1,2 };
	auto* test_CSR = new CSRMatrix<double>(rows, cols, 9, values, row_position, col_index);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* LU = new Solver<double>();
	LU->solverLuDecom(*test_CSR, b_vec, x_vec);

	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  CSR LU TEST SUCCESS!        ***" << endl;;

	}
	else
	{
		cout << "***  CSR LU TEST FAIL!           *** " << endl;

	}


	delete LU;
	delete[] x_vec;



};
void Test::CSR_SOR_test()
{
	int rows = 3, cols = 3;
	double values[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	int row_position[4] = { 0, 3, 6,9 };
	int col_index[9] = { 0,1,2,0,1,2,0,1,2 };
	auto* test_CSR = new CSRMatrix<double>(rows, cols, 9, values, row_position, col_index);
	double b_vec[3] = { 24464, 23281, 491 };

	double* x_vec = new double[rows];
	for (int i = 0; i < rows; i++)
	{
		x_vec[i] = 0;

	}
	int counter = 0;
	auto* SOR = new Solver<double>();
	SOR->solverSor(200, 1e-5, 2, *test_CSR, b_vec, x_vec, counter);
	double x_result[3] = { 0.458748,0.350392,-0.335172 };

	if (IsAlmostEqual(x_result, x_vec, 3, 2))
	{
		cout << "***  CSR SOR TEST SUCCESS!       ***" << endl;;

	}
	else
	{
		cout << "***  CSR SOR TEST FAIL!          *** " << endl;

	}


	delete SOR;
	delete[] x_vec;


};
void Test::CSR_SET_function_test()
{
	bool all_function_fine = true;
	int rows = 3, cols = 3;
	double values1[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	int row_position1[4] = { 0, 3, 6,9 };
	int col_index1[9] = { 0,1,2,0,1,2,0,1,2 };
	auto* test_CSR1 = new CSRMatrix<double>(rows, cols, 9, values1, row_position1, col_index1);
	test_CSR1->setValue(0, 2, 100);

	double values2[8] = { 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	int row_position2[4] = { 0, 2, 5,8 };
	int col_index2[8] = { 1,2,0,1,2,0,1,2 };
	auto* test_CSR2 = new CSRMatrix<double>(rows, cols, 9, values2, row_position2, col_index2);

	test_CSR2->setValue(0, 0, 100);
	if (test_CSR1->getValue(0, 2) != 100 || test_CSR2->getValue(0, 0) != 100)
	{
		cout << "***   CSR SET FUNCTION TEST FAIL  ***" << endl;

	}

	else
	{

		cout << "***   CSR SET FUNCTION TEST SUCCESS    ***" << endl;
	}

	delete test_CSR1;
	delete test_CSR2;

};
void Test::CSR_GET_DIAG_test()
{
	int rows = 3, cols = 3;
	double values1[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	int row_position1[4] = { 0, 3, 6,9 };
	int col_index1[9] = { 0,1,2,0,1,2,0,1,2 };
	auto* test_CSR1 = new CSRMatrix<double>(rows, cols, 9, values1, row_position1, col_index1);

	auto* output1 = new CSRMatrix<double>(rows, cols, 3, true);

	test_CSR1->getMatDiag(*output1);

	if (output1->getValue(0, 0) != test_CSR1->getValue(0, 0) ||
		output1->getValue(1, 1) != test_CSR1->getValue(1, 1)
		|| output1->getValue(2, 2) != test_CSR1->getValue(2, 2))
	{
		cout << "***  CSR GET DIAG TEST FAIL! ***" << endl;

	}

	else
	{

		cout << "***  CSR GET DIAG TEST SUCCESS!  ***" << endl;
	}

	delete test_CSR1;
	delete output1;

}
void Test::CSR_GET_REM_test()
{
	int rows = 3, cols = 3;
	double values1[9] = { 43850, 18467, 6334, 26500, 46788, 15724, 11478, 29358, 44937 };
	int row_position1[4] = { 0, 3, 6,9 };
	int col_index1[9] = { 0,1,2,0,1,2,0,1,2 };
	auto* test_CSR1 = new CSRMatrix<double>(rows, cols, 9, values1, row_position1, col_index1);

	auto* output1 = new CSRMatrix<double>(rows, cols, 9, true);

	test_CSR1->getMatRemn(*output1);
	if (output1->getValue(0, 0) != 0
		|| output1->getValue(1, 1) != 0
		|| output1->getValue(2, 2) != 0)
	{
		cout << "***  CSR GET REM TEST FAIL! ***" << endl;

	}

	else
	{

		cout << "***  CSR GET REM TEST SUCCESS!   ***" << endl;
	}

	delete test_CSR1;
	delete output1;



};
void Test::CSR_MATMUL_test()
{
	int rows = 3, cols = 3;
	double values[9] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
	int row_position[4] = { 0, 3, 6,9 };
	int col_index[9] = { 0,1,2,0,1,2,0,1,2 };
	auto* test_CSR = new CSRMatrix<double>(rows, cols, 9, values, row_position, col_index);
	auto* output = new CSRMatrix<double>(rows, cols, 9, true);
	test_CSR->matMatMult(*test_CSR, *output);
	double result_value[9] = { 30,36,42,66,81,96,102,126,150 };
	bool all_same = true;
	for (int i = 0; i < test_CSR->nnzs; i++)
	{
		if (output->values[i] != result_value[i])
		{
			all_same = false;
		}

	}
	delete test_CSR;
	delete output;
	if (all_same)
	{
		cout << "***  CSR MATMUL TEST SUCCESS!    ***" << endl;

	}
	else
	{
		cout << "***  CSR MATMUL TEST FAIL!    ***" << endl;
	}

};

double * Test::Performance_Test(int size,int runtimes) 
{
	clock_t startTime, endTime;
	int minimum = 3;
	int maximum = 100;
	double* values = new double[size*size];
	/*size = 4;
	double values[16] = { 2,0,0,1,5,7,0,5,4,0,6,0,1,0,0,5 };*/
	double* b_vec = new double[size];
	double* x_vec = new double[size];
	double* timelist = new double[10];
	srand((unsigned int)time(NULL));
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double temp = (double)rand() / RAND_MAX;
			values[i * size + j] = minimum + temp * (maximum - minimum);
			if (values[i * size + i] <= (maximum - minimum) / 2) {
				values[i * size + i] += (maximum - minimum) / 1.8;
			}
		}
	}
	for (int i = 0; i < size; i++) {
		double temp = (double)rand() / RAND_MAX;
		b_vec[i] = minimum + temp * (maximum - minimum);
	}
	for (int i = 0; i < size; i++) 
	{
		x_vec[i] = 0;
	}
	auto* test_csr = new CSRMatrix<double>(size, size, values);
	auto* test_matrix = new Matrix<double>(size, size, values);
	int outer_counter = 0;
	int inner_counter = 0;
	auto* Jaco = new Solver<double>();
	auto* GS = new Solver<double>();
	auto* GE = new Solver<double>();
	auto* LU = new Solver<double>();
	auto* SOR = new Solver<double>();
	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes) 
	{
		
		Jaco->solverJacobi(200, 1e-10, 2, *test_matrix, b_vec, x_vec, inner_counter);
		outer_counter += 1;
	}
	endTime = clock();
	timelist[0] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	

	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes)
	{

		GS->solverGausSeid(200, 1e-10, 2, *test_matrix, b_vec, x_vec, inner_counter);
		outer_counter += 1;
	}
	endTime = clock();
	timelist[1] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	

	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes)
	{

		GE->solverGausElim(*test_matrix, b_vec, x_vec);
		outer_counter += 1;
	}
	endTime = clock();
	timelist[2] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	

	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes)
	{

		LU->solverLuDecom(*test_matrix, b_vec, x_vec);
		outer_counter += 1;
	}
	endTime = clock();
	timelist[3] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	
	
	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes)
	{

		SOR->solverSor(200, 1e-10, 2, *test_matrix, b_vec, x_vec, inner_counter);
		outer_counter += 1;
	}
	endTime = clock();
	
	timelist[4] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	
	
	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes)
	{

		Jaco->solverJacobi(200, 1e-10, 2, *test_csr, b_vec, x_vec, inner_counter);
		outer_counter += 1;
	}
	endTime = clock();
	timelist[5] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	

	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes)
	{

		GS->solverGausSeid(200, 1e-10, 2, *test_csr, b_vec, x_vec, inner_counter);
		outer_counter += 1;
	}
	endTime = clock();
	timelist[6] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	

	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes)
	{

		GE->solverGausElim(*test_csr, b_vec, x_vec);
		outer_counter += 1;
	}
	endTime = clock();
	timelist[7] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	

	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes)
	{

		LU->solverLuDecom(*test_csr, b_vec, x_vec);
		outer_counter += 1;
	}
	endTime = clock();
	timelist[8] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	

	outer_counter = 0;
	startTime = clock();
	while (outer_counter < runtimes)
	{

		SOR->solverSor(200, 1e-10, 2, *test_csr, b_vec, x_vec, inner_counter);
		outer_counter += 1;
	}
	endTime = clock();

	timelist[9] = (double)(endTime - startTime) / CLOCKS_PER_SEC;
	
	return timelist;
	delete[]values;
	delete[]b_vec;
	delete[]x_vec;
	delete test_csr;
	delete test_matrix;
	delete Jaco;
	delete GS;
	delete GE;
	delete LU;
	delete SOR;
};

void Test::Performance_Test_all() 
{
	cout << "NOTICE:\t This test may take few minutes to run because of large set of tests"
		<< endl << ">\t Testing order:"
		<< endl << ">\t Jacobi -> GS -> GE -> LU -> SOR (1 - 5) with DES matrix"
		<< endl << ">\t Jacobi to SOR with CSR matrix ги6 - 10)" 
		<< endl << endl << endl;

	double* testlist1;
	cout << "TEST for SIZE 10" << endl;
	cout << "Loading......." << endl;
	testlist1 = this->Performance_Test(10, 15);
	for (int i = 0; i < 10; i++)
	{
		cout << "RUN Time: ";
		cout << testlist1[i] << "  ";
		cout << "s" << endl;

	}
	
	double* testlist2;
	cout << "TEST for SIZE 50" << endl;
	cout << "Loading......." << endl;
	testlist2 = this->Performance_Test(50, 15);
	for (int i = 0; i < 10; i++)
	{
		cout << "RUN Time: ";
		cout << testlist2[i] << "  ";
		cout << "s" << endl;

	}
	double* testlist3;
	cout << "TEST for SIZE 80" << endl;
	cout << "Loading......." << endl;
	testlist3 = this->Performance_Test(80, 15);
	for (int i = 0; i < 10; i++)
	{
		cout << "RUN Time: ";
		cout << testlist3[i] << "  ";
		cout << "s" << endl;

	}
	cout << "TEST for SIZE 100" << endl;
	cout << "Loading......." << endl;
	double* testlist4;
	testlist4 = this->Performance_Test(100, 15);
	for (int i = 0; i < 10; i++)
	{
		cout << "RUN Time: ";
		cout << testlist4[i] << "  ";
		cout << "s" << endl;

	}
	delete []testlist1;
	delete[]testlist2;
	delete[]testlist3;
	delete[]testlist4;


};


template <typename T>
bool IsAlmostEqual(T* x, T* y, int length, int ulp)
{
	// the machine epsilon has to be scaled to the magnitude of the values used
	// and multiplied by the desired precision in ULPs (units in the last place)
	for (int i = 0; i < length; i++)
	{

		if (!(abs(x[i] - y[i]) < 1e-3))
		{
			return false;
		}
	}

	return true;
}

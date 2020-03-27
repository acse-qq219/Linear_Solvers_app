#pragma once

/**
	Implementation of user interface (class)

	This class defines the way to test the function of the project. Including solver and matrix function.

*/
class Test
{
public:
	
	/**Constructor of Test */
	Test();
	/**Destructor of Test*/
	~Test();
	
	/**
	The set of all tests
	*/
	void main_test();
	/**
	test Jacob with dense matrix 
	*/
	void Jacob_test();
	/**
	test GS with dense matrix
	*/
	void GausSeid_test();
	/**
	test GE with dense matrix
	*/
	void GausElim_test();
	/**
	test LU with dense matrix
	*/
	void LuDecom_test();
	/**
	test SOR with dense matrix
	*/
	void SOR_test();
	/**
	The set of all solver test
	*/
	void Solver_test();
	/**
	The set of all csr test
	*/
	void CSR_test();
	/**
	test Jacob with CSR matrix
	*/
	void CSR_jacob_test();
	/**
	test GS with CSR matrix
	*/
	void CSR_GausSeid_test();
	/**
	test GE with CSR matrix
	*/
	void CSR_GausElim_test();
	/**
	test LU with CSR matrix
	*/
	void CSR_LUDecom_test();
	/**
	test SOR with CSR matrix
	*/
	void CSR_SOR_test();
	/**
	test setValue function in CSR
	*/
	void CSR_SET_function_test();
	/**
	test the get diag function 
	*/
	void CSR_GET_DIAG_test();
	/**
	test the get remain function (without diag)
	*/
	void CSR_GET_REM_test();
	/**
	test MAT*MAT function on CSR
	*/
	void CSR_MATMUL_test();
	/**
	test Performance
	*/
	double* Performance_Test(int size,int runtimes);
	/**
	test all Performance
	*/
	void Performance_Test_all();
	
};

#pragma once
#include <string>
#include <iostream>
#include <ctime>

/**
	Implementation of user interface (class)

	This class defines the way to get inputs from users. The main process looks like:
	Introducton -> Select Input Data Source -> Select Solevr -> Select Output
	Only the `interfaceInvalid()` method takes `message` arguments.

	@tparam T The type of data stored in the table
*/
template <class T>
class Interface {

public:
	/** Constructor of Interface

		Create a Interface object.
	*/
	Interface();

	/** Deconstructor of Interface
	*/
	~Interface();

private:
	/** Print error message on the screen to users

		This method is often called when input is wrong
		or the status of function is wrong

		@throws IoError Thrown if input is wrong or status is wrong
	*/
	void interfaceInvalid(std::string message);

	/** Print introduction of the system

		This method will be called by Interface Constructor, it will
		print introduction information on the screen. Users can input
		'y' or 'n' to continue or exit. If 'y', it will call
		`interfaceSelectData()` method; if 'n', it will exit; Otherwise,
		it will call `interfaceInvalid()` with "Invalid Selection !" as 
		argument.
	*/
	void interfaceIntro();

	/** Select data source from users

		This method will be called by `interfaceIntro()`, it will
		print data selection options on the screen. Users can input
		'1', '2', 'b' or 'x' to generate a matrix with random number,
		read data from a txt file, back or exit. If '1', it will call
		`generateRandMat()` method; if '2', it will call 
		`readDataFromFile()`; if 'b', it will call `interfaceIntro()`.
	*/
	void interfaceSelectData();

	/** Generate a matrix with random values

		This method will be called by `interfaceSelectData()`, it will
		ask users to enter size of the matrix and min & max value as 
		random value thershold. those value will be stored to DESMatrix 
		attribute variables. The minimum allowed size is 2x2 and maximum 
		one is 1000 x 1000. It will determine if the matrix is a DESMatrix
		or CSRMatrix for later use. After generated successfully, it will 
		call `generateRandVec()` method.
	*/
	void generateRandMat();

	/** Generate a vector with random values

		This method will be called by `generateRandMat()`, it will
		ask users to enter the minimum and maximum value as random value 
		thershold. Those value will be stored to Vector attribute variables.
		After generated successfully, it will call `interfaceSelectSolver()`
		method.
	*/
	void generateRandVec();

	/** Read data from a file

		This method will be called by `interfaceSelectData()`, it will
		ask users to enter the path of a txt file. Those value will be stored 
		to Matrix and Vector attribute variables. Same as `generateRandMat()`,
		it will determine DESMatrix and CSRMatrix. After generated successfully, 
		it will call `interfaceSelectSolver()` method.
	*/
	void readDataFromFile();

	/** Select solvers from users

		This method will be called by `generateRandVec()` or `readDataFromFile()`, 
		it will ask users to enter '1', '2', '3', '4', '5', 'b' or 'x'. If '1', it will
		call `interfaceJacobi()`; if '2', it will call `interfaceGasSed()`; if '3', it 
		will call `interfaceGasElim()`; if '4', it will call `interfaceLUDecom()`; if '5',
		it will call `interfaceSOR()`. These methods correspond to different member in
		`Solver()` calss.
	*/
	void interfaceSelectSolver();

	/** Use Jacobi Method to solve the linear system
		
		This method will be called by `interfaceSelectSolver()`. It will ask users to 
		enter iteration times, initial gusess for solutions and allowed precision. These
		value will be passed to solver member jacobi method. The running time will be 
		recorded as well.
	*/
	void interfaceJacobi();

	/** Use Gauss-Sediel Method to solve the linear system

		This method will be called by `interfaceSelectSolver()`. It will ask users to
		enter iteration times, initial gusess for solutions and allowed precision. These
		value will be passed to solver member Gauss-Sediel method. The running time 
		will be recorded as well.
	*/
	void interfaceGasSed();

	/** Use SOR Method to solve the linear system

		This method will be called by `interfaceSelectSolver()`. It will print the matrix
		and vector first. Then it will ask users to enter iteration times, initial gusess
		for solutions and allowed precision. These value will be passed to solver member
		SOR method. The running time will be recorded as well.
	*/
	void interfaceSOR();

	/** Use Gaussian Elimination Method to solve the linear system

		This method will be called by `interfaceSelectSolver()`. It will print the matrix
		and vector first. Then it will call solver Gaussian Elimination member method.
		The running time will be recorded as well.
	*/
	void interfaceGasElim();

	/** Use LU Decomposition Method to solve the linear system

		This method will be called by `interfaceSelectSolver()`. It will print the matrix
		and vector first. Then it will call solver LU Decomposition member method. 
		The running time will be recorded as well.
	*/
	void interfaceLUDecom();

	/** Function to test all solvers

		This method will be called by `interfaceSelectSolver()`. It will print the size
		of test linear system and runnning time of all solvers.
	*/
	void interfaceComTest();

	/// Select output from user and output some stuff
	void interfaceSelectOutput();
	void printSolutionVec();
	void saveDataToFile();

	/// bool flag to determine if input matrix is a CSRMatrix
	bool is_csr_mat = false;		///< is CSRMatrix or not
	bool is_solver_right = true;	///< is the solver status right or not
	bool have_output_file = false;	///< have the system output a txt file 

	/// CSRMatrix attribute
	int nnzs_main = -1;						///< number of non-zeros in the matrix
	int* row_pos_csrmat_main = nullptr;     ///< pointer to row position array
	int* col_pos_csrmat_main = nullptr;		///< pointer to col position array
	T* value_csr_mat_main = nullptr;		///< pointer to csr matrix non-zero value array

	/// DESMatrix attribute
	int row_desmat_main = -1;				///< number of rows in the matrix
	int col_desmat_main = -1;				///< number of columns in the matrix
	T* value_mat_main = nullptr;			///< pointer to matrix A value array

	/// Vector attribute
	T* value_b_vec_main = nullptr;			///< pointer to vector b value array
	T* value_x_vec_main = nullptr;			///< pointer to vector x value array

	/// Time attribute
	clock_t start_time = 0;					///< number of the start time of solver
	clock_t end_time = 0;					///< number of the end time of solver

	/// Iteration attribute
	int iter_times = 1;						///< iteartion times

	/// Precision attribute
	int prec = 4;							///< precious to specific decimal places
};


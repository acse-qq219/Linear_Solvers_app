#pragma once
#include "Matrix.h"
#include "CSRMatrix.h"

/**
	Implementation of solvers (class)

	This class defines five solvers which calculate solutions array for linear systems
	entered from interface. The solver methods take `Matrix<T>` (or `CSRMatrix<T>`),
	`b_right_value` and `b_right_size` as most important arguments.

	@tparam T The type of data stored in the table
*/
template <class T>
class Solver {

public:
	/** Constructor of Solver

		Create a Solver object.
	*/
	Solver();

	/** Deconstructor of Solver
	*/
	~Solver();

	/** Use Iterative Jacobi Method to calculate solutions array with DESMatrix.

		@param iterations			  Number of Iterations for this algorithm
		@param allowed_convergence    Required precision (to decimal places)
		@param init_guess			  Initial guess for solutions array
		@param Matrix<T> a_left		  (dense or sparse) matrix with datatype T
		@param b_right_value		  Pointer to 1-D array b
		@param x_ans				  Pointer to 1-D array x
		@param counter				  The number of iterations that satisfy the condition

		@return status of the convergence, true if reached the required precision in
		specified iterations, otherwise false.

		@throws DiagElemZeroError Thrown if diagonal elements exist zero value.
	*/
	bool solverJacobi(int iterations, double allowed_convergence, T init_guess, Matrix<T>& a_left, T* b_right_value, T* x_ans, int& counter);
	
	/** Use Iterative Jacobi Method to calculate solutions array with CSRMatrix.

		@overload
	*/
	bool solverJacobi(int iterations, double allowed_convergence, T init_guess, CSRMatrix<T>& a_left, T* b_right_value, T* x_ans, int& counter);


	/** Use Iterative Gauss-Seidel Method to calculate solutions array with DESMatrix.

		@param iterations			  Number of Iterations for this algorithm
		@param allowed_convergence    Required precision (to decimal places)
		@param init_guess			  Initial guess for solutions array
		@param Matrix<T> a_left		  (dense or sparse) matrix with datatype T
		@param b_right_value		  Pointer to 1-D array b
		@param x_ans				  Pointer to 1-D array x
		@param counter				  The number of iterations that satisfy the condition

		@return status of the convergence, true if reached the required precision in
		specified iterations, otherwise false.

		@throws DiagElemZeroError Thrown if diagonal elements exist zero value.
	*/
	bool solverGausSeid(int iterations, double allowed_convergence, T init_guess, Matrix<T>& a_left, T* b_right_value, T* x_ans, int& counter);

	/** Use Iterative Gauss-Seidel Method to calculate solutions array with CSRMatrix.

		@overload
	*/
	bool solverGausSeid(int iterations, double allowed_convergence, T init_guess, CSRMatrix<T>& a_left, T* b_right_value, T* x_ans, int& counter);


	/** Use Iterative Successive Over-relaxation Method to calculate solutions array with DESMatrix.

		@param iterations			  Number of Iterations for this algorithm
		@param allowed_convergence    Required precision (to decimal places)
		@param Matrix<T> a_left		  (dense or sparse) matrix with datatype T
		@param b_right_value		  Pointer to 1-D array b
		@param x_ans				  Pointer to 1-D array x
		@param counter				  The number of iterations that satisfy the condition

		@return status of the convergence, true if reached the required precision in
		specified iterations, otherwise false.

		@throws DiagElemZeroError Thrown if diagonal elements exist zero value.
	*/
	bool solverSor(int iterations, double allowed_convergence, T init_guess, Matrix<T>& a_left, T* b_right_value, T* x_ans, int& counter);
	
	/** Use Iterative Successive Over-relaxation Method to calculate solutions array with CSRMatrix.

		@overload
	*/
	bool solverSor(int iterations, double allowed_convergence, T init_guess, CSRMatrix<T>& a_left, T* b_right_value, T* x_ans, int& counter);


	/** Use Direct Gaussian Elimination Method to calculate solutions array with DESMatrix.

		@param Matrix<T> a_left		  (dense or sparse) matrix with datatype T
		@param b_right_value		  Pointer to 1-D array b
		@param x_ans				  Pointer to 1-D array x

		@return status of the convergence, true if reached the required precision in
		specified iterations, otherwise false.

		@throws DiagElemZeroError Thrown if diagonal elements exist zero value.
	*/
	bool solverGausElim(Matrix<T>& a_left, T* b_right_value, T* x_ans);

	/** Use Direct Gaussian Elimination Method to calculate solutions array with CSRMatrix.

		@overload
	*/
	bool solverGausElim(CSRMatrix<T>& a_left, T* b_right_value, T* x_ans);


	/** Use Direct Gaussian Elimination Method to calculate solutions array with DESMatrix.

		@param Matrix<T> a_left		  (dense or sparse) matrix with datatype T
		@param b_right_value		  Pointer to 1-D array b
		@param x_ans				  Pointer to 1-D array x

		@return status of the convergence, true if reached the required precision in
		specified iterations, otherwise false.

		@throws DiagElemZeroError Thrown if diagonal elements exist zero value.
	*/
	bool solverLuDecom(Matrix<T>& a_left, T* b_right_value, T* x_ans);

	/** Use Direct Gaussian Elimination Method to calculate solutions array with CSRMatrix.

		@overload
	*/
	bool solverLuDecom(CSRMatrix<T>& a_left, T* b_right_value, T* x_ans);
};
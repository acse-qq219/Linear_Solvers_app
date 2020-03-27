#pragma once
#include "Matrix.h"

/**
	Implementation of CSRMatrix (class)

	This is the Matrix class for CSR Format stored matrix. This class
	includes basic operations like addition, multiplication, etc.
	and advanced operations like convert to Original Format.

	@tparam T The type of data stored in the table
*/
template <class T>
class CSRMatrix : public Matrix<T> {

public:
	/** Constructor of CSRMatrix

		Create a CSRMatrix object with empty values. If preallocate is
		true, allocate memory to this object.

		@param rows				number of rows of the matrix
		@param cols				number of columns of the matrix
		@param preallocate		determine if allocate memory to it
		@param nnzs				number of nonezero value
	*/
	CSRMatrix(int rows, int cols, int nnzs, bool preallocate);

	/** Constructor of CSRMatrix

		Create a CSRMatrix object with given values.

		@param rows				number of rows of the matrix
		@param cols				number of columns of the matrix
		@param values_ptr		1-D arry of values
		@param row_position		1-D arry of row_position indicator value (row off-set)
		@param col_index		1-D arry of col indicator 
	*/
	/** Constructor of CSRMatrix

		Create a CSRMatrix object with given values.Use the same format like Matrix

		@param rows				number of rows of the matrix
		@param cols				number of columns of the matrix
		@param values_ptr		1-D arry of values
		
	*/
	CSRMatrix(int rows, int cols, T* values_ptr);
	CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index);
	/** Destructor of CSRMatrix

		Clear temporary memory allocation.
	*/
	~CSRMatrix();

	/** 
	Print the whole matrix 
	*/
	virtual void printMatrix();

	/** Get the diagonal elements of the CSRmatrix

		Other elements are 0.

		@param output			the CSRMatrix to record diagonal elements

		@throws	DimensionNotMatchError Thrown if value array of `output` matrix
		is NULL and the size of original matrix is not same with that of `output`.
	*/
	void getMatDiag(CSRMatrix<T>& output);
	/** Get the remain elements of the CSRmatrix

		Diagonal elements are 0.

		@param output			the CSRMatrix to record remain elements

		@throws	DimensionNotMatchError Thrown if value array of `output` matrix
		is NULL and the size of original matrix is not same with that of `output`.
	*/
	void getMatRemn(CSRMatrix<T>& output);
	
	/** Get inverse of a CSRmatrix
	*/
	void calMatInver(CSRMatrix<T>& output);
	
	/** Get transpose of a matrix
	*/
	void calMatTrans(CSRMatrix<T>& output);

	/** CSRMatrix-vector dot multiplication

		@param vec_size			the number of vector size
		@param vec_right		the pointer to 1-D RHS vector value array
		@param output			the pointer to 1-D output vector value array
	*/
	void matVecMult(int vec_size, T* vec_right, T* output);
	
	/** CSRMatrix-CSRmatrix dot multiplication

		@param mat_right		the Matrix of RHS matrix
		@param output			the Matrix to record results
	*/
	void matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output);
	
	/** CSRMatrix-CSRmatrix addition
	*/
	void matMatAdd(CSRMatrix<T>& mat_right, CSRMatrix<T>& output);

	/** CSRMatrix-CSRmatrix subtraction
	*/
	void matMatSub(CSRMatrix<T>& mat_right, CSRMatrix<T>& output);

	/** Get the value of spcified position in a CSR matrix

		Find the index accoring to row and column index (i, j)

		@param row			the row index of the CSR matrix
		@param col			the column index of the CSR matrix
	*/
	T getValue(int row, int col);

	/** Set the value of spcified position in  a CSR matrix

		Find the index accoring to row and column index (i, j)

		@param row			the row index of the CSR matrix
		@param col			the column index of the CSR matrix
		@param value		the value number to be set
	*/
	void setValue(int row, int col,T value);
	
	/** Get the value of spcified position in a CSR matrix

		Find the index according to one argument
		(consider value stored as 1-D array)

		@overload
	*/
	T getValue(int position);
	
	/** Set the value of spcified position in a CSR matrix

		Find the index according to one argument
		(consider value stored as 1-D array)

		@overload
	*/
	void setValue(int position, T value);
	
	/** Print the CSR matrix as a normal one. For example:

		x x x
		x x x
		x x x
	*/
	void format_print();
	
	// Variables
	int* row_position = nullptr;
	int* col_index = nullptr;
	int nnzs = -1;
};
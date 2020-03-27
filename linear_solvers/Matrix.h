#pragma once

/**
	Implementation of Matrix (class)

	This is the Matrix class for densely stored matrix. This class
	includes basic operations like addition, multiplication, etc.
	and advanced operations like calculate inverse and transpose of
	a DESMatrix.

	@tparam T The type of data stored in the table
*/
template <class T>
class Matrix {

public:
	/** Constructor of Matrix

		Create a Matrix object with empty values. If preallocate is 
		true, allocate memory to this object.

		@param rows				number of rows of the matrix
		@param cols				number of columns of the matrix
		@param preallocate		determine if allocate memory to it
	*/
	Matrix(int rows, int cols, bool preallocate);

	/** Constructor of Matrix

		Create a Matrix object with given values. 

		@param rows				number of rows of the matrix
		@param cols				number of columns of the matrix
		@param values_ptr		1-D arry of values
	*/
	Matrix(int rows, int cols, T* values_ptr);

	/** Destructor of Matrix

		Clear temporary memory allocation.
	*/
	virtual ~Matrix();

	/** Print each value of the matrix within a line
	*/
	virtual void printValues();

	/** Print the matrix as a normal one. For example:

		x x x
		x x x
		x x x
	*/
	virtual void printMatrix();

	/** Get the size of the matrix

		@return size_of_values get the private variable size_of_values
	*/
	int getMatSize() const { return size_of_values; }

	/** Set the size of the matrix

		@param new_size			the new size we want to change
	*/
	void setMatSize(int new_size) { size_of_values = new_size; }

	/** Get the diagonal elements of the matrix

		Other elements are 0.

		@param output			the Matrix to record diagonal elements
		
		@throws	DimensionNotMatchError Thrown if value array of `output` matrix
		is NULL and the size of original matrix is not same with that of `output`.
	*/
	virtual void getMatDiag(Matrix<T>& output);

	/** Get the remain elements of the matrix

		Diagonal elements are 0.

		@param output			the Matrix to record remain elements

		@throws	DimensionNotMatchError Thrown if value array of `output` matrix
		is NULL and the size of original matrix is not same with that of `output`.
	*/
	virtual void getMatRemn(Matrix<T>& output);

	/** Get the lower triangular elements of the matrix

		Upper triangular elements are 0.

		@param output			the Matrix to record lower triangular elements

		@throws	DimensionNotMatchError Thrown if value array of `output` matrix
		is NULL and the size of original matrix is not same with that of `output`.
	*/
	virtual void getMatLower(Matrix<T>& output);

	/** Get the upper triangular elements of the matrix

		Lower triangular elements are 0.

		@param output			the Matrix to record upper triangular elements

		@throws	DimensionNotMatchError Thrown if value array of `output` matrix
		is NULL and the size of original matrix is not same with that of `output`.
	*/
	virtual void getMatUpper(Matrix<T>& output);
	
	/** Calculate the inverse of the matrix

		@param output			the Matrix to record the inverse of the matrix

		@throws	DimensionNotMatchError Thrown if value array of `output` matrix
		is NULL and the size of original matrix is not same with that of `output`.
		@throws NoSquareMatrixError Thrown if the original matrix is not a square
		matrix.
		@throws DiagElemZeroError Thrown if diagonal elements exist zero value.
	*/
	virtual void calMatInver(Matrix<T>& output);

	/** Calculate the transpose of the matrix

		@param output			the Matrix to record the transpose of the matrix

		@throws	DimensionNotMatchError Thrown if value array of `output` matrix
		is NULL and the size of original matrix is not same with that of `output`.
	*/
	virtual void calMatTrans(Matrix<T>& output);

	/** Matrix-vector dot multiplication

		@param vec_size			the number of vector size
		@param vec_right		the pointer to 1-D RHS vector value array 
		@param output			the pointer to 1-D output vector value array 
	*/
	virtual void matVecMult(int vec_size, T* vec_right, T* ouput);

	/** Matrix-matrix dot multiplication

		@param mat_right		the Matrix of RHS matrix
		@param output			the Matrix to record results
	*/
	virtual void matMatMult(Matrix<T>& mat_right, Matrix<T>& output);

	/** Matrix-matrix addition

		@param mat_right		the Matrix of RHS matrix
		@param output			the Matrix to record results
	*/
	virtual void matMatAdd(Matrix<T>& mat_right, Matrix<T>& output);

	/** Matrix-matrix substraction

		@param mat_right		the Matrix of RHS matrix
		@param output			the Matrix to record results
	*/
	virtual void matMatSub(Matrix<T>& mat_right, Matrix<T>& output);

	/** Get the value of spcified position in a matrix

		Find the index accoring to row and column index (i, j)

		@param row			the row index of the matrix
		@param col			the column index of the matrix
	*/
	virtual T getValue(int row, int col);

	/** Get the value of spcified position in a matrix

		Find the index according to one argument
		(consider value stored as 1-D array)

		@overload
	*/
	virtual T getValue(int position);

	/** Set the value of spcified position in a matrix 
	
		Find the index accoring to row and column index (i, j)

		@param row			the row index of the matrix
		@param col			the column index of the matrix
		@param value		the value number to be set
	*/
	virtual void setValue(int row, int col, T value);

	/** Set the value of spcified position in a matrix

		Find the index according to one argument
		(consider value stored as 1-D array)

		@overload
	*/
	virtual void setValue(int position, T value);

	T* values = nullptr;			///< pointer to Matrix value array
	int rows = -1;					///< number of rows in the matrix
	int cols = -1;					///< number of columns in the matrix

protected:
	bool preallocated = false;		///< boolean variable to determine if allocate memory

private:
	int size_of_values = -1;		///< number of size of the matrix
};
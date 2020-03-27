#include <iostream>
#include "CSRMatrix.h"

using namespace std;


template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate) : Matrix<T>(rows, cols, false), nnzs(nnzs) {
	// If don't pass false in the initialisation list base constructor, 
	// it would allocate values to size (rows * cols) in base matrix class
	// So then set it to the real value we had passed in
	this->preallocated = preallocate;

	// If preallocate memory to it
	if (this->preallocated) {
		this->values = new T[this->nnzs];
		this->row_position = new int[this->rows + 1];
		this->col_index = new int[this->nnzs];
	}
}

// Constructor - setting the value with (T) pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index) : Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index)
{}
// Constructor -use the same format like Matrix
template<class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, T* values_ptr) : Matrix<T>(rows, cols, false)
{
	this->row_position = new int[this->rows + 1];
	int* temp_col_index = new int[this->rows * this->cols];
	T* temp_values = new T[this->rows * this->cols];

	int nnzs = 0;
	for (int i = 0; i < this->rows; i++) {
		bool has_been_checked = false;
		bool header_checker = true;

		for (int j = 0; j < this->cols; j++) {
			if (values_ptr[i * this->cols + j] != 0)
			{
				temp_values[nnzs] = values_ptr[i * this->cols + j];
				temp_col_index[nnzs] = j;

				if (header_checker)
				{
					this->row_position[i] = nnzs;
					header_checker = false;
					has_been_checked = true;
				}
				nnzs += 1;

			}
		}
		if (!has_been_checked)
		{
			this->row_position[i] = nnzs;

		}
	}
	this->row_position[this->rows] = nnzs;
	this->nnzs = nnzs;
	this->values = new T[nnzs];
	this->col_index = new int[nnzs];
	for (int i = 0; i < nnzs; i++) {
		this->values[i] = temp_values[i];
		this->col_index[i] = temp_col_index[i];
	}
	delete[]temp_col_index;
	delete[]temp_values;


};
template <class T>
CSRMatrix<T>::~CSRMatrix() {
	// Delete values arrays if preallocated memory
	if (this->preallocated) {
		delete[] this->row_position;
		delete[] this->col_index;
	}
	// Destructor is called after finish
}

template <class T>
void CSRMatrix<T>::printMatrix() {
	cout << "Printing matrix" << endl << endl;
	cout << "values: ";
	for (int i = 0; i < this->nnzs; i++) {
		cout << this->values[i] << " ";
	}
	cout << endl;
	cout << "row_position: ";
	for (int i = 0; i < this->rows + 1; i++) {
		cout << this->row_position[i] << " ";
	}
	cout << endl;
	cout << "col_index: ";
	for (int i = 0; i < this->nnzs; i++)
	{
		cout << this->col_index[i] << " ";
	}
	cout << endl;
}


// Get diagonal elements of a matrix
template <class T>
void CSRMatrix<T>::getMatDiag(CSRMatrix<T>& output) {
	if (output.values != nullptr) {
		// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	for (int i = 0; i < this->rows; i++)
	{
		output.values[i] = 0;
		output.col_index[i] = i;
		output.row_position[i] = i;

	}
	output.row_position[this->rows] = this->rows;
	for (int i = 0; i < this->rows; i++)
	{
		output.setValue(i, i, this->getValue(i, i));

	}

}

// Get remain elements of a matrix
template <class T>
void CSRMatrix<T>::getMatRemn(CSRMatrix<T>& output) {
	if (output.values != nullptr) {
		// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	for (int i = 0; i < this->nnzs; i++)
	{
		output.values[i] = this->values[i];
		output.col_index[i] = this->col_index[i];
	}
	for (int i = 0; i < this->rows + 1; i++)
	{
		output.row_position[i] = this->row_position[i];
	}
	for (int i = 0; i < this->rows; i++)
	{
		output.setValue(i, i, 0);

	}

}

// Get inverse of a matrix
template <class T>
void CSRMatrix<T>::calMatInver(CSRMatrix<T>& output) {
	//// Not implement yet ////
}

// Get transpose of a matrix
template <class T>
void CSRMatrix<T>::calMatTrans(CSRMatrix<T>& output) {
	//// Not implement yet ////
}

// Do matrix-vector dot multiplication (sparsely stored matrix)
// output = this .* vec_right
template <class T>
void CSRMatrix<T>::matVecMult(int vec_size, T* vec_right, T* output) {
	if (vec_right == nullptr || output == nullptr)
	{
		cerr << "RHS or output vector haven't been created" << endl;
		return;
	}

	// Check our dimensions match
	if (this->rows != vec_size) {
		cerr << "Input LHS matrix and RHS vector dimensions do not match" << endl;
		return;
	}
	for (int i = 0; i < this->rows; i++)
	{
		output[i] = 0;
	}
	for (int i = 0; i < this->rows; i++)
	{
		// Loop over all the entries in this col
		for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
		{
			// This is an example of indirect addressing
			// Can make it harder for the compiler to vectorise!
			output[i] += this->values[val_index] * vec_right[this->col_index[val_index]];
		}
	}
}

// Do matrix-matrix dot multiplication (sparsely stored matrix)
// output = this .* mat_right
template <class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output) {

	// Check our dimensions match
	if (this->cols != mat_right.rows) {
		cerr << "Input dimensions for RHS matrix don't match" << endl;
		return;
	}

	// Check if our output matrix has had space allocated to it
	if (output.values != nullptr) {
		// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			cerr << "Input dimensions for output matrix don't match" << endl;
			return;
		}
	}
	// The output hasn't been preallocated, so we are going to do that
	else {
		cerr << "Output hasn't been allocated" << std::endl;
	}

	int* row_position = new int[this->rows + 1];
	int* col_index = new int[this->rows * output.cols];
	T* value = new T[this->rows * output.cols];
	int counter = 0;
	// HOW DO WE SET THE SPARSITY OF OUR OUTPUT MATRIX HERE??
	

	for (int i = 0; i < this->rows; i++)
	{
		bool headchecker = true;
		bool have_not_been_checked = true;
		for (int j = 0; j < mat_right.cols; j++)
		{
			T temp = 0;
			for (int k = 0; k < mat_right.rows; k++)
			{
				temp += this->getValue(i, k) * mat_right.getValue(k, j);

			}
			if (temp != 0)
			{

				value[counter] = temp;
				col_index[counter] = j;
				if (headchecker)
				{
					row_position[i] = counter;
					headchecker = false;
					have_not_been_checked = false;
				}
				counter = counter + 1;
			}

		}
		if(have_not_been_checked)
		{
			row_position[i] = counter;
		}

	}
	row_position[this->rows] = counter;
	output.nnzs = counter;
	for (int i = 0; i < counter; i++)
	{

		output.col_index[i] = col_index[i];
		output.values[i] = value[i];
	}
	output.row_position = row_position;
}

// Do matrix-matrix addition 
template <class T>
void CSRMatrix<T>::matMatAdd(CSRMatrix<T>& mat_right, CSRMatrix<T>& output) {

	if (this->rows != mat_right.rows || this->cols != mat_right.cols) {
		cerr << "Input dimensions for RHS matrix don't match" << endl;
		return;
	}

}

// Do matrix-matrix subtraction
template <class T>
void CSRMatrix<T>::matMatSub(CSRMatrix<T>& mat_right, CSRMatrix<T>& output) {

	if (this->rows != mat_right.rows || this->cols != mat_right.cols) {
		cerr << "Input dimensions for RHS matrix don't match" << endl;
		return;
	}
	//// not implement yet ////
}

template <class T>
T CSRMatrix<T>::getValue(int row, int col)
{
	bool notin = true;
	for (int k = this->row_position[row]; k < this->row_position[row + 1]; k++) {
		if (this->col_index[k] == col) {
			return this->values[k];
			notin = false;
		}
	}
	if (notin) {
		return 0;
	}


}

template <class T>
void CSRMatrix<T>::setValue(int row, int col, T value)
{
	bool inside = false;
	int position = this->row_position[row];
	for (int k = this->row_position[row]; k < this->row_position[row + 1]; k++) {
		if (col > this->col_index[k])
		{
			position = k + 1;
		}
		else if (col == this->col_index[k]) {
			inside = true;
			position = k;
		}
	}
	if (inside) {

		this->values[position] = value;

	}
	else
	{

		for (int i = row + 1; i < this->rows + 1; i++)
		{
			this->row_position[i] += 1;
		}
		T* temp_value = new T[this->nnzs + 1];
		int* temp_col = new int[this->nnzs + 1];
		for (int i = 0; i < position; i++) {
			temp_value[i] = this->values[i];
			temp_col[i] = this->col_index[i];
		}
		temp_value[position] = value;
		temp_col[position] = col;
		for (int i = position + 1; i < this->nnzs + 1; i++) {
			temp_value[i] = this->values[i - 1];
			temp_col[i] = this->col_index[i - 1];
		}
		/*delete[] this->values;
		delete[] this->col_index;*/
		this->col_index = new int[this->nnzs + 1];
		this->values = new T[this->nnzs + 1];
		for(int i = 0;i<this->nnzs+1;i++)
		{
			this->values[i]= temp_value[i];
			this->col_index[i]= temp_col[i];
		
		}
		this->nnzs = this->nnzs + 1;
		delete[] temp_value;
		delete[] temp_col;
	}

}
template <class T>
T CSRMatrix<T>::getValue(int position)
{
	int row = position / this->cols;
	int col = position % this->cols;
	return this->getValue(row, col);

};
template <class T>
void CSRMatrix<T>::setValue(int position, T value)
{
	int row = position / this->cols;
	int col = position % this->cols;
	this->setValue(row, col, value);

};

template <class T>
void CSRMatrix<T>::format_print()
{
	cout << "Print Matrix" << endl;
	for (int i = 0; i < this->rows; i++) {
		cout << endl;
		for (int j = 0; j < this->cols; j++) {
			bool notin = true;
			for (int k = this->row_position[i]; k < this->row_position[i + 1]; k++) {
				if (this->col_index[k] == j) {
					cout << this->values[k] << " ";
					notin = false;
				}
			}
			if (notin) {
				cout << "0" << " ";
			}
		}

	}
	cout << endl;

}



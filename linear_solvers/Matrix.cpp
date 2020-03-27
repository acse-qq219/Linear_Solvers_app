#include <iostream>
#include <vector>
#include <iomanip>
#include "Matrix.h"

/// Constructor - using an initialisation list here
template <class T>
Matrix<T>::Matrix(int rows, int cols, bool preallocate) : rows(rows), cols(cols), size_of_values(rows* cols), preallocated(preallocate) {
	/// If we want to handle memory ourselves
	if (this->preallocated) {
		this->values = new T[size_of_values];
	}
}

/// Constructor - setting the value of the T pointer
template <class T>
Matrix<T>::Matrix(int rows, int cols, T* values_ptr) : rows(rows), cols(cols), size_of_values(rows* cols), values(values_ptr)
{}

/// Destructor
template <class T>
Matrix<T>::~Matrix() {
	/// Delete the values array
	if (this->preallocated) {
		delete[] this->values;
	}
}

/// Print all vaules of the matrix in one line
template <class T>
void Matrix<T>::printValues() {
	std::cout << "Printing values" << std::endl;
	for (int i = 0; i < this->size_of_values; i++) {
		std::cout << this->values[i] << " ";
	}
	std::cout << std::endl;
}

/// Print what the matrix looks like
template <class T>
void Matrix<T>::printMatrix() {
	std::cout << std::endl << ">\t Print matrix:" << std::endl;
	for (int j = 0; j < this->rows; j++) {
		std::cout << std::endl << "row "<< j + 1 << "\t";
		for (int i = 0; i < this->cols; i++) {
			/// We have explicitly assumed row-major ordering here
			std::cout << this->values[i + j * this->cols] << " ";
		}
	}

}

template <class T>
void Matrix<T>::getMatDiag(Matrix<T>& output) {
	/// Check if the output matrix has had space allocated to it
	if (output.values != nullptr) {
		/// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	/// The output hasn't been preallocated
	else {
		output.values = new T[this->size_of_values];
		output.preallocated = true;
	}

	/// Set values to zero beforehand
	for (int i = 0; i < this->size_of_values; i++) {
		output.values[i] = 0;
	}

	for (int i = 0; i < this->rows; i++) {
		output.values[i * output.cols + i] = this->values[i * this->cols + i];
	}
}

template <class T>
void Matrix<T>::getMatRemn(Matrix<T>& output) {
	/// Check if the output matrix has had space allocated to it
	if (output.values != nullptr) {
		/// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	/// The output hasn't been preallocated
	else {
		output.values = new T[this->size_of_values];
		output.preallocated = true;
	}

	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			if (i != j) {
				output.values[i * output.cols + j] = this->values[i * this->cols + j];
			}
			else {
				output.values[i * output.cols + j] = 0;
			}
		}
	}
}

/// Get lower triangular component (NOT lower triangular matrix)
template <class T>
void Matrix<T>::getMatLower(Matrix<T>& output) {
	/// Check if the output matrix has had space allocated to it
	if (output.values != nullptr) {
		/// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	/// The output hasn't been preallocated
	else {
		output.values = new T[this->size_of_values];
		output.preallocated = true;
	}


	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			if (j <= i) {
				output.values[i * output.cols + j] = this->values[i * this->cols + j];
			}
			else {
				output.values[i * output.cols + j] = 0;
			}
		}
	}
}

template <class T>
void Matrix<T>::getMatUpper(Matrix<T>& output) {
	/// Check if the output matrix has had space allocated to it
	if (output.values != nullptr) {
		/// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	/// The output hasn't been preallocated
	else {
		output.values = new T[this->size_of_values];
		output.preallocated = true;
	}

	for (int i = 0; i < this->rows; i++) {
		for (int j = 0; j < this->cols; j++) {
			if (j > i) {
				output.values[i * output.cols + j] = this->values[i * this->cols + j];
			}
			else {
				output.values[i * output.cols + j] = 0;
			}
		}
	}
}

/// Calculate inverse of a matrix
template <class T>
void Matrix<T>::calMatInver(Matrix<T>& output) {
	if (this->rows != this->cols) {
		std::cerr << "The matrix A shoule be a square matrix" << std::endl;
		return;
	}

	for (int i = 0; i < this->rows; i++) {
		if (this->values[i * this->rows + i] == 0) {
			std::cerr << "The diagonal of input matrix cannot have zero" << std::endl;
			return;
		}
	}

	/// Check if the output matrix has had space allocated to it
	if (output.values != nullptr) {
		/// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	/// The output hasn't been preallocated
	else {
		output.values = new T[this->size_of_values];
		output.preallocated = true;
	}

	int n = this->cols;

	auto* mat_L = new Matrix<T>(n, n, true);
	auto* mat_U = new Matrix<T>(n, n, true);
	auto* mat_L_INV = new Matrix<T>(n, n, true);
	auto* mat_U_INV = new Matrix<T>(n, n, true);

	/// Set values to zero beforehand
	for (int i = 0; i < this->size_of_values; i++) {
		output.values[i] = 0;
		mat_L->values[i] = 0;
		mat_U->values[i] = 0;
		mat_L_INV->values[i] = 0;
		mat_U_INV->values[i] = 0;
	}

	for (int i = 0; i < n; i++) {
		mat_U->values[i] = this->values[i];
	}
	for (int i = 1; i <= n; i++) {
		mat_L->values[i * n - n] = this->values[i * n - n] / mat_U->values[0];
	}
	for (int r = 1; r < n; r++) {
		for (int i = r; i < n; i++) {
			double sum1 = 0;
			for (int k = 0; k < r; k++) {
				sum1 += mat_L->values[r * n + k] * mat_U->values[k * n + i];
			}
			mat_U->values[r * n + i] = this->values[r * n + i] - sum1;
		}
		if (r != n) {
			for (int i = r; i < n; i++) {
				double sum2 = 0;
				for (int k = 0; k < r; k++) {
					sum2 += mat_L->values[i * n + k] * mat_U->values[k * n + r];
				}
				mat_L->values[i * n + r] = (this->values[i * n + r] - sum2) / mat_U->values[r * n + r];
			}
		}
	}

	for (int i = 0; i < mat_U->rows; i++) {
		if (mat_U->values[i * mat_U->rows + i] == 0) {
			std::cerr << "Inverse matrix doesn't exist" << std::endl;
			return;
		}
	}

	for (int i = 0; i < n; i++) {
		mat_U_INV->values[i * n + i] = 1 / mat_U->values[i * n + i];
		for (int k = i - 1; k >= 0; k--) {
			double s = 0.00;
			for (int j = k + 1; j <= i; j++) {
				s += mat_U->values[k * n + j] * mat_U_INV->values[j * n + i];
			}
			mat_U_INV->values[k * n + i] = -s / mat_U->values[k * n + k];
		}
	}

	for (int i = 0; i < n; i++) {
		mat_L_INV->values[i * n + i] = 1 / mat_L->values[i * n + i];
		for (int k = i + 1; k < n; k++) {
			for (int j = i; j <= k - 1; j++) {
				mat_L_INV->values[k * n + i] -= mat_L->values[k * n + j] * mat_L_INV->values[j * n + i];
			}
		}
	}
	for (int i = 0; i < mat_U_INV->rows; i++) {
		for (int k = 0; k < mat_U_INV->cols; k++) {
			for (int j = 0; j < mat_L_INV->cols; j++) {
				output.values[i * output.cols + j] += mat_U_INV->values[i * mat_U_INV->cols + k] * mat_L_INV->values[k * mat_L_INV->cols + j];
			}
		}
	}

	delete mat_L;
	delete mat_U;
	delete mat_L_INV;
	delete mat_U_INV;
}

/// Calculate transpose of a matrix
template <class T>
void Matrix<T>::calMatTrans(Matrix<T>& output) {
	/// Check if the output matrix has had space allocated to it
	if (output.values != nullptr) {
		/// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	/// The output hasn't been preallocated
	else {
		output.values = new T[this->size_of_values];
		output.preallocated = true;
	}

	/// Set values to zero beforehand
	for (int i = 0; i < this->size_of_values; i++) {
		output.values[i] = 0;
	}

	int n = this->rows;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			output.values[j * n + i] = this->values[i * n + j];
		}
	}
}

/// Do matrix-vector dot multiplication (Densely stored matrix)
/// output = this .* vec_right
template <class T>
void Matrix<T>::matVecMult(int vec_size, T* vec_right, T* output) {
	/// Check if RHS and output vector have been assigned value to it
	if (vec_right == nullptr || output == nullptr) {
		std::cerr << "Input or output vector haven't been created" << std::endl;
		return;
	}
	/// Check our dimensions match
	if (this->rows != vec_size) {
		std::cerr << "Input LHS matrix and RHS vector dimensions do not match" << std::endl;
		return;
	}

	for (int i = 0; i < vec_size; i++) {
		output[i] = 0;
	}

	/// Do matrix-vector multiplication
	for (int i = 0; i < this->rows; i++) {
		for (int k = 0; k < this->cols; k++) {
			output[i] += this->values[i * this->cols + k] * vec_right[k];
		}
	}
}

/// Do matrix-matrix dot multiplication (Densely stored matrix)
/// output = this .* mat_right
template <class T>
void Matrix<T>::matMatMult(Matrix<T>& mat_right, Matrix<T>& output) {
	/// Check our dimensions match
	if (this->rows != mat_right.rows) {
		std::cerr << "Input dimensions do not match" << std::endl;
		return;
	}

	/// Check if the output matrix has had space allocated to it
	if (output.values != nullptr) {
		/// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	/// The output hasn't been preallocated
	else {
		output.values = new T[this->size_of_values];
		output.preallocated = true;
	}

	/// Set values to zero beforehand
	for (int i = 0; i < this->size_of_values; i++) {
		output.values[i] = 0;
	}

	for (int i = 0; i < this->rows; i++) {
		for (int k = 0; k < this->cols; k++) {
			for (int j = 0; j < mat_right.cols; j++) {
				output.values[i * output.cols + j] += this->values[i * this->cols + k] * mat_right.values[k * mat_right.cols + j];
			}
		}
	}
}

/// Do matrix-matrix addition (Densely stored matrix)
/// output = this + mat_right
template <class T>
void Matrix<T>::matMatAdd(Matrix<T>& mat_right, Matrix<T>& output) {
	/// Check our dimensions match
	if (this->rows != mat_right.rows || this->cols != mat_right.cols) {
		std::cerr << "Input dimensions do not match" << std::endl;
		return;
	}

	/// Check if the output matrix has had space allocated to it
	if (output.values != nullptr) {
		/// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	/// The output hasn't been preallocated
	else {
		output.values = new T[this->size_of_values];
		output.preallocated = true;
	}

	/// Set values to zero beforehand
	for (int i = 0; i < this->size_of_values; i++) {
		output.values[i] = 0;
	}

	for (int i = 0; i < this->size_of_values; i++) {
		output.values[i] += this->values[i] + mat_right.values[i];
	}
}

/// Do matrix-matrix subtraction (Densely stored matrix)
/// output = this - mat_right
template <class T>
void Matrix<T>::matMatSub(Matrix<T>& mat_right, Matrix<T>& output) {
	/// Check our dimensions match
	if (this->rows != mat_right.rows || this->cols != mat_right.cols) {
		std::cerr << "Input dimensions do not match" << std::endl;
		return;
	}

	/// Check if the output matrix has had space allocated to it
	if (output.values != nullptr) {
		/// Check our dimensions match
		if (this->rows != output.rows || this->cols != output.cols) {
			std::cerr << "Input dimensions for output matrix do not match" << std::endl;
			return;
		}
	}
	/// The output hasn't been preallocated
	else {
		output.values = new T[this->size_of_values];
		output.preallocated = true;
	}

	/// Set values to zero beforehand
	for (int i = 0; i < this->size_of_values; i++) {
		output.values[i] = 0;
	}

	for (int i = 0; i < this->size_of_values; i++) {
		output.values[i] += this->values[i] - mat_right.values[i];
	}
}

template<class T>
T Matrix<T>::getValue(int row, int col) 
{
	return this->values[row * this->cols + col];

};
template<class T>

void Matrix<T>::setValue(int row, int col, T value) 
{

	this->values[row * this->cols + col] = value;

};
template<class T>
T Matrix<T>::getValue(int position) 
{

	return this->values[position];
};
template<class T>
void Matrix<T>::setValue(int position, T value) 
{
	this->values[position] = value;

};
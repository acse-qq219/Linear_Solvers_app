#include <iostream>
#include <string>
#include <vector>
#include "Solver.h"

#define COUNT_TIME_MAX 10001

using namespace std;

template <class T>
Solver<T>::Solver() {}

template <class T>
Solver<T>::~Solver() {}

template <class T>
bool Solver<T>::solverJacobi(int iterations, double allowed_convergence, T init_guess, Matrix<T>& a_left, T* b_right_value, T* x_ans, int& counter) {
	// Assign all solutions for the initial guess
	for (int i = 0; i < a_left.rows; i++) {
		x_ans[i] = init_guess;
	}

	auto* diag_mat = new Matrix<T>(a_left.rows, a_left.cols, true);
	auto* remain_mat = new Matrix<T>(a_left.rows, a_left.cols, true);

	a_left.getMatDiag(*diag_mat);
	a_left.getMatRemn(*remain_mat);

	// Get the inverse of diagonal matrix
	for (int i = 0; i < a_left.rows; i++) {
		if (a_left.getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value"
				<< "\n No solution";
			return false;
		}
		else {
			diag_mat->setValue(i, i, 1. / diag_mat->getValue(i, i));
		}
	}

	// Result of inverse of diagonal matrix dot multiply (vector) b
	// output_db = inv(D) .* b
	T* output_db = new T[diag_mat->rows];
	diag_mat->matVecMult(diag_mat->rows, b_right_value, output_db);
	// Get negative of inverse of diagonal matrix for later use
	for (int i = 0; i < diag_mat->rows; i++) {
		if (diag_mat->getValue(i, i) != 0) {
			diag_mat->setValue(i, i, (-diag_mat->getValue(i, i)));
		}
	}

	counter = 0;
	do {

		// Result of negative inverse of diagonal matrix dot multiple remainder matrix
		// output_dr = -inv(D) .* R
		auto* output_dr = new Matrix<T>(diag_mat->rows, remain_mat->cols, true);
		diag_mat->matMatMult(*remain_mat, *output_dr);

		// Result of output_dr matrix dot multiply (vector) x
		// output_drx = output_dr .* x
		T* output_drx = new T[output_dr->rows];
		output_dr->matVecMult(a_left.rows, x_ans, output_drx);

		// Result of output_drx add output_db
		// output_add = output_drx + output_db
		T* output_add = new T[output_dr->rows];
		for (int i = 0; i < output_dr->rows; i++) {
			output_add[i] = output_drx[i] + output_db[i];
		}

		T err = 0;
		int err_counter = 0;
		for (int i = 0; i < a_left.rows; i++) {
			err = abs(output_add[i] - x_ans[i]);
			x_ans[i] = output_add[i];
			if (err <= allowed_convergence) {
				err_counter++;
			}
		}

		counter++;
		if (err_counter == a_left.rows) {
			delete diag_mat;
			delete remain_mat;
			delete output_dr;
			delete[] output_drx;
			delete[] output_db;
			delete[] output_add;

			return true;
		}

		delete output_dr;
		delete[] output_drx;
		delete[] output_add;

	} while (counter < iterations);

	delete[] output_db;
	delete diag_mat;
	delete remain_mat;

	return false;
}

template <class T>
bool Solver<T>::solverGausSeid(int iterations, double allowed_convergence, T init_guess, Matrix<T>& a_left, T* b_right_value, T* x_ans, int& counter) {
	// Assign all solutions for the initial guess
	T* temp_x = new T[a_left.rows];
	for (int i = 0; i < a_left.rows; i++) {
		x_ans[i] = init_guess;
		temp_x[i] = 0;
	}

	auto* diag_mat = new Matrix<T>(a_left.rows, a_left.cols, true);
	a_left.getMatDiag(*diag_mat);

	// Get the negative inverse of diagonal matrix
	for (int i = 0; i < diag_mat->rows; i++) {
		if (diag_mat->getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value"
				<< "\n No solution";
			return false;
		}
		else {
			diag_mat->setValue(i, i, 1. / diag_mat->getValue(i, i));
		}
	}

	T temp_L, temp_U, temp_sum, err;
	counter = 0;
	int err_counter = 0; // if err is less than allowed convergence, err_counter + 1
	do {
		for (int i = 0; i < a_left.rows; i++) {

			temp_L = temp_U = temp_sum = err = 0;

			for (int j = 0; j < i; j++) {
				temp_L += a_left.getValue(i, j) * x_ans[j];

			}
			for (int k = a_left.cols - 1; k >= i + 1; k--) {
				temp_U += a_left.getValue(i, k) * x_ans[k];
			}
			temp_sum = b_right_value[i] - temp_L - temp_U;
			temp_x[i] = diag_mat->getValue(i, i) * temp_sum;

			err = abs(temp_x[i] - x_ans[i]);
			x_ans[i] = temp_x[i];

			if (err <= allowed_convergence) {
				err_counter++;
			}
		}

		counter++;
		if (err_counter == a_left.rows) {
			delete diag_mat;
			return true;
		}

	} while (counter < iterations);

	delete diag_mat;
	return false;
}

template <class T>
bool Solver<T>::solverGausElim(Matrix<T>& a_left, T* b_right_value, T* x_ans) {

	for (int i = 0; i < a_left.rows; i++) {
		if (a_left.getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value" << endl
				<< "No solution" << endl << endl;
			system("pause");
			return false;
		}
	}

	// Elimination
	T temp = 0.0;
	for (int k = 0; k < a_left.rows - 1; k++) {
		for (int i = k + 1; i < a_left.rows; i++) {
			temp = a_left.getValue(i, k) / a_left.getValue(k, k);

			for (int j = k; j < a_left.rows; j++) {
				a_left.setValue(i, j, a_left.getValue(i, j) - temp * a_left.getValue(k, j));
			}
			b_right_value[i] -= temp * b_right_value[k];
		}
	}

	//Back substitution
	for (int k = a_left.rows - 1; k >= 0; k--) {
		temp = 0.0;
		for (int j = k + 1; j < a_left.rows; j++) {
			temp += a_left.getValue(k, j) * x_ans[j];

		}

		x_ans[k] = (b_right_value[k] - temp) / a_left.getValue(k, k);
	}

	return true;
}

template <class T>
bool Solver<T>::solverLuDecom(Matrix<T>& a_left, T* b_right_value, T* x_ans) {

	for (int i = 0; i < a_left.rows; i++) {
		if (a_left.getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value" << endl
				<< "No solution" << endl << endl;
			system("pause");
			return false;
		}
	}

	T* mat_U = new T[a_left.rows * a_left.cols];
	T* mat_L = new T[a_left.rows * a_left.cols];
	T* mat_y = new T[a_left.rows];

	for (int i = 0; i < a_left.rows * a_left.cols; i++) {
		mat_L[i] = 0;
		mat_U[i] = 0;
	}

	for (int i = 0; i < a_left.rows; i++) {
		mat_y[i] = 0;
	}

	for (int i = 0; i < a_left.rows; i++) {
		mat_U[i] = a_left.getValue(0, i);
	}

	for (int i = 1; i < a_left.rows; i++) {
		mat_L[i * a_left.rows] = a_left.getValue(i * a_left.rows) / mat_U[0];
	}

	for (int i = 1; i < a_left.rows; i++) {
		for (int k = i; k < a_left.rows; k++) {
			T sum1 = 0;
			for (int j = 0; j < i; j++) {
				sum1 += mat_L[i * a_left.rows + j] * mat_U[j * a_left.rows + k];
			}
			mat_U[i * a_left.rows + k] = a_left.getValue(i, k) - sum1;
		}
		if (i != a_left.rows - 1) {
			for (int k = i; k < a_left.rows; k++) {
				T sum2 = 0;
				for (int j = 0; j < i; j++) {
					sum2 += mat_L[k * a_left.rows + j] * mat_U[j * a_left.rows + i];
				}
				mat_L[k * a_left.rows + i] = (a_left.getValue(k, i) - sum2) / mat_U[i * a_left.rows + i];
			}
		}
	}
	for (int i = 0; i < a_left.rows; i++) {
		mat_L[i * a_left.rows + i] = 1;
	}

	mat_y[0] = b_right_value[0];
	for (int i = 1; i < a_left.rows; i++) {
		T sum3 = 0;
		for (int k = 0; k < i; k++) {
			sum3 += mat_L[i * a_left.rows + k] * mat_y[k];
		}
		mat_y[i] = b_right_value[i] - sum3;
	}

	x_ans[a_left.rows - 1] = mat_y[a_left.rows - 1] / mat_U[a_left.rows * a_left.rows - 1];
	for (int i = a_left.rows - 2; i >= 0; i--) {
		T sum4 = 0;
		for (int k = i + 1; k < a_left.rows; k++) {
			sum4 += mat_U[i * a_left.rows + k] * x_ans[k];
		}
		x_ans[i] = (mat_y[i] - sum4) / mat_U[i * a_left.rows + i];
	}

	delete[] mat_U;
	delete[] mat_L;
	delete[] mat_y;

	return true;
}


template <class T>
bool Solver<T>::solverSor(int iterations, double allowed_convergence, T init_guess, Matrix<T>& a_left, T* b_right_value, T* x_ans, int& counter) {

	for (int i = 0; i < a_left.rows; i++) {
		if (a_left.getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value" << endl
				<< "No solution" << endl << endl;
			system("pause");
			return false;
		}
	}

	T w = 1.46;
	T* vec_y = new T[a_left.rows];
	for (int i = 0; i < a_left.rows; i++) {
		vec_y[i] = 0;
	}

	for (int i = 0; i < a_left.rows; i++) {
		x_ans[i] = init_guess;
	}

	counter = 0;
	do {
		T err = 0;
		for (int i = 0; i < a_left.rows; i++) {
			T s = 0;
			for (int j = 0; j < a_left.rows; j++)
				if (j != i)
					s += a_left.getValue(i, j) * x_ans[j];
			vec_y[i] = (b_right_value[i] - s) / a_left.getValue(i, i);
			vec_y[i] = (1 - w) * x_ans[i] + w * vec_y[i];
			err = abs(x_ans[i] - vec_y[i]);
			x_ans[i] = vec_y[i];
			if (err <= allowed_convergence) break;
		}
		counter++;
	} while (counter < iterations);

	if (counter == iterations) {
		delete[] vec_y;
		return false;
	}

	delete[] vec_y;
	return true;
}


template <class T>
bool Solver<T>::solverJacobi(int iterations, double allowed_convergence, T init_guess, CSRMatrix<T>& a_left, T* b_right_value, T* x_ans, int& counter) {
	// Assign all solutions for the initial guess
	for (int i = 0; i < a_left.rows; i++) {
		x_ans[i] = init_guess;
	}

	auto* diag_mat = new CSRMatrix<T>(a_left.rows, a_left.cols, a_left.rows, true);
	auto* remain_mat = new CSRMatrix<T>(a_left.rows, a_left.cols, a_left.nnzs, true);

	a_left.getMatDiag(*diag_mat);
	a_left.getMatRemn(*remain_mat);
	// Get the inverse of diagonal matrix
	for (int i = 0; i < a_left.rows; i++) {
		if (a_left.getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value"
				<< "\n No solution";
			return false;
		}
		else {
			diag_mat->setValue(i, i, 1. / diag_mat->getValue(i, i));
		}
	}

	// Result of inverse of diagonal matrix dot multiply (vector) b
	// output_db = inv(D) .* b
	T* output_db = new T[diag_mat->rows];
	diag_mat->matVecMult(diag_mat->rows, b_right_value, output_db);
	// Get negative of inverse of diagonal matrix for later use
	for (int i = 0; i < diag_mat->rows; i++) {
		if (diag_mat->getValue(i, i) != 0) {
			diag_mat->setValue(i, i, (-diag_mat->getValue(i, i)));
		}
	}

	counter = 0;
	do {

		// Result of negative inverse of diagonal matrix dot multiple remainder matrix
		// output_dr = -inv(D) .* R
		auto* output_dr = new CSRMatrix<T>(diag_mat->rows, remain_mat->cols, a_left.nnzs, true);
		diag_mat->matMatMult(*remain_mat, *output_dr);
		// Result of output_dr matrix dot multiply (vector) x
		// output_drx = output_dr .* x
		T* output_drx = new T[output_dr->rows];

		output_dr->matVecMult(a_left.rows, x_ans, output_drx);

		// Result of output_drx add output_db
		// output_add = output_drx + output_db
		T* output_add = new T[output_dr->rows];
		for (int i = 0; i < output_dr->rows; i++) {
			output_add[i] = output_drx[i] + output_db[i];
		}

		T err = 0;
		int err_counter = 0;
		for (int i = 0; i < a_left.rows; i++) {
			err = abs(output_add[i] - x_ans[i]);
			x_ans[i] = output_add[i];
			if (err <= allowed_convergence) {
				err_counter++;
			}
		}

		counter++;
		if (err_counter == a_left.rows) {
			delete diag_mat;
			delete remain_mat;
			delete output_dr;
			delete[] output_drx;
			delete[] output_db;
			delete[] output_add;

			return true;
		}

		delete output_dr;
		delete[] output_drx;
		delete[] output_add;

	} while (counter < iterations);

	delete[] output_db;
	delete diag_mat;
	delete remain_mat;

	return false;
}

template <class T>
bool Solver<T>::solverGausSeid(int iterations, double allowed_convergence, T init_guess, CSRMatrix<T>& a_left, T* b_right_value, T* x_ans, int& counter) {
	// Assign all solutions for the initial guess
	T* temp_x = new T[a_left.rows];
	for (int i = 0; i < a_left.rows; i++) {
		x_ans[i] = init_guess;
		temp_x[i] = 0;
	}
	auto* diag_mat = new CSRMatrix<T>(a_left.rows, a_left.cols, a_left.cols, true);
	a_left.getMatDiag(*diag_mat);
	// Get the negative inverse of diagonal matrix
	for (int i = 0; i < diag_mat->rows; i++) {
		if (diag_mat->getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value"
				<< "\n No solution";
			return false;
		}
		else {
			diag_mat->setValue(i, i, 1. / diag_mat->getValue(i, i));
		}
	}

	T temp_L, temp_U, temp_sum, err;
	counter = 0;
	int err_counter = 0; // if err is less than allowed convergence, err_counter + 1
	do {
		for (int i = 0; i < a_left.rows; i++) {

			temp_L = temp_U = temp_sum = err = 0;

			for (int j = 0; j < i; j++) {
				temp_L += a_left.getValue(i, j) * x_ans[j];

			}
			for (int k = a_left.cols - 1; k >= i + 1; k--) {
				temp_U += a_left.getValue(i, k) * x_ans[k];
			}
			temp_sum = b_right_value[i] - temp_L - temp_U;
			temp_x[i] = diag_mat->getValue(i, i) * temp_sum;

			err = abs(temp_x[i] - x_ans[i]);
			x_ans[i] = temp_x[i];

			if (err <= allowed_convergence) {
				err_counter++;
			}
		}

		counter++;
		if (err_counter == a_left.rows) {
			delete diag_mat;
			return true;
		}

	} while (counter < iterations);

	delete diag_mat;
	return false;
}

template <class T>
bool Solver<T>::solverGausElim(CSRMatrix<T>& a_left, T* b_right_value, T* x_ans) {

	for (int i = 0; i < a_left.rows; i++) {
		if (a_left.getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value" << endl
				<< "No solution" << endl << endl;
			system("pause");
			return false;
		}
	}

	// Elimination
	T temp = 0.0;
	for (int k = 0; k < a_left.rows - 1; k++) {
		for (int i = k + 1; i < a_left.rows; i++) {
			temp = a_left.getValue(i, k) / a_left.getValue(k, k);

			for (int j = k; j < a_left.rows; j++) {
				a_left.setValue(i, j, a_left.getValue(i, j) - temp * a_left.getValue(k, j));
			}
			b_right_value[i] -= temp * b_right_value[k];
		}
	}

	//Back substitution
	for (int k = a_left.rows - 1; k >= 0; k--) {
		temp = 0.0;
		for (int j = k + 1; j < a_left.rows; j++) {
			temp += a_left.getValue(k, j) * x_ans[j];

		}

		x_ans[k] = (b_right_value[k] - temp) / a_left.getValue(k, k);
	}

	return true;
}

template <class T>
bool Solver<T>::solverLuDecom(CSRMatrix<T>& a_left, T* b_right_value, T* x_ans) {

	for (int i = 0; i < a_left.rows; i++) {
		if (a_left.getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value" << endl
				<< "No solution" << endl << endl;
			system("pause");
			return false;
		}
	}

	T* mat_U = new T[a_left.rows * a_left.cols];
	T* mat_L = new T[a_left.rows * a_left.cols];
	T* mat_y = new T[a_left.rows];

	for (int i = 0; i < a_left.rows * a_left.cols; i++) {
		mat_L[i] = 0;
		mat_U[i] = 0;
	}

	for (int i = 0; i < a_left.rows; i++) {
		mat_y[i] = 0;
	}

	for (int i = 0; i < a_left.rows; i++) {
		mat_U[i] = a_left.getValue(0, i);
	}

	for (int i = 1; i < a_left.rows; i++) {
		mat_L[i * a_left.rows] = a_left.getValue(i * a_left.rows) / mat_U[0];
	}

	for (int i = 1; i < a_left.rows; i++) {
		for (int k = i; k < a_left.rows; k++) {
			T sum1 = 0;
			for (int j = 0; j < i; j++) {
				sum1 += mat_L[i * a_left.rows + j] * mat_U[j * a_left.rows + k];
			}
			mat_U[i * a_left.rows + k] = a_left.getValue(i, k) - sum1;
		}
		if (i != a_left.rows - 1) {
			for (int k = i; k < a_left.rows; k++) {
				T sum2 = 0;
				for (int j = 0; j < i; j++) {
					sum2 += mat_L[k * a_left.rows + j] * mat_U[j * a_left.rows + i];
				}
				mat_L[k * a_left.rows + i] = (a_left.getValue(k, i) - sum2) / mat_U[i * a_left.rows + i];
			}
		}
	}

	for (int i = 0; i < a_left.rows; i++) {
		mat_L[i * a_left.rows + i] = 1;
	}

	mat_y[0] = b_right_value[0];
	for (int i = 1; i < a_left.rows; i++) {
		T sum3 = 0;
		for (int k = 0; k < i; k++) {
			sum3 += mat_L[i * a_left.rows + k] * mat_y[k];
		}
		mat_y[i] = b_right_value[i] - sum3;
	}

	x_ans[a_left.rows - 1] = mat_y[a_left.rows - 1] / mat_U[a_left.rows * a_left.rows - 1];
	for (int i = a_left.rows - 2; i >= 0; i--) {
		T sum4 = 0;
		for (int k = i + 1; k < a_left.rows; k++) {
			sum4 += mat_U[i * a_left.rows + k] * x_ans[k];
		}
		x_ans[i] = (mat_y[i] - sum4) / mat_U[i * a_left.rows + i];
	}

	delete[] mat_U;
	delete[] mat_L;
	delete[] mat_y;

	return true;
}
template <class T>
bool Solver<T>::solverSor(int iterations, double allowed_convergence, T init_guess, CSRMatrix<T>& a_left, T* b_right_value, T* x_ans, int& counter) {

	for (int i = 0; i < a_left.rows; i++) {
		if (a_left.getValue(i, i) == 0) {
			cerr << "Main diagonal elements of LHS A matrix exist 0 value" << endl
				<< "No solution" << endl << endl;
			system("pause");
			return false;
		}
	}

	T w = 1.46;
	T* vec_y = new T[a_left.rows];
	for (int i = 0; i < a_left.rows; i++) {
		vec_y[i] = 0;
	}

	for (int i = 0; i < a_left.rows; i++) {
		x_ans[i] = init_guess;
	}

	counter = 0;
	do {
		T err = 0;
		for (int i = 0; i < a_left.rows; i++) {
			T s = 0;
			for (int j = 0; j < a_left.rows; j++)
				if (j != i)
					s += a_left.getValue(i, j) * x_ans[j];
			vec_y[i] = (b_right_value[i] - s) / a_left.getValue(i, i);
			vec_y[i] = (1 - w) * x_ans[i] + w * vec_y[i];
			err = abs(x_ans[i] - vec_y[i]);
			x_ans[i] = vec_y[i];
			if (err <= allowed_convergence) break;
		}
		counter++;
	} while (counter < iterations);

	if (counter == iterations) {
		delete[] vec_y;
		return false;
	}

	delete[] vec_y;
	return true;
}
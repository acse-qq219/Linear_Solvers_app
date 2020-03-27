#include <fstream>
#include <iomanip>
#include "Interface.h"
#include "Matrix.h"
#include "CSRMatrix.h"
#include "Solver.h"
#include "Test.h"

#define MAX_ITERATION	1e5
#define MAX_MAT_ROW		1e6
#define MAX_PRECISION	251

using namespace std;

template <class T>
Interface<T>::Interface() {
	interfaceIntro();
}

template <class T>
Interface<T>::~Interface() {
	delete[] value_mat_main;
	delete[] value_b_vec_main;
}

template <class T>
void Interface<T>::interfaceInvalid(string message) {
	cin.clear();
	cin.ignore(numeric_limits<streamsize>::max(), '\n');
	cerr << endl << endl
		<< "ERROR:\t " << message
		<< endl << endl;
	system("pause");
}

template <class T>
void Interface<T>::interfaceIntro() {

	system("CLS");

	cout << endl << endl << endl
		<< ">\t TEAM FOR JUEJUE" << endl
		<< ">\t LINEAR SOLVER SYSTEM" << endl
		<< ">\t DESIGNED BY: QIUCHEN QIAN, YUSEN ZHOU, HE ZHU"
		<< endl << endl;

	cout << endl << endl << endl
		<< ">\t This software is designed to solve linear systems" << endl
		<< ">\t e.g. Ax = b and we have solvers of" << endl
		<< ">\t Direct: Gaussian Elimination, LU Decomposition " << endl
		<< ">\t Iterative: Jacobi Method, Gauss-Seidel Method"
		<< endl << endl;

	cout << endl << "NOTICE:\t A is a n x n square matrix !" << endl
		<< "NOTICE:\t b is a n size vector !"
		<< endl << endl;

	cout << endl << ">\t Would you like to continue?" << endl
		<< ">\t a: SELF-TEST-OF-ALL-SOLVERS" << endl
		<< ">\t y: Select input data" << endl
		<< ">\t x: exit" << endl << endl;
	string select;
	cin >> select;

	char select_char = select[0];

	if (select.size() != 1) {
		select_char = '#';
	}

	switch (select_char)
	{
	case 'a': system("CLS"); interfaceComTest(); break;
	case 'y': break;
	case 'n': exit(0);
	default: interfaceInvalid("Invalid Selection !"); interfaceIntro();
	}

	interfaceSelectData();
}

template <class T>
void Interface<T>::interfaceComTest() {
	Test test;
	test.Performance_Test_all();

	system("pause");

	interfaceIntro();
}

template <class T>
void Interface<T>::interfaceSelectData() {

	system("CLS");

	cout << endl << ">\t Choose how to generate the LHS matrix A:"
		<< endl << endl
		<< "\t 1: Generate a matrix with random number" << endl
		<< "\t 2: Read from data file" << endl
		<< "\t b: Back" << endl
		<< "\t x: Exit"
		<< endl << endl
		<< "INPUT:\t";
	string select;
	cin >> select;

	char select_char = select[0];

	if (select.size() != 1) {
		select_char = '#';
	}

	switch (select_char)
	{
	case '1': generateRandMat(); break;
	case '2': readDataFromFile(); break;
	case 'b': interfaceIntro(); break;
	case 'x': exit(0);
	default: interfaceInvalid("Invalid Selection !"); interfaceSelectData();
	}
}

template <class T>
void Interface<T>::generateRandMat() {

	system("CLS");

	cout << endl << ">\t Enter row number of the LHS matrix A: (it is the column number as well)"
		<< endl << "EXAMPLE:\t 15"
		<< endl << "NOTICE:\t minimum size: 2x2, maximum size: 100000x100000"
		<< endl << endl
		<< "INPUT:\t ";
	int row;
	cin >> row;
	if (cin.fail()) {
		interfaceInvalid("Please Enter an integer only !");
		generateRandMat();
	}
	if (row <= 1 || row >= MAX_MAT_ROW) {
		interfaceInvalid("Out of range, minimum size: 2x2, maximum size: 100000x100000");
		generateRandMat();
	}

	this->row_desmat_main = row;
	this->col_desmat_main = row;

	cout << endl << endl
		<< ">\t Enter the minimum value of random number:"
		<< endl << "EXAMPLE:\t 2.333"
		<< endl << endl
		<< "INPUT:\t ";
	double minimum;
	cin >> minimum;
	while (cin.fail()) {
		cerr << endl << "ERROR:\t Please enter a number !" << endl << endl;
		cout << "INPUT:\t ";
		cin.clear();
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin >> minimum;
	}

	cout << endl << endl
		<< ">\t Enter the maximum value of random number:"
		<< endl << "EXAMPLE:\t 9.999"
		<< endl << endl
		<< "INPUT:\t ";
	double maximum;
	cin >> maximum;
	while (cin.fail() || maximum <= minimum) {
		cerr << endl << "ERROR:\t Please enter a larger number !" << endl << endl;
		cout << "INPUT:\t ";
		cin.clear();
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin >> maximum;
	}

	this->value_mat_main = new T[row * row];

	cout << endl << endl
		<< ">\t Generating..."
		<< endl << endl;
	srand((unsigned int)time(NULL));

	int zero_counter = 0;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < row; j++) {
			T temp = (T)rand() / RAND_MAX;
			this->value_mat_main[i * row + j] = minimum + temp * (maximum - minimum);
			if (this->value_mat_main[i * row + j] == 0) {
				zero_counter++;
			}
			if (this->value_mat_main[i * row + i] <= (maximum - minimum) / 2) {
				this->value_mat_main[i * row + i] += (maximum - minimum) / 1.8;
			}
		}
	}

	int non_zero_counter = 0;
	if (zero_counter >= row * row / 2) {
		this->is_csr_mat = true;
		this->nnzs_main = row * row - zero_counter;
		this->value_csr_mat_main = new T[this->nnzs_main];
		this->row_pos_csrmat_main = new int[row + 1];
		this->col_pos_csrmat_main = new int[this->nnzs_main];
		this->row_pos_csrmat_main[0] = 0;
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < row; j++) {
				if (this->value_mat_main[i * row + j] != 0) {
					this->value_csr_mat_main[non_zero_counter] = this->value_mat_main[i * row + j];
					this->col_pos_csrmat_main[non_zero_counter] = j;
					non_zero_counter++;
				}
			}
			this->row_pos_csrmat_main[i + 1] = non_zero_counter;
		}
	}
	else {
		this->is_csr_mat = false;
	}

	cout << ">\t Done."
		<< endl << endl;
	system("pause");

	generateRandVec();
}

template <class T>
void Interface<T>::generateRandVec() {

	system("CLS");

	cout << endl << ">\t The RHS vector b has same size as rows of matrix A."
		<< endl << ">\t Enter the minimum value of random number:"
		<< endl << "EXAMPLE:\t 2.333"
		<< endl << endl
		<< "INPUT:\t ";
	double minimum;
	cin >> minimum;
	if (cin.fail()) {
		interfaceInvalid("Please Enter a number only !");
		generateRandVec();
	}

	cout << endl << endl
		<< ">\t Enter the maximum value of random number:"
		<< endl << "EXAMPLE:\t 9.999"
		<< endl << endl
		<< "INPUT:\t ";
	double maximum;
	cin >> maximum;
	while (cin.fail() || maximum <= minimum) {
		cerr << endl << "ERROR:\t Please enter a larger number !" << endl << endl;
		cout << "INPUT:\t ";
		cin.clear();
		cin.ignore(numeric_limits<streamsize>::max(), '\n');
		cin >> maximum;
	}

	this->value_b_vec_main = new T[this->row_desmat_main];

	cout << endl << endl << ">\t Generating..."
		<< endl << endl;

	srand((unsigned int)time(NULL));

	for (int i = 0; i < this->row_desmat_main; i++) {
		T temp = (T)rand() / RAND_MAX;
		this->value_b_vec_main[i] = minimum + temp * (maximum - minimum);
	}

	cout << ">\t Done."
		<< endl << endl;
	system("pause");

	interfaceSelectSolver();
}

template <class T>
void Interface<T>::readDataFromFile() {

	system("CLS");

	fstream myfile;
	string filename;
	cout << endl << endl
		<< ">\t Please enter the path to access the txt file:" << endl
		<< "EXAMPLE:\t ./data/input.txt"
		<< endl << endl
		<< "INPUT:\t ";
	cin >> filename;

	myfile.open(filename, fstream::in);

	if (!myfile) {
		interfaceInvalid("Please enter a right path !");
		readDataFromFile();
	}

	string input;
	getline(myfile, input);

	cout << endl << endl
		<< ">\t Reading file from   " << filename << "   ..."
		<< endl << endl;

	int rows, cols;
	rows = cols = 0;
	myfile >> rows >> cols;
	if (myfile.fail()) {
		interfaceInvalid("The data format is wrong !");
		myfile.close();
		interfaceSelectData();
	}

	int zero_counter = 0;
	int test = 0;
	this->value_mat_main = new T[rows * cols];
	this->value_b_vec_main = new T[rows];

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			myfile >> this->value_mat_main[i * cols + j];

			if (myfile.fail()) {
				interfaceInvalid("The data format for matrix is wrong !");
				myfile.close();
				interfaceSelectData();
			}
			if (this->value_mat_main[i * cols + j] == 0) {
				zero_counter++;
			}
		}
	}

	for (int i = 0; i < rows; i++) {
		myfile >> this->value_b_vec_main[i];
		if (myfile.fail()) {
			interfaceInvalid("The data format for vector is wrong !");
			myfile.close();
			interfaceSelectData();
		}
	}

	myfile.close();

	int non_zero_counter = 0;
	if (zero_counter >= rows * cols / 2) {
		this->is_csr_mat = true;
		this->nnzs_main = rows * cols - zero_counter;
		this->value_csr_mat_main = new T[this->nnzs_main];
		this->row_pos_csrmat_main = new int[rows + 1];
		this->col_pos_csrmat_main = new int[this->nnzs_main];
		this->row_pos_csrmat_main[0] = 0;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < rows; j++) {
				if (this->value_mat_main[i * rows + j] != 0) {
					this->value_csr_mat_main[non_zero_counter] = this->value_mat_main[i * rows + j];
					this->col_pos_csrmat_main[non_zero_counter] = j;
					non_zero_counter++;
				}
			}
			this->row_pos_csrmat_main[i + 1] = non_zero_counter;
		}
	}
	else {
		this->is_csr_mat = false;
	}

	this->row_desmat_main = rows;
	this->col_desmat_main = cols;

	cout << ">\t Done."
		<< endl << endl;
	system("pause");

	interfaceSelectSolver();
}

template <class T>
void Interface<T>::interfaceSelectSolver() {

	system("CLS");

	cout << endl << ">\t Choose how to solve this linear system:"
		<< endl << endl
		<< ">\t 1: Jacobi Method" << endl
		<< ">\t 2: Gauss-Seidel Method" << endl
		<< ">\t 3: Gaussian Elimination Method" << endl
		<< ">\t 4: LU Decomposition Method" << endl
		<< ">\t 5: SOR Method" << endl
		<< ">\t b: back" << endl
		<< ">\t x: exit" << endl
		<< endl << endl
		<< "INPUT:\t ";
	string select;
	cin >> select;

	char select_char = select[0];

	if (select.size() != 1) {
		select_char = '#';
	}

	switch (select_char)
	{
	case '1': system("CLS"); interfaceJacobi(); break;
	case '2': system("CLS"); interfaceGasSed(); break;
	case '3': system("CLS"); interfaceGasElim(); break;
	case '4': system("CLS"); interfaceLUDecom(); break;
	case '5': system("CLS"); interfaceSOR(); break;
	case 'b':
		if (!this->have_output_file) {
			delete[] this->value_b_vec_main;
			delete[] this->value_mat_main;
			delete[] this->value_x_vec_main;
		}
		interfaceSelectData(); break;
	case 'x': exit(0);
	default: interfaceInvalid("Invalid Selection !"); interfaceSelectSolver();
	}
}

template <class T>
void Interface<T>::interfaceJacobi() {

	this->have_output_file = false;
	this->start_time = this->end_time = 0;
	this->iter_times = 1;
	auto* solver_main = new Solver<T>();

	if (this->is_csr_mat) {
		auto* a_main = new CSRMatrix<T>(this->row_desmat_main, this->col_desmat_main, this->nnzs_main,
			this->value_csr_mat_main, this->row_pos_csrmat_main, this->col_pos_csrmat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl <<
				">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << "<\t Please enter times you want to iterate: (1 - 10000)"
			<< endl << "EXAMPLE:\t 500"
			<< endl << endl
			<< "INPUT:\t ";
		int iter;
		cin >> iter;

		if (cin.fail()) {
			interfaceInvalid("Please Enter an integer only !");
			interfaceJacobi();
		}
		if (iter == 0 || iter >= MAX_ITERATION) {
			interfaceInvalid("Out of range, please enter an integer from 1 to 10000 only !");
			interfaceJacobi();
		}

		cout << endl << endl
			<< "Please enter a non-zero initial guess for all solutions:"
			<< endl << "EXAMPLE:\t 2.333"
			<< endl << endl
			<< "INPUT:\t ";
		double init;
		cin >> init;
		while (cin.fail() || init == 0.0) {
			cerr << endl << "ERROR:\t Please enter a non-zero number !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> init;
		}

		cout << endl << endl
			<< "Please enter an integer for precision (to decimal places) (1 - 200):"
			<< endl << "EXAMPLE:\t 5"
			<< endl << endl
			<< "INPUT:\t ";
		int user_prec;
		cin >> user_prec;
		while (cin.fail() || user_prec <= 0 || user_prec >= MAX_PRECISION) {
			cerr << endl << "ERROR:\t Please enter a positive integer with correct range !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> user_prec;
		}

		this->prec = user_prec;
		double exact_prec = pow(10, -user_prec);

		cout << endl << endl
			<< ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverJacobi(iter, exact_prec, init, *a_main,
			this->value_b_vec_main, this->value_x_vec_main, this->iter_times);
		this->end_time = clock();

		delete a_main;
		delete[] this->row_pos_csrmat_main;
		delete[] this->col_pos_csrmat_main;
		delete[] this->value_csr_mat_main;
	}
	else {
		auto* a_main = new Matrix<T>(this->row_desmat_main, this->col_desmat_main, this->value_mat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl <<
				">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << "<\t Please enter times you want to iterate: (1 - 10000)"
			<< endl << "EXAMPLE:\t 500"
			<< endl << endl
			<< "INPUT:\t ";
		int iter;
		cin >> iter;

		if (cin.fail()) {
			interfaceInvalid("Please Enter an integer only !");
			interfaceJacobi();
		}
		if (iter == 0 || iter >= MAX_ITERATION) {
			interfaceInvalid("Out of range, please enter an integer from 1 to 10000 only !");
			interfaceJacobi();
		}

		cout << endl << endl
			<< "Please enter a non-zero initial guess for all solutions:"
			<< endl << "EXAMPLE:\t 2.333"
			<< endl << endl
			<< "INPUT:\t ";
		double init;
		cin >> init;
		while (cin.fail() || init == 0.0) {
			cerr << endl << "ERROR:\t Please enter a non-zero number !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> init;
		}

		cout << endl << endl
			<< "Please enter an integer for precision (to decimal places) (1 - 200):"
			<< endl << "EXAMPLE:\t 5"
			<< endl << endl
			<< "INPUT:\t ";
		int user_prec;
		cin >> user_prec;
		while (cin.fail() || user_prec <= 0 || user_prec >= MAX_PRECISION) {
			cerr << endl << "ERROR:\t Please enter a positive integer with correct range !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> user_prec;
		}

		this->prec = user_prec;
		double exact_prec = pow(10, -user_prec);

		cout << endl << endl
			<< ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverJacobi(iter, exact_prec, init, *a_main,
			this->value_b_vec_main, this->value_x_vec_main, this->iter_times);
		this->end_time = clock();

		delete a_main;
	}

	delete solver_main;

	cout << ">\t Done."
		<< endl << endl;
	system("pause");

	interfaceSelectOutput();
}

template <class T>
void Interface<T>::interfaceGasSed() {

	this->have_output_file = false;
	this->start_time = this->end_time = 0;
	this->iter_times = 1;
	auto* solver_main = new Solver<T>();

	if (this->is_csr_mat) {
		auto* a_main = new CSRMatrix<T>(this->row_desmat_main, this->col_desmat_main, this->nnzs_main,
			this->value_csr_mat_main, this->row_pos_csrmat_main, this->col_pos_csrmat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl << ">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << "<\t Please enter times you want to iterate: (1 - 10000)"
			<< endl << "EXAMPLE:\t 2333"
			<< endl << endl
			<< "INPUT:\t ";
		int iter;
		cin >> iter;

		if (cin.fail()) {
			interfaceInvalid("Please Enter an integer only !");
			interfaceJacobi();
		}
		if (iter == 0 || iter > MAX_ITERATION) {
			interfaceInvalid("Out of range, please enter an integer from 1 to 10000 only !");
			interfaceJacobi();
		}

		cout << endl << endl
			<< "Please enter a non-zero initial guess for all solutions:"
			<< endl << "EXAMPLE:\t 2.333"
			<< endl << endl
			<< "INPUT:\t ";
		double init;
		cin >> init;
		while (cin.fail() || init == 0.0) {
			cerr << endl << "Please enter a non-zero number !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> init;
		}

		cout << endl << endl
			<< "Please enter an integer for precision (to decimal places):"
			<< endl << "EXAMPLE:\t 5"
			<< endl << endl
			<< "INPUT:\t ";
		int user_prec;
		cin >> user_prec;
		while (cin.fail() || user_prec <= 0 || user_prec >= MAX_PRECISION) {
			cerr << endl << "ERROR:\t Please enter a positive integer with correct range !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> user_prec;
		}

		this->prec = user_prec;
		double exact_prec = pow(10, -user_prec);

		cout << endl << endl
			<< ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverGausSeid(iter, exact_prec, init, *a_main,
			this->value_b_vec_main, this->value_x_vec_main, this->iter_times);
		this->end_time = clock();

		delete a_main;
		delete[] this->row_pos_csrmat_main;
		delete[] this->col_pos_csrmat_main;
		delete[] this->value_csr_mat_main;
	}
	else {
		auto* a_main = new Matrix<T>(this->row_desmat_main, this->col_desmat_main, this->value_mat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl << ">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << "<\t Please enter times you want to iterate: (1 - 10000)"
			<< endl << "EXAMPLE:\t 2333"
			<< endl << endl
			<< "INPUT:\t ";
		int iter;
		cin >> iter;

		if (cin.fail()) {
			interfaceInvalid("Please Enter an integer only !");
			interfaceJacobi();
		}
		if (iter == 0 || iter > MAX_ITERATION) {
			interfaceInvalid("Out of range, please enter an integer from 1 to 10000 only !");
			interfaceJacobi();
		}

		cout << endl << endl
			<< "Please enter a non-zero initial guess for all solutions:"
			<< endl << "EXAMPLE:\t 2.333"
			<< endl << endl
			<< "INPUT:\t ";
		double init;
		cin >> init;
		while (cin.fail() || init == 0.0) {
			cerr << endl << "Please enter a non-zero number !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> init;
		}

		cout << endl << endl
			<< "Please enter an integer for precision (to decimal places):"
			<< endl << "EXAMPLE:\t 5"
			<< endl << endl
			<< "INPUT:\t ";
		int user_prec;
		cin >> user_prec;
		while (cin.fail() || user_prec <= 0 || user_prec >= MAX_PRECISION) {
			cerr << endl << "ERROR:\t Please enter a positive integer with correct range !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> user_prec;
		}

		this->prec = user_prec;
		double exact_prec = pow(10, -user_prec);

		cout << endl << endl
			<< ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverGausSeid(iter, exact_prec, init, *a_main,
			this->value_b_vec_main, this->value_x_vec_main, this->iter_times);
		this->end_time = clock();

		delete a_main;
	}

	delete solver_main;

	cout << ">\t Done."
		<< endl << endl;
	system("pause");

	interfaceSelectOutput();
}

template <class T>
void Interface<T>::interfaceSOR() {

	this->have_output_file = false;
	this->start_time = this->end_time = 0;
	this->iter_times = 1;
	auto* solver_main = new Solver<T>();

	if (this->is_csr_mat) {
		auto* a_main = new CSRMatrix<T>(this->row_desmat_main, this->col_desmat_main, this->nnzs_main,
			this->value_csr_mat_main, this->row_pos_csrmat_main, this->col_pos_csrmat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl << ">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << "<\t Please enter times you want to iterate: (1 - 10000)"
			<< endl << "EXAMPLE:\t 2333"
			<< endl << endl
			<< "INPUT:\t ";
		int iter;
		cin >> iter;

		if (cin.fail()) {
			interfaceInvalid("Please Enter an integer only !");
			interfaceJacobi();
		}
		if (iter == 0 || iter > MAX_ITERATION) {
			interfaceInvalid("Out of range, please enter an integer from 1 to 10000 only !");
			interfaceJacobi();
		}

		cout << endl << endl
			<< "Please enter a non-zero initial guess for all solutions:"
			<< endl << "EXAMPLE:\t 2.333"
			<< endl << endl
			<< "INPUT:\t ";
		double init;
		cin >> init;
		while (cin.fail() || init == 0.0) {
			cerr << endl << "Please enter a non-zero number !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> init;
		}

		cout << endl << endl
			<< "Please enter an integer for precision (to decimal places):"
			<< endl << "EXAMPLE:\t 5"
			<< endl << endl
			<< "INPUT:\t ";
		int user_prec;
		cin >> user_prec;
		while (cin.fail() || user_prec <= 0 || user_prec >= MAX_PRECISION) {
			cerr << endl << "ERROR:\t Please enter a positive integer with correct range !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> user_prec;
		}

		this->prec = user_prec;
		double exact_prec = pow(10, -user_prec);

		cout << endl << endl
			<< ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverSor(iter, exact_prec, init, *a_main,
			this->value_b_vec_main, this->value_x_vec_main, this->iter_times);
		this->end_time = clock();

		delete a_main;
		delete[] this->row_pos_csrmat_main;
		delete[] this->col_pos_csrmat_main;
		delete[] this->value_csr_mat_main;
	}
	else {
		auto* a_main = new Matrix<T>(this->row_desmat_main, this->col_desmat_main, this->value_mat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl << ">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << "<\t Please enter times you want to iterate: (1 - 10000)"
			<< endl << "EXAMPLE:\t 2333"
			<< endl << endl
			<< "INPUT:\t ";
		int iter;
		cin >> iter;

		if (cin.fail()) {
			interfaceInvalid("Please Enter an integer only !");
			interfaceJacobi();
		}
		if (iter == 0 || iter > MAX_ITERATION) {
			interfaceInvalid("Out of range, please enter an integer from 1 to 10000 only !");
			interfaceJacobi();
		}

		cout << endl << endl
			<< "Please enter a non-zero initial guess for all solutions:"
			<< endl << "EXAMPLE:\t 2.333"
			<< endl << endl
			<< "INPUT:\t ";
		double init;
		cin >> init;
		while (cin.fail() || init == 0.0) {
			cerr << endl << "Please enter a non-zero number !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> init;
		}

		cout << endl << endl
			<< "Please enter an integer for precision (to decimal places):"
			<< endl << "EXAMPLE:\t 5"
			<< endl << endl
			<< "INPUT:\t ";
		int user_prec;
		cin >> user_prec;
		while (cin.fail() || user_prec <= 0 || user_prec >= MAX_PRECISION) {
			cerr << endl << "ERROR:\t Please enter a positive integer with correct range !"
				<< endl << endl;
			cout << "INPUT:\t ";
			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cin >> user_prec;
		}

		this->prec = user_prec;
		double exact_prec = pow(10, -user_prec);

		cout << endl << endl
			<< ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverSor(iter, exact_prec, init, *a_main,
			this->value_b_vec_main, this->value_x_vec_main, this->iter_times);
		this->end_time = clock();

		delete a_main;
	}

	delete solver_main;

	cout << ">\t Done."
		<< endl << endl;
	system("pause");

	interfaceSelectOutput();
}

template <class T>
void Interface<T>::interfaceGasElim() {

	this->have_output_file = false;
	this->start_time = this->end_time = 0;
	this->iter_times = 1;
	auto* solver_main = new Solver<T>();

	if (this->is_csr_mat) {
		auto* a_main = new CSRMatrix<T>(this->row_desmat_main, this->col_desmat_main, this->nnzs_main,
			this->value_csr_mat_main, this->row_pos_csrmat_main, this->col_pos_csrmat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl
				<< ">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverGausElim(*a_main, this->value_b_vec_main, this->value_x_vec_main);
		this->end_time = clock();

		delete a_main;
		delete[] this->row_pos_csrmat_main;
		delete[] this->col_pos_csrmat_main;
		delete[] this->value_csr_mat_main;
	}
	else {
		auto* a_main = new Matrix<T>(this->row_desmat_main, this->col_desmat_main, this->value_mat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl
				<< ">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverGausElim(*a_main, this->value_b_vec_main, this->value_x_vec_main);
		this->end_time = clock();

		delete a_main;
	}

	delete solver_main;

	cout << ">\t Done."
		<< endl << endl;
	system("pause");

	interfaceSelectOutput();
}

template <class T>
void Interface<T>::interfaceLUDecom() {

	this->have_output_file = false;
	this->start_time = this->end_time = 0;
	this->iter_times = 1;
	auto* solver_main = new Solver<T>();

	if (this->is_csr_mat) {
		auto* a_main = new CSRMatrix<T>(this->row_desmat_main, this->col_desmat_main, this->nnzs_main,
			this->value_csr_mat_main, this->row_pos_csrmat_main, this->col_pos_csrmat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl
				<< ">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverLuDecom(*a_main, this->value_b_vec_main, this->value_x_vec_main);
		this->end_time = clock();

		delete a_main;
		delete[] this->row_pos_csrmat_main;
		delete[] this->col_pos_csrmat_main;
		delete[] this->value_csr_mat_main;
	}
	else {
		auto* a_main = new Matrix<T>(this->row_desmat_main, this->col_desmat_main, this->value_mat_main);

		if (this->row_desmat_main <= 30) {
			a_main->printMatrix();
			cout << endl << endl
				<< ">\t The vector b is:" << endl
				<< " [ ";
			for (int i = 0; i < this->row_desmat_main; i++) {
				cout << this->value_b_vec_main[i] << " ";
			}
			cout << "]" << endl << endl;
		}
		else {
			cout << endl << endl
				<< ">\t The Matrix and vector is too large to show !"
				<< endl << endl;
		}

		this->value_x_vec_main = new T[this->row_desmat_main];
		for (int i = 0; i < this->row_desmat_main; i++) this->value_x_vec_main[i] = 0;

		cout << ">\t Processing..."
			<< endl << endl;

		this->start_time = clock();
		this->is_solver_right = solver_main->solverLuDecom(*a_main, this->value_b_vec_main, this->value_x_vec_main);
		this->end_time = clock();

		delete a_main;
	}

	delete solver_main;

	cout << ">\t Done."
		<< endl << endl;
	system("pause");

	interfaceSelectOutput();
}

template <class T>
void Interface<T>::interfaceSelectOutput() {

	system("CLS");

	cout << endl << ">\t Choose how to deal with data:"
		<< endl << endl
		<< ">\t 1: print solutions and running time of the solver" << endl
		<< ">\t 2: output data to a CSV file" << endl
		<< ">\t b: back" << endl
		<< ">\t x: exit" << endl
		<< endl << endl
		<< "INPUT:\t ";
	string select;
	cin >> select;

	char select_char = select[0];

	if (select.size() != 1) {
		select_char = '#';
	}

	switch (select_char)
	{
	case '1': printSolutionVec(); break;
	case '2': saveDataToFile(); break;
	case 'b': interfaceSelectData(); break;
	case 'x': exit(0);
	default: interfaceInvalid("Invalid selection !"); interfaceSelectOutput(); break;
	}
}

template<class T>
void Interface<T>::printSolutionVec() {

	if (!this->is_solver_right) {
		interfaceInvalid("This linear system has no acceptable solutions ! Try solver 3 or 4 please !");
		interfaceSelectOutput();
	}

	if (this->prec <= 15) {
		cout << fixed << setprecision(this->prec);
	}
	else {
		cout << endl << endl
			<< "NOTICE:\t The solution is too long to print, only display first 15 places";
		cout << fixed << setprecision(15);
	}

	cout << endl << endl
		<< ">\t The solution is:" << endl
		<< " [ ";
	for (int i = 0; i < this->row_desmat_main; i++) cout << this->value_x_vec_main[i] << " ";
	cout << "]";

	cout << endl << endl
		<< ">\t End at the " << this->iter_times << "-th iteration.";

	cout << endl << endl
		<< ">\t Running time of the solver is: "
		<< fixed << setprecision(6)
		<< (double)(this->end_time - this->start_time) / (double)(CLOCKS_PER_SEC) * 1000
		<< " ms" << endl << endl;

	system("pause");

	interfaceSelectOutput();
}

template<class T>
void Interface<T>::saveDataToFile() {

	system("CLS");

	if (this->have_output_file) {
		interfaceInvalid("Do not output repeated file !");
		interfaceSelectOutput();
	}

	fstream myFile;
	string path = "./data/output.txt";
	cout << endl << ">\t Please enter the path you want to save:"
		<< endl << "EXAMPLE:\t ./data/output.txt"
		<< endl << endl
		<< "INPUT:\t ";
	cin >> path;
	myFile.open(path, fstream::out);
	if (myFile.fail()) {
		interfaceInvalid("Cannot open the file !");
		interfaceSelectOutput();
	}

	cout << endl << endl << "Saving to " << path << " ...";

	myFile << "matrix A" << setprecision(this->prec);
	for (int i = 0; i < this->row_desmat_main; i++) {
		myFile << endl;
		for (int j = 0; j < this->col_desmat_main; j++) {
			myFile << this->value_mat_main[i * this->row_desmat_main + j] << " ";
		}
	}

	myFile << endl << endl << "vector b" << endl;
	for (int i = 0; i < this->row_desmat_main; i++) {
		myFile << this->value_b_vec_main[i] << " ";
	}

	myFile << endl << endl << "solution" << endl;
	for (int i = 0; i < this->row_desmat_main; i++) {
		myFile << this->value_x_vec_main[i] << " ";
	}

	myFile << endl << endl << "end in " << this->iter_times << "-th iteration";

	myFile << endl << endl << "running time: "
		<< setprecision(6)
		<< (double)(this->end_time - this->start_time) / (double)(CLOCKS_PER_SEC) * 1000 << " ms";

	myFile.close();
	delete[] this->value_mat_main;
	delete[] this->value_b_vec_main;
	delete[] this->value_x_vec_main;

	this->have_output_file = true;

	cout << endl << endl
		<< ">\t Completed."
		<< endl << endl;
	system("pause");

	interfaceSelectData();
}
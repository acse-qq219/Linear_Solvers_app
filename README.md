# Linear-Solver-System

A linear solver project for solving linear equations (like Ax = b) provided by users.

## Getting Started

* Dependencies:  At least C++ MSVC 2017

* Download:

1. Create a directory

```
  mkdir linear-solver
````

2. Enter this directory

```
  cd linear-solver/
```

2. Clone this repository

```
  git clone https://github.com/acse-2019/acse-5-assignment-acse-5-team-for-juejue.git
```

* Deployment: 

1. Open your MSVC and open this local folder
2. Double-click the `linear_solver.sln` in your Solution Explorer (find in "view" bar or `Ctrl + Alt + L`)
3. Build the solution or `Ctrl + Alt + F7`
4. Run `Ctrl + F5`

## User Guidance

* The main process looks like: Introducton (or self-test) -> Select Input Data Source -> Select Solevr -> Select Output

* Users can check the documentation of this project [index](html_Documentation/index.html) to find detail information of all classes and theirs members:

1. Open `/html_Documentation` folder
2. Find `index.html` and double click
3. Check usage of all classes and member functions as you wish   :)

* Users can have a overview of five solvers performance by selecting option "SELF-TEST-OF-ALL-SOLVERS" in the introduction page of the user interface. This will run a program which may take few minutes. Running time of five solvers with different test sets will be displayed to users. Note that test sets are generated randmomly with minimum value 3 and maximum value 100.

* Users may select to generate a specified size matrix with random number between wanted minimum and maximum value
or read data from a `.txt` file. Note that if users want to **solve the linear system through sparsely stored matrix**, they
have to read data source from a external file. 


* The default input data file is stored in `linear-solver\acse-5-assignment-acse-5-team-for-juejue\linear_solvers\data`
as [input.txt](linear_solvers/data/input.txt). The default format of `input.txt` is:

```
     First line: Description/Source/Reference/Author/etc.
     Second line: 10 10      (rows and columns the matrix owns)
     Third line: data data data data data data data data data data (data of the matrix)
                                  .
                                  .
                                  .
     13th line: data data data data data data data data data data (data of the vector)
```

* Users can select their own data file as data source but the data format **MUST BE THE SAME AS ABOVE**. You may find an example file called [input.txt](linear_solvers/data/input.txt) in `\data` folder

* Users can select four linear solvers to solve the linear system:

1. Iterative Jacobi Method -- Users need to input iteration times, initial guess for all solutions and wanted precision
2. Iterative Gauss-Seidel Method -- Users need to input iteration times, initial guess for all solutions and wanted precision
3. Direct Gaussian Elimination Method -- No need for extra input
4. Direct LU Decomposition -- No need for extra input
5. Iterative SOR Method -- Users need to input iteration times, initial guess for all solutions and wanted precision

* Users can select three output:

1. Print all solutions as x vector with wanted precision and solver running time. Note that if no solution or the solution cannot meet precision requirement, the system will not output the result, but users can still choose to output results to a `.txt` file (default [output.txt](linear_solvers/data/output.txt) in `\data` folder). Moreover, the system may show 0 ms sometimes. This is due to small size of the matrix (usually smaller than 50 x 50).
2. Output all data to a txt file. Users can save the output data file (default in `\data` folder).

* An example of format of `output.txt` is: (10 x 10 matrix)

```
     First line: matrix A (string)
     Second line: data data data data data data data data data data (data of the matrix)
                                  .
                                  .
                                  .
     13th line: vector b (string)
     14th line: data data data data data data data data data data (data of the vector b)
     15th line: (whitespace)
     16th line: solution (string)
     17th line: data data data data data data data data data data (data of the vector x)
     18th line: (whitespace)
     19th line: end in data-th iteration
     20th line: (whitespace)
     21th line: running time: data ms
```

## Authors

* Qiuchen Qian -- Interface, DESMatrix, Jacobi and Gauss-Seidel Solver
* Yusen Zhou -- CSRMatrix, test and Intergation
* He Zhu -- Readfile (Interface), Gauss Elimination, LU Decomposition and SOR solver

## License

This project is licensed under the MIT License -- see the [LICENSE.md](LICENSE) file for details

## Acknowledgements

* JueJue

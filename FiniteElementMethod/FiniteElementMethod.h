#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/Householder>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

namespace FiniteElementMethod{
double F(double x);

double ExactSolution(double);

pair<MatrixXd,MatrixXd> SolveGlobalMatrix(const size_t);

void Solution(MatrixXd&,size_t);

double ApproximateSolution(double x, const VectorXd& U, int numElements, double h);

double ErrorInfinite(MatrixXd&,size_t);

double SimpsonMethod(const vector<double>& ,double);
double ErrorL2(MatrixXd&,size_t);

double CondMatrix(MatrixXd&);

void Time(ofstream& );
}
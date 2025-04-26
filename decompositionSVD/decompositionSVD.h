#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <iomanip>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <utility>
#include <Eigen/Dense>
#include <Eigen/Cholesky>

using namespace std;
using namespace Eigen;

namespace SVD{
double P(double,double);

bool fill_matrix(MatrixXd*, MatrixXd* );

pair<double, VectorXd> powerIteration(MatrixXd A, int maxIterations, double tolerance);

vector<pair<double, VectorXd>> findEigenvaluesAndVectors(MatrixXd A, int maxIterations, double tolerance);
}

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Householder>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <regex>
#include <string>

using namespace std;
using namespace Eigen;



double f(double, vector<double>, size_t); //правая часть
double w(double);                       //точное решение

MatrixXd krylovAcceleration(const vector<MatrixXd> &, size_t);

void FillLocalMatrix(MatrixXd &, MatrixXd, MatrixXd &, size_t);
void CollocationEquations(MatrixXd &, MatrixXd &, size_t);
void MatchingEquations(MatrixXd &, MatrixXd &, MatrixXd &, size_t);
void BoundaryConditions(MatrixXd &, MatrixXd &, size_t);

double ChebyshevRoot(size_t);
double ChebyshevPolynom(size_t, double);

double derivative_1(size_t, double);
double derivative_2(size_t, double);
double derivative_3(size_t, double);
double derivative_4(size_t, double);

void QR_method(MatrixXd &, MatrixXd &, MatrixXd &);
void Rotation_method(MatrixXd &, MatrixXd &);

void GaussMethodBack(const MatrixXd &, const MatrixXd &, MatrixXd &, bool);
void QR_method(MatrixXd &, MatrixXd &, MatrixXd &);
void qrHouseholder(const MatrixXd &, MatrixXd &, MatrixXd &);

vector<pair<double, double>> SegmentDivision(double, double);
vector<double> SegmentCenter(double, double);
vector<vector<double>> UniformGrid();

pair<double, double> Error(MatrixXd &);
void Solution(MatrixXd &);
void ErrorInPoint(MatrixXd &);

double CondMatrix(MatrixXd &);

void Time(ofstream &);

void TranspositionRows(MatrixXd &, size_t, size_t);
size_t MaxCoffColls(MatrixXd &, size_t);
double Signum(double x);

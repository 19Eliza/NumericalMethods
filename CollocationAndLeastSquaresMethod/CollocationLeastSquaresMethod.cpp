#include"CollocationLeastSquaresMethod.h"

using namespace std;
using namespace Eigen;

static const int N = 4;
static const int k = 20;

static const double l = 1.0;
static const double g1_0 = 0.0, g1_l = 0.0, g2_0 = 0.0, g2_l = 0.0;
static const double h_2 = l / k;
static const double h = h_2 / 2.0;

static const double p_c = pow(h, 4);

static const double pi = 3.141592653589793;

static const double Q = 100;

const string out1 = {"resultKrylov.txt"};
const string out2 = {"Krylov.txt"};
const string out3 = {"result.txt"};

double f(double y, vector<double> x_c, size_t numCell) {
  double x = y * h + x_c[numCell];
  return (x * x * x * x + 14.0 * x * x * x + 49.0 * x * x + 32.0 * x - 12.0) *
         exp(x);
}

double w(double x) {
  return x * x * (1.0 - x) * (1.0 - x) * exp(x);
}

MatrixXd krylovAcceleration(const vector<MatrixXd> &c, size_t num) {

  if (num < 2 || c.size() < num) {
    throw invalid_argument(
        "Insufficient number of matrices in c or invalid num.");
  }

  size_t n = num - 1; // Number of residuals
  vector<MatrixXd> residuals(n);
  for (size_t i = 0; i < n; ++i) residuals[i] = (c[i + 1]) - (c[i]);
  
  MatrixXd b = -residuals[n - 1];
  MatrixXd A(k * (N + 1), n - 1);
  for (size_t i = 0; i < (n - 1); i++) A.col(i) = residuals[i + 1] - residuals[i];
  
  MatrixXd alpha = A.colPivHouseholderQr().solve(b);
  for (int i = 0; i < (n - 1); ++i) (c[n]) += alpha(i) * residuals[i + 1];

  return c[n];
}

void FillLocalMatrix(MatrixXd &A, MatrixXd &b, MatrixXd &C, size_t numberCell) {
  CollocationEquations(A, b, numberCell);
  if (k > 1)
    MatchingEquations(A, b, C, numberCell);
  if (numberCell == 0 || numberCell == (k - 1))
    BoundaryConditions(A, b, numberCell);
}

void CollocationEquations(MatrixXd &A, MatrixXd &b, size_t numberCell) {
  vector<double> x_c = SegmentCenter(0, l);
  
  for (size_t i = 0; i <= N; i++) {
    double y = ChebyshevRoot(i);
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = derivative_4(j, y);
    }
    b(i, 0) = p_c * f(y, x_c, numberCell);
  }
}

void MatchingEquations(MatrixXd &A, MatrixXd &b, MatrixXd &C,
                        size_t numberCell) {
  double left = -1.0;
  double right = 1.0;

  size_t i = N + 1;

  if (numberCell == 0) {
    // j+1
    b(i, 0) = 0.0;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = ChebyshevPolynom(j, right) + derivative_1(j, right);
      b(i, 0) += C(N * (numberCell + 1) + j + numberCell + 1, 0) *
                    (ChebyshevPolynom(j, left) + derivative_1(j, left));
    }
    ++i;
    b(i, 0) = 0.0;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = derivative_2(j, right) + derivative_3(j, right);
      b(i, 0) += C(N * (numberCell + 1) + j + numberCell + 1, 0) *
                    (derivative_2(j, left) + derivative_3(j, left));
    }
    
  }
  if (numberCell != 0 && numberCell != (k - 1)) {
    // j-1
    b(i, 0) = 0.0;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = ChebyshevPolynom(j, left) - derivative_1(j, left);
      b(i, 0) += C(N * (numberCell - 1) + j + numberCell - 1, 0) *
                    (ChebyshevPolynom(j, right) - derivative_1(j, right));
    }
    ++i;
    b(i, 0) = 0.0;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = derivative_2(j, left) - derivative_3(j, left);
      b(i, 0) += C(N * (numberCell - 1) + j + numberCell - 1, 0) *
                    (derivative_2(j, right) - derivative_3(j, right));
    }
    // j+1
    ++i;
    b(i, 0) = 0.0;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = ChebyshevPolynom(j, right) + derivative_1(j, right);
      b(i, 0) += C(N * (numberCell + 1) + j + numberCell + 1, 0) *
                    (ChebyshevPolynom(j, left) + derivative_1(j, left));
    }
    ++i;
    b(i, 0) = 0.0;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = derivative_2(j, right) + derivative_3(j, right);
      b(i, 0) += C(N * (numberCell + 1) + j + numberCell + 1, 0) *
                    (derivative_2(j, left) + derivative_3(j, left));
    }
    
  }
  if (numberCell == (k - 1)) {
    // j-1
    b(i, 0) = 0.0;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = ChebyshevPolynom(j, left) - derivative_1(j, left);
      b(i, 0) += C(N * (numberCell - 1) + j + numberCell - 1, 0) *
                    (ChebyshevPolynom(j, right) - derivative_1(j, right));
    }
    ++i;
    b(i, 0) = 0.0;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = derivative_2(j, left) - derivative_3(j, left);
      b(i, 0) += C(N * (numberCell - 1) + j + numberCell - 1, 0) *
                    (derivative_2(j, right) - derivative_3(j, right));
    }
  
  }
 
}

void BoundaryConditions(MatrixXd &A, MatrixXd &b, size_t numberCell) {
  size_t i = 0;
  if (k > 1)
    i = N + 3;
  else
    i = N + 1;

  double left = -1.0;
  double right = 1.0;

  if (numberCell == 0) {
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = ChebyshevPolynom(j, left);
    }
    b(i, 0) = 0;
    ++i;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = derivative_1(j, left);
    }
    b(i, 0) = 0;
  }
  if (k == 1)
    i++;
  if (numberCell == (k - 1)) {
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = ChebyshevPolynom(j, right);
    }
    b(i, 0) = 0;
    ++i;
    for (size_t j = 0; j <= N; j++) {
      A(i, j) = derivative_1(j, right);
    }
    b(i, 0) = 0;
  }

}

double ChebyshevRoot(size_t i) {
  return cos(((2.0 * (i + 1) - 1.0) / (2.0 * (N + 1))) * pi);
}

double ChebyshevPolynom(size_t n, double y) {

  if (n == 0)
    return 1.0;
  if (n == 1)
    return y;
  return (2.0 * y * ChebyshevPolynom(n - 1, y) - ChebyshevPolynom(n - 2, y));
 
}

double derivative_1(size_t n, double y) {
  if (n == 0)
    return 0.0;
  if (n == 1)
    return 1.0;
  return (2.0 * y * derivative_1(n - 1, y) - derivative_1(n - 2, y) +
          2.0 * ChebyshevPolynom(n - 1, y));
  
}
double derivative_2(size_t n, double y) {
  if (n == 0)
    return 0.0;
  if (n == 1)
    return 0.0;
  return (2.0 * y * derivative_2(n - 1, y) - derivative_2(n - 2, y) +
          4.0 * derivative_1(n - 1, y));

}

double derivative_3(size_t n, double y) {

  if (n == 0)
    return 0.0;
  if (n == 1)
    return 0.0;
  return (2.0 * y * derivative_3(n - 1, y) - derivative_3(n - 2, y) +
          6.0 * derivative_2(n - 1, y));

}

double derivative_4(size_t n, double y) {
  if (n == 0)
    return 0.0;
  if (n == 1)
    return 0.0;
  return (2.0 * y * derivative_4(n - 1, y) - derivative_4(n - 2, y) +
          8.0 * derivative_3(n - 1, y));
}

void GaussMethodBack(const MatrixXd &A,const MatrixXd &b, MatrixXd &y,
                       bool viewMatrix) {

  if (!viewMatrix) {

    for (size_t i = 0; i < A.rows(); i++) {

      double sum = 0;
      for (size_t j = 0; j + 1 <= i; j++) {
        sum += A(i, j) * (y(j));
      }

      y(i) = (1.0 / A(i, i)) * (b(i)-sum);
    }

  } else {

    for (size_t i = A.rows(); i > 0; i--) {

      double sum = 0;
      for (size_t j = i; j <= A.rows() - 1; j++) {
        sum += A(i - 1, j) * (y(j));
      }

      y(i)  = (1.0 / A(i - 1, i - 1)) * (b(i - 1) - sum)+1;
    }

  }
}

void QR_method(MatrixXd &A, MatrixXd &b, MatrixXd &C) {

  MatrixXd A0 = A;

  MatrixXd Q, R;

  qrHouseholder(A0, Q, R);

  R.conservativeResize(A.cols(), A.cols());

  bool R_view = true; // true, если верхнетруегольная матрица слева

  MatrixXd b0 = Q.transpose() * b;

  b0.conservativeResize(A.cols(), 1);

  GaussMethodBack(R, b0, C, R_view);
}

void qrHouseholder(const MatrixXd &A, MatrixXd &Q, MatrixXd &R) {
  int m = A.rows();
  int n = A.cols();

  MatrixXd H = MatrixXd::Identity(m, m); // Q
  R = A;

  for (int k = 0; k < n; ++k) {
    VectorXd x = R.col(k).tail(m - k);

    VectorXd v = x;
    v(0) += (x(0) >= 0 ? 1 : -1) * x.norm();
    v.normalize();

    MatrixXd H_k = MatrixXd::Identity(m - k, m - k) - 2 * v * v.transpose();

    MatrixXd H_full = MatrixXd::Identity(m, m);
    H_full.block(k, k, m - k, m - k) = H_k;
    R = H_full * R;

    H = H * H_full.transpose();
  }

  Q = H.transpose();
}

vector<pair<double, double>> SegmentDivision(double a, double b) {
  vector<pair<double, double>> cellBoundaries(k);
  for (size_t i = 0; i < k; i++) {
    cellBoundaries[i].first = a + i * h_2;
    cellBoundaries[i].second = a + (i + 1) * h_2;
  }
  return cellBoundaries;
}

vector<double> SegmentCenter(double a, double b) {
  vector<double> cellCenter(k);
  vector<pair<double, double>> cellBoundaries = SegmentDivision(0, l);
  for (size_t i = 0; i < k; i++) {
    cellCenter[i] =
        (cellBoundaries[i].first + cellBoundaries[i].second) / 2.0;
  }
  return cellCenter;
}
vector<vector<double>> UniformGrid() {
  vector<pair<double, double>> cell_boundaries = SegmentDivision(0, l);
  vector<vector<double>> x_m(k);
  for (size_t i = 0; i < cell_boundaries.size(); i++) {
    double a = cell_boundaries[i].first;
    double b = cell_boundaries[i].second;
    double h_loc = (b - a) / (Q - 1);
    for (size_t j = 0; j <= (Q - 1); j++) {
      x_m[i].push_back(a + j * h_loc);
    }
  }
  return x_m;
}

pair<double, double> Error(MatrixXd &C) {
  vector<vector<double>> x_m = UniformGrid(); //размер k
  vector<double> x_c = SegmentCenter(0, l);
  vector<double> diff_in_cell(k);
  vector<double> max_w_in_cell(k);
  for (size_t i = 0; i < k; i++) {
    for (size_t j = 0; j < x_m[i].size(); j++) {
      double y = (x_m[i][j] - x_c[i]) / h;
      double v = 0.0;
      for (size_t s = 0; s <= N; s++) {
        v += C((N + 1) * i + s, 0) * ChebyshevPolynom(s, y);
      }
      if (abs(v - w(x_m[i][j])) > diff_in_cell[i])
        diff_in_cell[i] = abs(v - w(x_m[i][j]));

      if (abs(w(x_m[i][j])) > max_w_in_cell[i])
        max_w_in_cell[i] = abs(w(x_m[i][j]));
    }
  }
  auto iter1 = max_element(diff_in_cell.begin(), diff_in_cell.end());
  auto iter2 = max_element(max_w_in_cell.begin(), max_w_in_cell.end());

  double relative_error = (*iter1) / (*iter2);
  double absolute_error = *iter1;
  return {relative_error, absolute_error};
}

void Solution(MatrixXd &C) {
  ofstream f1{out1};
  ofstream f2{out2};
  vector<vector<double>> x_m = UniformGrid(); //размер k
  vector<double> x_c = SegmentCenter(0, l);
  for (size_t i = 0; i < k; i++) {
    for (size_t j = 0; j < x_m[i].size(); j++) {
      double y = (x_m[i][j] - x_c[i]) / h;
      double v = 0.0;
      for (size_t s = 0; s <= N; s++) {
        v += C((N + 1) * i + s, 0) * ChebyshevPolynom(s, y);
      }
      f1 << x_m[i][j] << " " << v << endl;
      f2 << x_m[i][j] << " " << w(x_m[i][j]) << endl;
    }
  }
  f1.close();
  f2.close();
}

void ErrorInPoint(MatrixXd &C) {
  ofstream f1{out1};
  ofstream f2{out2};
  vector<vector<double>> x_m = UniformGrid(); //размер k
  vector<double> x_c = SegmentCenter(0, l);
  for (size_t i = 0; i < k; i++) {
    for (size_t j = 0; j < x_m[i].size(); j++) {
      double y = (x_m[i][j] - x_c[i]) / h;
      double v = 0.0;
      for (size_t s = 0; s <= N; s++) {
        v += C((N + 1) * i + s, 0) * ChebyshevPolynom(s, y);
      }
      f1 << x_m[i][j] << " " << abs(v - w(x_m[i][j])) << endl;
    }
  }
}

double CondMatrix(MatrixXd &A) {
  JacobiSVD<MatrixXd> svd(A);
  double cond_A = svd.singularValues()(0) /
                  svd.singularValues()(svd.singularValues().size() - 1);
  return cond_A;
}

void Time(ofstream &of) {
  of << "Time: " << (double)clock() / CLOCKS_PER_SEC << endl;
}

void Rotation_method(MatrixXd &A, MatrixXd &b){
  size_t i;
  for (i = 0; i < A.cols(); i++) {

    size_t j = MaxCoffColls(A, i);
    TranspositionRows(A, i, j);
    TranspositionRows(b, i, j);

    for (size_t l = i + 1; l < A.rows(); l++) {
      double cos;
      double sin;
      if (A(i, i) * A(i, i) + A(l, i) * A(l, i) == 0) {
        cos = 1.0;
        sin = 0.0;
      } else {
        double sqr =
            (1.0 / sqrt(A(i, i) * A(i, i) + A(l, i) * A(l, i)));
        cos = A(i, i) * sqr;
        sin = -A(l, i) * sqr;
      }

      double p;
      if (abs(sin) < abs(cos)) {
        p = Signum(cos) * sin;
      } else {
        p = Signum(sin) / cos;
      }

      if (abs(p) < 1.0) {
        sin = p;
        cos = sqrt(1.0 - p * p);
      } else {
        cos = 1.0 / p;
        sin = sqrt(1.0 - cos * cos);
      }

      for (size_t k = 0; k < A.cols(); k++) {
        double tmp1 = A(i, k);
        double tmp2 = A(l, k);
        A(i, k) = tmp1 * cos - tmp2 * sin;
        A(l, k) = tmp1 * sin + tmp2 * cos;
      }

      for (size_t k = 0; k < b.cols(); k++) {
        double tmp1 = b(i, k);
        double tmp2 = b(l, k);
        b(i, k) = tmp1 * cos - tmp2 * sin;
        b(l, k) = tmp1 * sin + tmp2 * cos;
      }
    }
  }
}

void TranspositionRows(MatrixXd &A, size_t i, size_t j) {
  size_t k = 0;
  while (k < A.cols()) {
    swap(A(i, k), A(j, k));
    k++;
  }
}

size_t MaxCoffColls(MatrixXd &A, size_t i) {
  size_t k;
  double max_coff = -1.0;
  size_t max_index;
  for (k = i; k < A.rows(); k++) {
    if (abs(A(k, i)) > max_coff) {
      max_coff = A(k, i);
      max_index = k;
    }
  }
  return max_index;
}

double Signum(double x) {
  if (x > 0) {
    return 1;
  } else {
    if (x == 0) {
      return 0;
    }
    return -1;
  }
}
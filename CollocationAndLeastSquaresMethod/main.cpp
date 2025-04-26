#include"CollocationLeastSquaresMethod.h"

static const int N = 4;
static const int k = 20;// количество ячеек
static const int M = N + 5;

static const int N1 = N + 1;
static const int n = (N + 1) * k;
static const int m = (N + 5) * k;

static const double eps = pow(10, -12);

static const bool defKrylov="false";

const string result = {"result.txt"};
const string pseudoError = {"pseudoError.txt"};

int main() {


  MatrixXd C(n, 1);
  for (int i = 0; i < n; ++i) C(i, 0) = 0.4;//начальное значение

  MatrixXd C_prev;

  double error;

  size_t N_iter = 0;

  vector<double> cond_Ai(k);

  int num = 81;
  vector<MatrixXd> df_c;

  ofstream f{result};
  f.precision(2);

  ofstream f2{pseudoError};
  f2.precision(2);

  do {
    if (defKrylov && N_iter > 0) { // усорение Крылова
      if (df_c.size() < num) df_c.push_back(C);
      else {
        C = krylovAcceleration(df_c, num); // новое приближение
        df_c.clear();
      }
    }

    C_prev = C;
    for (size_t i = 0; i <= k - 1; ++i) {
      MatrixXd A(M, N1);
      MatrixXd b(M, 1);
      MatrixXd c_local(N1, 1);
      FillLocalMatrix(A, b, C, i);
      cond_Ai[i] = CondMatrix(A);
      QR_method(A, b, c_local);
      for (size_t j = 0; j <= N; ++j) {
        C(N * i + j + i, 0) = c_local(j, 0);
      }
    }

    vector<double> diff_C(n);
    for (size_t l = 0; l < n; ++l)diff_C[l] = abs(C(l, 0) - C_prev(l, 0));

    auto iterMaxElement = max_element(diff_C.begin(), diff_C.end());
    error = *iterMaxElement;
    ++N_iter;
    f2 << N_iter << " " << log2(error) << endl; //псевдопогрешность
  } while (!(error < eps));

  if (k == 10)
    Solution(C);

  pair<double, double> errorRelativeAbsolute = Error(C); //относительная-абсолютная погрешность

  f << "Absolute_error   "
    << "Relanive_error  "
    << "N_iter  " << endl;
  f << errorRelativeAbsolute.second << " " << errorRelativeAbsolute.first << " " << N_iter << endl;

  Time(f);

  for (size_t i = 0; i < k; i++) f << "cond(A_" << i << ")=" << cond_Ai[i] << endl;

  return 0;
}
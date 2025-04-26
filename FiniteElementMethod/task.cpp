// #include <iostream>
// #include <fstream>
// #include <Eigen/Dense>
// #include <cmath>
// #include <vector>
// #include <ctime>

// using namespace std;
// using namespace Eigen;

// static const string str={"/root/projects/task_3_dop/data_out.txt"};
// ofstream file{str};


// double u_exact(double x) {
//     return (2 * cos(1 - x) - sin(x)) / cos(1) + x * x - 2;
// }

// double f(double x) {
//     return x * x;
// }

// double u_exact_square(double x) {
//     return pow(u_exact(x), 2);
// }


// double ErL2(const Eigen::VectorXd& Uh, const Eigen::VectorXd& Utr) {
//     double num = 0.0, denom = 0.0;
//     for (int i = 0; i < Uh.size(); ++i) {
//         num += pow(Uh[i] - Utr[i], 2);
//         denom += pow(Utr[i], 2);
//     }
//     return std::sqrt(num / denom);
// }

// double Error_Inf(const VectorXd& Uh, const VectorXd& Utr) {
//     double max1 = 0.0, max2 = 0.0;
//     for (int i = 0; i < Uh.size(); ++i) {
//         max1 = max(max1, abs(Uh[i] - Utr[i]));
//         max2 = max(max2, abs(Utr[i]));
//     }
//     return max1 / max2;
// }

// void print_results(const VectorXd& Uh, const VectorXd& Utr, double time, double cond, double L2prev, double L_Inf_prev,ofstream& file) {
//     //double L2 = ErL2(Uh, Utr);
//     double L_Inf = Error_Inf(Uh, Utr);
//     //double RL2 = std::abs(std::log2(L2 / L2prev));
//     double R_Inf = abs(log2(L_Inf / L_Inf_prev));
//     file<<"Время: " << time << "\nError_Inf = " << L_Inf<< "\nСходимость Inf " << R_Inf << endl;
// }


// VectorXd find_solution(int N, double Le) {
//     MatrixXd K(2 * N + 1, 2 * N + 1);
//     MatrixXd P(2 * N + 1, 2 * N + 1);
//     K.setZero();
//     P.setZero();

//     MatrixXd A_P(3, 3);
//     A_P << 4, 2, -1,
//            2, 16, 2,
//            -1, 2, 4;
//     A_P *= Le / 30;

//     MatrixXd A_K(3, 3);
//     A_K << 7, -8, 1,
//            -8, 16, -8,
//            1, -8, 7;
//     A_K = A_K / (3 * Le) - A_P;

//     for (int i = 0; i < N; ++i) {
//         P.block<3, 3>(2 * i, 2 * i) -= A_P;
//         K.block<3, 3>(2 * i, 2 * i) += A_K;
//     }
    
//     K.row(0).setZero();
//     P.row(0).setZero();
//     K(0, 0) = 1.0;

//     K(2*N, 2*N-1) = -1.0 / Le; // Вклад от производной
//     K(2*N, 2*N) = 1.0 / Le;
//     P(2*N) = 1.0; // Значение производной

//     VectorXd f_vals(2*N + 1);
//     for (int i = 0; i < 2 * N + 1; ++i) {
//         f_vals[i] = f(i * Le);
//     }

//     // Решаем систему K * u = P * f
//     VectorXd u = K.colPivHouseholderQr().solve(P * f_vals);
    
//     return u;
// }

// int main() {
//     vector<int> grid = {2, 4, 8, 16, 32, 64, 128};
//     double L2_prev = 0, Inf_prev = 0;

//     for (int N : grid) {
//         file<<"N="<<N<<endl;
//         double start_time = clock();
//         double Le = 1.0 / N;
//         VectorXd u = find_solution(N, Le);

//         VectorXd X = VectorXd::LinSpaced(N * 100 + 1, 0, 1);
//         VectorXd Xi = VectorXd::LinSpaced(2 * N + 1, 0, 1);
//         VectorXd ξ = VectorXd::LinSpaced(101, -1, 1);

//         VectorXd Uh(N * 100 + 1);
//         Uh.setZero();

//         // Собираем решение в Uh
//         for (int i = 0; i < N; ++i) {
//             int start = 100 * i;
//             int end = 100 * (i + 1);
//             for (int j = 0; j < 100; ++j) {
//                 double xi = ξ[j];
//                 Uh[start + j] = (xi * (xi - 1) / 2) * u[2 * i] + ((1 - xi) * (1 + xi)) * u[2 * i + 1] + (xi * (xi + 1) / 2) * u[2 * i + 2];
//             }
//         }

//         // Время выполнения
//         double time = (clock() - start_time) / (double)CLOCKS_PER_SEC;

//         // Вычисляем точное решение для сравнения
//         VectorXd Utr(N * 100 + 1);
//         for (int i = 0; i < N * 100 + 1; ++i) {
//             Utr[i] = u_exact(i * Le);
//         }

//         // Вывод результатов
//         print_results(Uh, Utr, time, 0, L2_prev, Inf_prev,file);
//         L2_prev = ErL2(Uh, Utr);
//         Inf_prev = Error_Inf(Uh, Utr);
//     }

//     return 0;
// }






// #include <iostream>
// #include <Eigen/Dense>
// #include <vector>
// #include <cmath>
// #include <fstream>

// using namespace Eigen;
// using namespace std;

// static const string str={"/root/projects/task_3_dop/data_out.txt"};
// ofstream f{str};

// // Точное решение задачи
// double exactSolution(double x) {
//     return (2 * cos(1 - x) - sin(x)) / cos(1) + x * x - 2;
// }

// // Производная точного решения
// double exactDerivative(double x) {
//     return sin(1 - x) / cos(1) + 2 * x - cos(x) / cos(1);
// }

// // Функции формы для квадратичных элементов
// Vector3d shapeFunctions(double xi) {
//     return {(1 - xi) * (1 - 2 * xi), 4 * xi * (1 - xi), xi * (2 * xi - 1)};
// }

// // Производные функций формы
// Vector3d shapeFunctionDerivatives(double xi) {
//     return {4 * xi - 3, 4 - 8 * xi, 4 * xi - 1};
// }
// // Сборка локальной матрицы жесткости
// Matrix3d localStiffnessMatrix(double h) {
//     Matrix3d K_local= Matrix3d::Zero();
//     Matrix3d K_local_1 = Matrix3d::Zero();
//     Matrix3d K_local_2 = Matrix3d::Zero();
//     K_local_1<< 7, -8, 1,
//                 -8, 16, -8,
//                 1, -8,  7;
//     K_local_1*=1/(3*h);
//     K_local_2<< 4, 2, -1,
//                 2, 16, 2,
//                 -1, 2,  4;
//     K_local_2*=h/30;
//     K_local=K_local_1-K_local_2;
//     return K_local;
// }

// // Сборка локального вектора правой части
// Vector3d localLoadVector(double h, double x0) {
//     Vector3d F_local = Vector3d::Zero();
//     for (double xi : {-sqrt(1.0 / 3), sqrt(1.0 / 3)}) {
//         auto N = shapeFunctions(xi);
//         double x = x0 + h * (xi + 1) / 2.0;
//         F_local += N * (x * x) * h / 2.0;
//     }
//     return F_local;
// }


// double cond_matrix(MatrixXd* A){
// 	JacobiSVD<MatrixXd> svd(*A);
// 	double cond_A = svd.singularValues()(0) /svd.singularValues()(svd.singularValues().size()-1);
// 	return cond_A;
// }

// int main() {
//     int N= 8;
//     double L = 1.0;
//     double h = L / N;
//     int n = 2 * N + 1;
//     VectorXd u = VectorXd::Zero(n);
//     VectorXd rhs = VectorXd::Zero(n);
//     MatrixXd stiffness = MatrixXd::Zero(n, n);

//     // Сборка глобальной матрицы жесткости и вектора правой части
//     for (int j = 0; j < N; ++j) {
//         int i0 = 2 * j;
//         int i1 = 2 * j+ 1;
//         int i2 = 2 * j + 2;

//         Matrix3d K_local = localStiffnessMatrix(h);
//         Vector3d F_local = localLoadVector(h, j* h);

//         stiffness.block<3, 3>(i0, i0) += K_local;
//         rhs.segment<3>(i0) += F_local;
//     }
// // Учет краевых условий
//     stiffness.row(0).setZero();
//     stiffness(0, 0) = 1.0;
//     rhs(0) = 0.0;

//     stiffness.row(nNodes - 1).setZero();
//     stiffness(nNodes - 1, nNodes - 1) = 1.0;
//     rhs(nNodes - 1) = 1.0;

//     // Решение системы
//     u = stiffness.colPivHouseholderQr().solve(rhs);

//     // Вывод численного и точного решения
//     f<< "x\tnumeric\texact\terror\n";
//     for (int i = 0; i < nNodes; ++i) {
//         double x = i * h / 2.0;
//         double u_exact = exactSolution(x);
//         f<< x << "\t" << u[i] << "\t" << u_exact << "\t" << fabs(u[i] - u_exact) << "\n";
//     }
//     f<<cond_matrix(&stiffness)<<endl;

//     return 0;
// }
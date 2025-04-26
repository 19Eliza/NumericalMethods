#include"FiniteElementMethod.h"

using namespace std;
using namespace Eigen;

static const double L=1.0;

static const double Q=10000;

static const string approximateSolution={"approximateSolution.txt"};
static const string exactSolution={"exactSolution.txt"};

namespace FiniteElementMethod{
double F(double x){
	return pow(x,2);
}

double ExactSolution(double x){
    return ((2.0*cos(1-x)-sin(x))/cos(1.0))+x*x-2.0;
}


pair<MatrixXd,MatrixXd> SolveGlobalMatrix(const size_t N){
	double h=L/N;
    MatrixXd K(2*N + 1, 2*N + 1);
    MatrixXd P(2*N + 1, 2*N + 1);
    K.setZero();
    P.setZero();

    MatrixXd A_P(3, 3);
    A_P << 4, 2, -1,
           2, 16, 2,
           -1, 2, 4;
    A_P *= h / 30.0;

    MatrixXd A_K(3, 3);
    A_K << 7, -8, 1,
           -8, 16, -8,
           1, -8, 7;
	A_K*=1.0/(3.0 * h);

    A_K = A_K - A_P;

    for (int i = 0; i < N; ++i) {
        P.block<3, 3>(2 * i, 2 * i) -= A_P;
        K.block<3, 3>(2 * i, 2 * i) += A_K;
    }
    
    K.row(0).setZero();
    K(0, 0) = 1.0;

    VectorXd f_vals(2*N + 1);
	VectorXd Xi = VectorXd::LinSpaced(2*N + 1, 0, 1);
    for (int i = 0; i < 2*N + 1; ++i) {
        f_vals[i] = F(Xi[i]);
    }

	MatrixXd V(2*N+1,1);
	V.setZero();
	V(2*N)=1.0;

	P=P*f_vals+V;
	P.row(0).setZero();
	
    // K * u = P
    VectorXd u = K.colPivHouseholderQr().solve(P);
    
    return {u,K};
}


void Solution(MatrixXd& u,size_t N){
	ofstream f1{approximateSolution};
	ofstream f2{exactSolution};
    int q=20;
	VectorXd x = VectorXd::LinSpaced(q+1,0,1);
	double hh=L/q;

    double h=L/N;
	
	for(size_t i=0;i<=q;i++){
        double x_m=x[i];
		f1<<x_m<<" "<<ApproximateSolution(x_m,u,N,h)<<endl;
		f2<<x_m<<" "<<ExactSolution(x_m)<<endl;
	}

}


double Approximate_Solution(double x, const VectorXd& U, int N, double h) {
 
    int element = static_cast<int>(x / h);
	if (element >= N) {
            element = N - 1;
    }

    double x0 = element * h;
    double x1 = x0+h;
	double x_c=(x0+x1)/2.0;

    double N0_val = (x-x_c)*(x-x1)/((x0-x_c)*(x0-x1));
    double N1_val = (x-x0)*(x-x1)/((x_c-x0)*(x_c-x1));
    double N2_val = (x-x0)*(x-x_c)/((x1-x0)*(x1-x_c));

	int i=2*element;

    double u_approx = N0_val * U(i) + N1_val * U(i+ 1) + N2_val * U(i + 2);

    return u_approx;
}


double ErrorInfinite(MatrixXd& u,size_t N){
	double h=L/N;
	VectorXd x = VectorXd::LinSpaced(Q,0,1);
	double max_diff=0.0;
	double max_val_u=0.0;
	for(size_t i=0;i<=(Q-1);i++){
		double x_m=x[i];
		max_diff=max(max_diff,abs(exact_u(x_m)-Approximate_Solution(x_m,u,N,h)));
		max_val_u=max(max_val_u,abs(exact_u(x_m)));
	}
	double relative_error=max_diff/max_val_u;
	return relative_error;
}


double SimpsonMethod(const vector<double>& y,double hh) {
    int n = y.size() - 1;
    if (n % 2 != 0) {
        throw invalid_argument("Количество интервалов должно быть четным.");
    }

    double result = 0.0;

    for (int i = 0; i < (n-1); i += 2) {
        result += y[i] + 4.0 * y[i + 1] + y[i + 2];
    }

    return (hh / 3.0) * result;
}


double ErrorL2(MatrixXd& u,size_t N){
	double h=L/N;
	VectorXd x = VectorXd::LinSpaced(Q+1,0,1);
	vector<double> diff(x.size());
	vector<double> val_u(x.size());
	for(size_t i=0;i<=Q;i++){
		double x_m=x[i];
		diff[i]=pow(abs(exact_u(x_m)-Approximate_Solution(x_m,u,N,h)),2);
		val_u[i]=pow(abs(exact_u(x_m)),2);
	}
	double numerator;
	double denominator;
	double hh=1/Q;
	numerator=Simpson_Method(diff,hh);
	denominator=Simpson_Method(val_u,hh);

	double error_2=sqrt(numerator/denominator);

	return error_2;
}


double CondMatrix(MatrixXd& A){
	JacobiSVD<MatrixXd> svd(A);
	double cond_A = svd.singularValues()(0) /svd.singularValues()(svd.singularValues().size()-1);
	return cond_A;
}


void Time(ofstream& of){
	of<<"Time: "<<(double)clock()/CLOCKS_PER_SEC<<endl;
}
}



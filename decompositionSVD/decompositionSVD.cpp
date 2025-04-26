#include"decompositionSVD.h"

using namespace std;
using namespace Eigen;

static constexpr char delim{ '\n' };

enum method{QR_Method,Normal_Equations_Method,q};

const string data ="data_3.txt";

namespace SVD{

double P(double x,double j){
	return pow(x,j);
}

pair<double, VectorXd> powerIteration(MatrixXd A, int maxIterations, double tolerance) {
    int n = A.rows();
    VectorXd x = VectorXd::Random(n).normalized();
    double eigenvalue = 0.0;

    for (int k = 0; k < maxIterations; ++k) {
        VectorXd next_x = A * x;
        next_x.normalize();
        
        double next_eigenvalue = x.dot(A * x); 
        if (abs(next_eigenvalue - eigenvalue) < tolerance) {
            eigenvalue = next_eigenvalue;
            break;
        }
        x = next_x;
        eigenvalue = next_eigenvalue;
    }

    return {eigenvalue, x};
}

vector<pair<double, VectorXd>> findEigenvaluesAndVectors(MatrixXd A, int maxIterations , double tolerance ) {
    int n = A.rows();
    vector<pair<double, VectorXd>> results;

    for (int i = 0; i < n; ++i) {
        auto [lambda, eigenvector] = powerIteration(A, maxIterations, tolerance);
        results.emplace_back(lambda, eigenvector);
        A = A - lambda * (eigenvector * eigenvector.transpose());
    }

    return results;
}


bool FillMatrix(MatrixXd& A, MatrixXd& b){

	ifstream f1{data};

	size_t k=0;

	while (!f1.eof()) {
		string s;
		getline(f1, s, delim);

		string s1,s2;
		size_t i,j;
		
		for(i=0;i<s.size();i++){
			if(s[i] && s[i]!=' '){
				s1+=s[i];
			}
			else{
				break;
			}
		}

		for(j=i+1;j<s.size();j++){
			if(s[j] && s[j]!=' '){
				s2+=s[j];
			}
			else{
				break;
			}
		}

		float eps=stof(s1);

		for(size_t j=0;j<A->cols();j++){
			A(k,j)=P(eps,j);
		}
		
		float sigm=stof(s2);
		b(k,0)=sigm;
		
		++k;
	}

	f1.close();

	return true;
}

}



#include"decompositionSVD.h"

static const int M{7806}; //data_3-7806//data_4-1000//data_5-1000
static const int n{12};//количство столбцов

using namespace std;
using namespace Eigen;
using namespace SVD;

static const string result={"result.txt"};

int main(){
    ofstream f{result};

    double time=0.0;

    for(int N=2;N<=n;N++){
        f<<"N="<<N<<endl;
        MatrixXd A(M,N);
        MatrixXd b(M,1);

        FillMatrix(A,b);

        MatrixXd AtA = A.transpose() * A;

        auto results = findEigenvaluesAndVectors(AtA,1000,1e-10);

        VectorXd eigenvalues(results.size());

        for(size_t i=0;i<eigenvalues.size();i++) eigenvalues(i)=results[i].first;
        

        size_t size=results[0].second.size();

        MatrixXd eigenvectors(size,results.size());

        for(size_t i=0;i<eigenvalues.cols();i++)  eigenvectors.col(i)=results[i].second;
        

        vector<int> indices(eigenvalues.size());
        for (int i = 0; i < indices.size(); ++i) indices[i] = i;
        

        sort(indices.begin(), indices.end(), [&eigenvalues](int a, int b) {
            return eigenvalues[a] > eigenvalues[b];
        });

        VectorXd sorted_eigenvalues(eigenvalues.size());
        MatrixXd sorted_eigenvectors(eigenvectors.rows(), eigenvectors.cols());
        for (int i = 0; i < indices.size(); ++i) {
            sorted_eigenvalues(i) = eigenvalues(indices[i]);
            sorted_eigenvectors.col(i) = eigenvectors.col(indices[i]);
        }

        MatrixXd V = sorted_eigenvectors; // Матрица V

        
        VectorXd singularValues = sorted_eigenvalues.cwiseSqrt();
        cout << "Сингулярные значения (Σ):\n" << singularValues << "\n\n";
    
        MatrixXd Sigma_inv = MatrixXd::Zero(A.cols(), A.rows());
        for (int i = 0; i < singularValues.size(); ++i) {
            if (singularValues(i)!=0.0) {
                Sigma_inv(i, i) = 1.0/singularValues(i);
            }
        }


        MatrixXd U = A * V * Sigma_inv;
        double maxSingularValue = singularValues(0);
        double minSingularValue = singularValues(singularValues.size() - 1);

        double cond_A=maxSingularValue/minSingularValue;
        double cond_AtA=cond_A*cond_A;

        MatrixXd c=V*Sigma_inv*U.transpose()*b;

        VectorXd diff=A*c-b;

        
        double sum=0.0;
        for(size_t i=0;i<M;i++){
            sum+=abs(diff(i))*abs(diff(i));
        }
        sum/=M;
        double NRMSE;
        NRMSE=sqrt(sum)/maxSingularValue;

        f<<scientific<<"Nrmse: "<<NRMSE<<"  "<<cond_A<<" "<<cond_AtA<<endl;
        time=Time(f,time);
        f<<endl;

    }

    return 0;
}















// int main() {
//     double time=0.0;
//     f.precision(2);
//     f<<"Nrmse:  "<<"  Cond_A:   "<<"  Cond_AtA: "<<endl;

//     for(int N=2;N<=n;N++)
//     {
//     f<<"N="<<N<<endl;
//     MatrixXd A(M,N);
// 	MatrixXd b(M,1);

// 	fill_matrix(&A,&b);

//     MatrixXd AtA = A.transpose() * A;

//     SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(AtA);
//       if (eigenSolver.info() != Eigen::Success) {
//         cerr << "Ошибка при вычислении собственных значений и векторов!" <<endl;
//         return -1;
//     }

//     VectorXd eigenvalues = eigenSolver.eigenvalues();
//     MatrixXd eigenvectors = eigenSolver.eigenvectors();

//     vector<int> indices(eigenvalues.size());
//     for (int i = 0; i < indices.size(); ++i) {
//         indices[i] = i;
//     }

//     sort(indices.begin(), indices.end(), [&eigenvalues](int a, int b) {
//         return eigenvalues[a] > eigenvalues[b];
//     });

//     VectorXd sorted_eigenvalues(eigenvalues.size());
//     MatrixXd sorted_eigenvectors(eigenvectors.rows(), eigenvectors.cols());
//     for (int i = 0; i < indices.size(); ++i) {
//         sorted_eigenvalues(i) = eigenvalues(indices[i]);
//         sorted_eigenvectors.col(i) = eigenvectors.col(indices[i]);
//     }

//     MatrixXd V = sorted_eigenvectors; // Матрица V

//     VectorXd singularValues = sorted_eigenvalues.cwiseSqrt();
//     std::cout << "Сингулярные значения (Σ):\n" << singularValues << "\n\n";
    
//     MatrixXd Sigma_inv = MatrixXd::Zero(A.cols(), A.rows());
//     for (int i = 0; i < singularValues.size(); ++i) {
//         if (singularValues(i)!=0.0) {
//             Sigma_inv(i, i) = 1.0/singularValues(i);
//         }
//     }


//     MatrixXd U = A * V * Sigma_inv;

//     double maxSingularValue = singularValues(0);
//     double minSingularValue = singularValues(singularValues.size() - 1);

//     double cond_A=maxSingularValue/minSingularValue;
//     double cond_AtA=cond_A*cond_A;

//     MatrixXd c=V*Sigma_inv*U.transpose()*b;

//     VectorXd diff=A*c-b;

    
//     double sum=0.0;
//     for(size_t i=0;i<M;i++){
//         sum+=abs(diff(i))*abs(diff(i));
//     }
//     sum/=M;
//     double NRMSE_;
//     NRMSE_=sqrt(sum)/maxSingularValue;

    
//     f<<scientific<<"Nrmse: "<<NRMSE_<<" "<<cond_A<<" "<<cond_AtA<<endl;
//     time=Time(f,time);
//     f<<endl;

//     if(N==3) Solution(A,c,b);

//     }

//     f.close();
//     return 0;
// }
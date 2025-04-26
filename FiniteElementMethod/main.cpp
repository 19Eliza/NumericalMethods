#include"FiniteElementMethod.h"

using namespace Eigen;
using namespace FiniteElementMethod;
using namespace std;

static const size_t NN=128;
static const double L=1.0;

static const string fileName={"result.txt"};

int main(){

    ofstream file{fileName};
    file.precision(2);

    double R_L2Prev=0.0;
    double R_InfPrev=0.0;

    double errorInfPrev=0.0;
    double errorL2Prev=0.0;

    file<<"Error_Inf:\tR_Inf\tError_L2\tR_L2\tCond(K)\n";

    for(size_t N=2;N<=NN;N*=2){
    
        file<<"Element number N="<<N<<endl;

        const size_t n=2*N+1;
        const double h=L/N;

        MatrixXd u(n,1);
    
        auto solution=SolveGlobalMatrix(N);

        u=solution.first;

        auto K=solution.second;
       
        double errorInf=ErrorInfinite(u,N);
        double errorL2=ErrorL2(u,N);

        file<<scientific<<error_inf<<"\t"<<log2(errorInfPrev/errorInf)<<"\t"<<
                        error_L2<<"\t"<<log2(errorL2Prev/errorL2)<<"\t"<<CondMatrix(K)<<endl;
                        
        errorInfPrev=errorInf;
        errorL2Prev=errorL2;

        if(N==2)Solution(u,N);

        Time(file);  
        
        file<<endl;

    }

    file.close();
    return 0;
}

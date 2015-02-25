#include"Sparse.h"

int main(){

    Sparse A = Sparse(3,3);
    //A(0,0) = 2;
    A.setVal(0,0,5);
    A.setVal(0,1,10); 
    A.setVal(1,1,4); 

    A.setVal(2,2,3); 
    std::vector<std::vector<double> > Adense;
    Adense = A.dense();

    for(int i=0; i<3; i++)
    {
        std::printf("[");
        for(int j=0; j<3; j++)
        {
            std::printf("%f", A(i,j));
            //std::printf("%f", Adense[i][j]);
            if( j != 2)
                std::printf(", ");
            else
                std::printf("]");
        }
        std::printf(endl);
    }
    
    std::vector<double> xx (3,1);
    
    for(int i=0; i<3; i++)
    {
        std::cout << xx[i] << endl;
    }
    std::vector<double> b = A*xx;
    for(int i=0; i<3; i++)
    {
        std::cout << b[i] << endl;
    }

    return 0;
}


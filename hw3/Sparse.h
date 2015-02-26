#include<iostream>
#include<cstdio>
#include<vector>
#include<string>
#include<map>

#define endl "\n"

class Sparse {

    private:
        
        // spase matrix indexes and elements
        std::vector<std::map<int,double> > _vals;

        // matrix dimensions
        int _M;
        int _N;

    public:

        Sparse(int M, int N);
        virtual ~Sparse();
        std::vector<double> operator * (const std::vector<double> x) const;
        void setVal(int i, int j, double value);
        double getVal(int i, int j);
        std::vector<std::vector<double> > dense(void);
        double operator() (int i, int j);
        void display(void);
        Sparse transpose(void);
        int size(void);
};

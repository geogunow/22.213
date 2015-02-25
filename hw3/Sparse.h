#include"utils.h"

/*
   Tuple structure for storing index, value pairs
   */
template <typename T>
struct Tuple
{
    int index;
    T value;
};

#define endl "\n"

class Sparse {

    private:
        
        // spase matrix indexes and elements
        std::vector<std::vector<Tuple<double> > > _vals;

        // matrix dimensions
        int _M;
        int _N;

    public:

        Sparse(int M, int N);
        virtual ~Sparse();
        std::vector<double> operator * (const std::vector<double> x) const;
        double & setVal(int i, int j, double value);
        double getVal(int i, int j);
        std::vector<std::vector<double> > dense(void);
        double operator() (int i, int j) const;
        double & operator() (int i, int j);
};

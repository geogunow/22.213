#ifndef __SOLUTIONS__
#define __SOLUTIONS__
#include<vector>
#include"utils.h"

class eigenSolution{
    public:
        std::vector<double> flux;
        std::vector< std::vector<double> > gFlux;
        std::vector<double> power;
        std::vector< std::vector<double> > gPower;
        double keff;
        int outer_iters;
        int inner_iters;
        
        eigenSolution();
        virtual ~eigenSolution();
        void indexArrays(Indexer index);
        std::vector<double> getFlux(int g);
        std::vector<double> getPower(int g);
};

class rkSolution{
    public:
        std::vector<std::vector<double> > flux;
        std::vector<std::vector<std::vector<double> > > gFlux;
        std::vector<std::vector<double> > powerProfile;
        std::vector<std::vector<std::vector<double> > > gPowerProfile;
        std::vector<double> power;
        int unused;

        rkSolution();
        virtual ~rkSolution();
        void indexArrays(Indexer index, int T);
        std::vector<double> getFlux(int t, int g);
        std::vector<double> getPowerProfile(int t, int g);
        std::vector<double> getRawPowerProfile(int t);
        std::vector<double> getPower();
};


#endif

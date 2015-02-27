#ifndef __XSDATA__
#define __XSDATA__

#include<iostream>
#include<vector>

#define endl "\n"

class XSdata{
    public:
        std::vector<double> D;
        std::vector<double> siga;
        std::vector<std::vector<double> > sigs;
        std::vector<double> nuSigf;
        std::vector<double> chi;
        int G;

        XSdata();
        virtual ~XSdata();
        void setD(double * xs, int ng);
        void setSiga(double * xs, int ng);
        void setSigs(double ** xs, int ng1, int ng2);
        void setNuSigf(double * xs, int ng);
        void setChi(double * xs, int ng);
        void printXS();
};

#endif

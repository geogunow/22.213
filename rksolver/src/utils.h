#ifndef __UTILS__
#define __UTILS__

#include<iostream>
#include<cstdio>
#include<vector>
#include<string>

#define endl "\n"

class Indexer{
    public:
        int _N;
        int _G;
        Indexer(int N, int G)
        {
            _N = N;
            _G = G;
        };
        virtual ~Indexer() {};
        int operator() (int n, int g){
            return _N*g + n;
        };
};


double dot(std::vector<double> x1, std::vector<double> x2);
void printVector(std::vector<double> x);

#endif

#ifndef MATRIX_H_
#define MATRIX_H_


#include <vector> 
#include "head.h"
using namespace std;

struct Matrix
{
    /* Matrix */
    vector<int> ia;
    vector<double> di;
    vector<double> al;
    vector<double> au;
};


vector<double> MultAOnq(Matrix &matr, vector<double>& q);
double calcNormE(const vector1D &x);


#endif
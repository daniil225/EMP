#ifndef SLAU_H_
#define SLAU_H_

#include "Matrix.h"
#include <vector>

using namespace std;



struct Slau
{
    Matrix Matr;
    /* f */
    vector<double> f;
};

void LUDecomposition(Slau &slau);
void GausForward(Slau &slau,vector<double> &y);
void GausBack(Slau &slau, vector<double> &x);
void SolveSlau(Slau &slau ,vector<double>&x);
void TestSlau();


#endif
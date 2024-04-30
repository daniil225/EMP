#ifndef SLAU_H_
#define SLAU_H_

#include "Matrix.h"
#include <vector>


/* SLAU с матрицей хранящейся в Разряженном блочном формате */
struct SLAU_BlockMatrix
{
    int32_t N = -1;
    int32_t size = -1;
    SparseBlockOfMatrix matrix;
    std::vector<std::vector<double>> b;
    std::vector<double> x;
};


/* SLAU с матрицей хранящейся в Разряженном формате */
struct SLAU_SparseMatrix
{
    int32_t N = -1;
    int32_t size = -1;
    SparseMatrix matrix;
    std::vector<double> b;
    std::vector<double> x;
};

/* SLAU с матрицей хранящейся в Профильном формате */
struct SLAU_ProfileMatrix
{
    int32_t N = -1;
    int32_t size = -1;
    ProfileMatrix Matr;
    std::vector<double> f;
    std::vector<double> x;
};


void LUDecomposition(SLAU_ProfileMatrix &slau);
void GausForward(SLAU_ProfileMatrix &slau,std::vector<double> &y);
void GausBack(SLAU_ProfileMatrix &slau, std::vector<double> &x);
void SolveSlau(SLAU_ProfileMatrix &slau ,std::vector<double>&x);

#endif
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
    ProfileMatrix matrix;
    std::vector<double> b;
    std::vector<double> x;
};

#endif
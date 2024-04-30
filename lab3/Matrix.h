#ifndef MATRIX_H_
#define MATRIX_H_

#include <cctype>
#include <vector>
#include <iostream>
/* Блок в матрице */
struct Block
{
    double pij = 0;
    double cij = 0;

    Block() = default;
    Block(double pij, double cij)
    {
        this->pij = pij;
        this->cij = cij;
    }

    ~Block() = default;
};



/* Блочная матрица в профильном формате  */
struct SparseBlockOfMatrix
{
    int32_t N = -1;          // Размер матрицы слау
    int32_t size = -1;       // Размер матриц  ggl, ggu, jg
    std::vector<Block> ggu;  // size
    std::vector<Block> ggl;  // size
    std::vector<Block> di;   // N
    std::vector<int32_t> ig; // N+1
    std::vector<int32_t> jg; // size

    SparseBlockOfMatrix() = default;
    SparseBlockOfMatrix(int32_t N, int32_t size);

    void AllocMemory(int32_t N, int32_t size);

    void PrintDenseMatrix();

    void InsertBlock(Block &elem, int32_t i, int32_t j);
    
    void AddBlock(Block &elem, int32_t i, int32_t j);
    bool IsBlockInMatrix(int32_t i, int32_t j);

    ~SparseBlockOfMatrix() = default;
};


struct SparseMatrix
{

    int32_t N = -1;
    int32_t size = -1;
    std::vector<double> di;
    std::vector<double> ggu;
    std::vector<double> ggl;
    std::vector<int32_t> ig;
    std::vector<int32_t> jg;

    SparseMatrix() = default;

    SparseMatrix(int32_t N, int32_t size);

    void AllocMemory(int32_t N, int32_t size);

    void PrintDenseMatrix();

    void InsertElem(double elem, int32_t i, int32_t j);

};

struct ProfileMatrix
{
    int32_t N = -1;
    int32_t size = -1;
    
    std::vector<int> ia;
    std::vector<double> di;
    std::vector<double> al;
    std::vector<double> au;

    ProfileMatrix() = default;

    ProfileMatrix(int32_t N, int32_t size);
    void AllocateMemory(int32_t N, int32_t size);

    void InsertElem(double elem, int32_t i, int32_t j);

    ~ProfileMatrix() = default;
};


Block GetBlockMatrix(SparseBlockOfMatrix &Matr, int32_t i, int32_t j);
double GetElementSparsematrix(SparseMatrix &Matr, int32_t i, int32_t j);
void BlockMatrix2SparseMatrix(SparseBlockOfMatrix &src_Matr, SparseMatrix& dst_Matr);

void SparseMatrix2ProfileMatrix(SparseMatrix &src_Matr, ProfileMatrix& dst_Matr);

#endif
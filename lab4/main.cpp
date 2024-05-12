#include <iostream>
#include "Matrix.h"
#include "SlauSolver.h"
#include "SLAU.h"

void SaveVec(const std::string &filename, int32_t precession)
{

}


int main()
{   


    SLAUSolvers::IterationSolvers::InitDateFile InitStruct;

    InitStruct.type = SLAUSolvers::MatrixType::SPARSE;
    InitStruct.di = "Test1/di.txt";
    InitStruct.ig = "Test1/ig.txt";
    InitStruct.jg = "Test1/jg.txt";
    InitStruct.ggl = "Test1/ggl.txt";
    InitStruct.ggu = "Test1/ggu.txt";
    InitStruct.kuslau = "Test1/kuslau.txt";
    InitStruct.f = "Test1/f.txt";
    
    SLAUSolvers::IterationSolvers::SLAU slau;

    SLAUSolvers::IterationSolvers::Load(InitStruct, slau, true);

    double x[4394];
    for(int i = 0; i < 4394; i++)
    {
        x[i] = i;
    }

    for(int i = 0; i < 4394; i++)
    {
        if(slau.matrix.di[i] > 0)
            slau.matrix.di[i] = 0.1;
    }

    SLAUSolvers::IterationSolvers::MultA(slau.matrix, x, slau.f);

    SLAUSolvers::IterationSolvers::SetSolveMode(slau, SLAUSolvers::IterationSolvers::Solvemode::LOS_NOSYMETRIC_LU_FACT);
    SLAUSolvers::IterationSolvers::SolveSLAU(slau, true);

    SLAUSolvers::IterationSolvers::SaveX(slau, "Test1/x.txt", "\n");

    /* Преобразуем внутренний формат матрицы решателя SLAU в формат для конвертации */
    // Matrix::SparseMatrix SparseMatrix(slau.N, slau.size);
    // int N = slau.N;
    // int size = slau.size;
    // std::copy(slau.matrix.di, slau.matrix.di + N, SparseMatrix.di.begin());
    // std::copy(slau.matrix.ggl, slau.matrix.ggl + size, SparseMatrix.ggl.begin());
    // std::copy(slau.matrix.ggu,slau.matrix.ggu + size, SparseMatrix.ggu.begin() );
    // std::copy(slau.matrix.ig,slau.matrix.ig + N+1, SparseMatrix.ig.begin() );
    // std::copy(slau.matrix.jg,slau.matrix.jg + size, SparseMatrix.jg.begin() );
    //SparseMatrix.PrintDenseMatrix("Sparse.txt");

    // SparseMatrix.PrintDenseMatrix("Sparse.txt");
    // Конвертим в Профильный формат 
    // SLAU_ProfileMatrix LUSolver;
    // //Matrix::SparseMatrix2ProfileMatrix(SparseMatrix, LUSolver.Matr);
    // LUSolver.x.resize(N);
    // LUSolver.f.resize(N);
    // LUSolver.N = N;
    // LUSolver.size = LUSolver.Matr.size;
    // std::copy(slau.f, slau.f+N, LUSolver.f.begin());

    // /* LU решатель */
    // SolveSlau(LUSolver, LUSolver.x);

    std::cout << "\n";

    return 0;
}
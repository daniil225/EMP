#include <iostream>

#include "Matrix.h"
#include "Grid3D_StreightQuadPrismatic.h"
#include "SLAU.h"
#include "GenerateMatrixPortrait.h"
#include "FEMSolver.h"
#include "Test.h"

ParamDE Test_1 = Test1();

int main()
{

    FEMSolver fem_solver("Area1.txt", Test_1);
    fem_solver.Calc();
    fem_solver.DivideGrid(2);
    
    //SparseBlockOfMatrix Matrix;

    // Matrix.AllocMemory(6,6);
    // Matrix.ig = {1,1,2,2,4,5,7};
    // Matrix.jg = {1,1,2,4,2,4};

    // for(auto &e: Matrix.ig)
    //     e--;
    
    // for(auto &e: Matrix.jg)
    //     e--;
    
    // Block elem;
    
    // for(int i = 0; i < Matrix.N; i++)
    // {
    //     for(int j = 0; j < Matrix.N; j++)
    //     {
    //         elem.pij = i+j;
    //         elem.cij = 2*(i+j);
    //         Matrix.InsertBlock(elem, i, j);
    //     }
    // }
    
    // Matrix.PrintDenseMatrix();
    // SparseMatrix matr;
    // BlockMatrix2SparseMatrix(Matrix, matr);

    // matr.PrintDenseMatrix();
   
    return 0;
}
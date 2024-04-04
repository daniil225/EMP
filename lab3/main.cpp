#include <iostream>

#include "Matrix.h"

int main()
{

    // Grid3D_StreightQuadPrismatic Grid;

    // GridStatus Status = Grid.Load("Area1.txt");
    // if(Status.GetState() == State::OK)
    // {
    //     Status = Grid.GenerateGrid();
    //     if(Status.GetState() == State::OK)
    //     {
    //         //Grid3D_Size GridSize = Grid.GetGridSize();
    //         Grid.PrintGridSlice(0);
    //         //InfoManeger::PrintInfo(Grid[3].info);
    //         Finit_Element_StreightQuadPrismatic Element = Grid.GetElement(0.0666666666668,0.5, 0.25);
    //         Element.PrintElement();
            
    //     }
    //     else
    //     {
    //         cout << Status.GetMsg();
    //     }
    // }
    // else
    // {
    //     cout << Status.GetMsg();
    // }

    SparseBlockOfMatrix Matrix;

    Matrix.AllocMemory(6,6);
    Matrix.ig = {1,1,2,2,4,5,7};
    Matrix.jg = {1,1,2,4,2,4};

    for(auto &e: Matrix.ig)
        e--;
    
    for(auto &e: Matrix.jg)
        e--;
    
    Block elem;
    
    for(int i = 0; i < Matrix.N; i++)
    {
        for(int j = 0; j < Matrix.N; j++)
        {
            elem.pij = i+j;
            elem.cij = 2*(i+j);
            Matrix.InsertBlock(elem, i, j);
        }
    }
    
    Matrix.PrintDenseMatrix();
    SparseMatrix matr;
    BlockMatrix2SparseMatrix(Matrix, matr);

    matr.PrintDenseMatrix();
   
    return 0;
}
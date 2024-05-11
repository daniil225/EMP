#include <iostream>
#include "Matrix.h"
#include "SlauSolver.h"

int main()
{   


    SLAUSolvers::IterationSolvers::InitDateFile InitStruct;

    InitStruct.type = SLAUSolvers::MatrixType::SPARSE;
    InitStruct.di = "Test/di.txt";
    InitStruct.ig = "Test/ig.txt";
    InitStruct.jg = "Test/jg.txt";
    InitStruct.ggl = "Test/ggl.txt";
    InitStruct.ggu = "Test/ggu.txt";
    InitStruct.kuslau = "Test/kuslau.txt";
    InitStruct.f = "Test/f.txt";
    
    SLAUSolvers::IterationSolvers::SLAU slau;

    SLAUSolvers::IterationSolvers::Load(InitStruct, slau, true);

    double x[6] = {1,2,3,4,5,6};
    SLAUSolvers::IterationSolvers::MultA(slau.matrix, x, slau.f);

    SLAUSolvers::IterationSolvers::SetSolveMode(slau, SLAUSolvers::IterationSolvers::Solvemode::LOS_NOSYMETRIC_LU_FACT);
    SLAUSolvers::IterationSolvers::SolveSLAU(slau, true);

    SLAUSolvers::IterationSolvers::SaveX(slau, "Test/x.txt", "\n");


    

    return 0;
}
#include <iostream>

#include "FEMSolver.h"
#include "Test.h"
#include "Integrate_3DStreightQuadPrismatic.h"
#include "SlauSolver.h"
#include <chrono>
#include <iomanip>
#include <unistd.h>

typedef std::function<double(double, double, double)> Func_3D;

BaseGrid3D_Integrate_StreightQuadPrismatic BaseGrid2BaseGridIntegrate(const BaseGrid3DStreightQuadPrismatic &baseGrid)
{
    /* Тупо копируем все содержимое в структуру и все */
    BaseGrid3D_Integrate_StreightQuadPrismatic res;

    res.Nx = baseGrid.Nx;
    res.Ny = baseGrid.Ny;
    res.Nz = baseGrid.Nz;

    res.BaseGridXZ.resize(res.Nz);
    for (int32_t i = 0; i < res.Nz; i++)
        res.BaseGridXZ[i].resize(res.Nx);

    for (int32_t i = 0; i < baseGrid.Nz; i++)
    {
        for (int32_t j = 0; j < baseGrid.Nx; j++)
        {
            res.BaseGridXZ[i][j].x = baseGrid.BaseGridXZ[i][j].x;
            res.BaseGridXZ[i][j].z = baseGrid.BaseGridXZ[i][j].z;
        }
    }

    res.BaseGridY.resize(res.Ny);
    res.BaseGridY = baseGrid.BaseGridY;

    res.DivideParam.resize(3);
    res.DivideParam[0].resize(res.Nx - 1);
    res.DivideParam[1].resize(res.Nz - 1);
    res.DivideParam[2].resize(res.Ny - 1);

    for (int32_t i = 0; i < baseGrid.DivideParam[0].size(); i++)
    {
        res.DivideParam[0][i].coef = baseGrid.DivideParam[0][i].coef;
        res.DivideParam[0][i].num = baseGrid.DivideParam[0][i].num;
    }

    for (int32_t i = 0; i < baseGrid.DivideParam[1].size(); i++)
    {
        res.DivideParam[1][i].coef = baseGrid.DivideParam[1][i].coef;
        res.DivideParam[1][i].num = baseGrid.DivideParam[1][i].num;
    }

    for (int32_t i = 0; i < baseGrid.DivideParam[2].size(); i++)
    {
        res.DivideParam[2][i].coef = baseGrid.DivideParam[2][i].coef;
        res.DivideParam[2][i].num = baseGrid.DivideParam[2][i].num;
    }

    res.isReadyToUse = true;

    return res;
}

ParamDE Test_1 = Test1();
ParamDE Test_2 = Test2();


int main()
{

    std::vector<double> omega = {10e5, 10e12};

    int idx = 2;

    Test_2.omega = omega[idx - 1];

    FEMSolver fem_solver("Area1.txt", Test_2);
    fem_solver.DivideGrid(6);
    fem_solver.ClearAll(); /* Очистили внутреннее состояние  */
    fem_solver.PreCalc();  /* Подготовили слау */


    SLAU_SparseMatrix slau_sparse = fem_solver.GetSparseSLAU();
    SLAU_ProfileMatrix slau_profile;
    //SparseMatrix2ProfileMatrix(slau_sparse.matrix, slau_profile.Matr);


    SLAUSolvers::IterationSolvers::SLAU iterslau;
    SLAUSolvers::IterationSolvers::Allocate_Memory_SLAU(iterslau, slau_sparse.N, slau_sparse.size, SLAUSolvers::MatrixType::SPARSE);
    /* Теперь инициализация самой слау ну тут тупо кипирование */
    std::copy(slau_sparse.matrix.di.begin(), slau_sparse.matrix.di.end(), iterslau.matrix.di);
    std::copy(slau_sparse.matrix.ggl.begin(), slau_sparse.matrix.ggl.end(), iterslau.matrix.ggl);
    std::copy(slau_sparse.matrix.ggu.begin(), slau_sparse.matrix.ggu.end(), iterslau.matrix.ggu);
    std::copy(slau_sparse.matrix.ig.begin(), slau_sparse.matrix.ig.end(), iterslau.matrix.ig);
    std::copy(slau_sparse.matrix.jg.begin(), slau_sparse.matrix.jg.end(), iterslau.matrix.jg);
    std::copy(slau_sparse.b.begin(), slau_sparse.b.end(), iterslau.f);
    iterslau.maxiter = 15000;
    iterslau.eps = 1e-10;

    std::vector<double> x(slau_sparse.N);
    double elem = -1.0;

    for(int32_t i = 0; i < slau_sparse.N; i++)
    {
        if(i % 25 == 0)
            elem += 1.0;

        x[i] = elem;
    }


    SLAUSolvers::IterationSolvers::MultA(iterslau.matrix, x.data(), iterslau.f);


    cout << "Save matrix\n";
    SLAUSolvers::IterationSolvers::InitDateFile Initfile;

    Initfile.di = "saveTest/di.txt";
    Initfile.f = "saveTest/f.txt";
    Initfile.kuslau = "saveTest/kuslau.txt";
    Initfile.ggl = "saveTest/ggl.txt";
    Initfile.ggu = "saveTest/ggu.txt";
    Initfile.ig = "saveTest/ig.txt";
    Initfile.jg = "saveTest/jg.txt";



    /* Сохраняем матрицу СЛАУ и правую часть ее мы предварительно сгенерируем */

    std::ofstream out;

    out.open(Initfile.di);
    for(auto e: slau_sparse.matrix.di)
        out << std::setprecision(16) <<  e << "\n";
    out.close();


    out.open(Initfile.f);
    for(int32_t i = 0; i < slau_sparse.N; i++)
        out << std::setprecision(16) << iterslau.f[i] << "\n";
    out.close();


    out.open(Initfile.ggl);
    for(auto e: slau_sparse.matrix.ggl)
        out << std::setprecision(16) << e << "\n";
    out.close();


    out.open(Initfile.ggu);
    for(auto e: slau_sparse.matrix.ggu)
        out << std::setprecision(16) << e << "\n";
    out.close();


    out.open(Initfile.kuslau);
    out <<  fem_solver.GetNodeCount()*2 << "\n";
    out << 15000 << "\n" << 1e-10 << "\n";
    out.close();


    out.open(Initfile.ig);
    for(auto e: slau_sparse.matrix.ig)
        out << e << "\n";
    out.close();



    out.open(Initfile.jg);
    for(auto e: slau_sparse.matrix.jg)
        out << e << "\n";
    out.close();

    out.open("saveTest/x.txt");
    for(auto e: x)
        out << std::setprecision(16) << e << "\n";
    out.close();

    return 0;
}
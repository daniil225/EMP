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

std::vector<double> Solve(ParamDE Test, int32_t dividecoef = 1.0, int32_t maxiter = 100, double eps = 1e-10, bool printInterSlau = true, bool printNodeError = true)
{
    std::vector<double> res(4);

    cout << "\n";
    cout << "------------------------------------------------------------------\n";
    cout << "Solve start\n";
    auto begin = std::chrono::steady_clock::now();

    FEMSolver fem_solver("Area1.txt", Test);
    BaseGrid3D_Integrate_StreightQuadPrismatic IntegGrid = BaseGrid2BaseGridIntegrate(fem_solver.GetBaseGrid());

    fem_solver.DivideGrid(dividecoef);
    fem_solver.ClearAll(); /* Очистили внутреннее состояние  */
    fem_solver.PreCalc();  /* Подготовили слау */
    fem_solver.Calc(maxiter, eps, printInterSlau, printNodeError);

    /* Функции для проверки сходимости */
    Func_3D f_norm_u_c = [&](double x, double y, double z) -> double
    { return std::pow(fem_solver.u_c(x, y, z) - Test.u_c(x, y, z), 2.0); };
    Func_3D f_norm_u_s = [&](double x, double y, double z) -> double
    { return std::pow(fem_solver.u_s(x, y, z) - Test.u_s(x, y, z), 2.0); };

    Integrate_3DStreightQuadPrismatic Integ_c(IntegGrid, f_norm_u_c);
    Integrate_3DStreightQuadPrismatic Integ_s(IntegGrid, f_norm_u_s);

    std::pair<double, double> fn_c = Integ_c.Quad("default", 1e-7, 4);
    std::pair<double, double> fn_s = Integ_s.Quad("default", 1e-7, 4);

    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);

    cout << "Solve end. Total time: " << fixed << std::setprecision(2) << elapsed_ms.count() / 1000.0 << "s.\n";
    cout << "Finit Element Count = " << fem_solver.GetNodeCount() << "\n";
    cout << "f_norm_u_c = " << std::scientific << std::setprecision(4) << std::sqrt(fn_c.first) << " ErrorInteg = " << fn_c.second << "\n";
    cout << "f_norm_u_s = " << std::sqrt(fn_s.first) << " ErrorInteg = " << fn_s.second << "\n";
    cout << "General Error = " << std::sqrt(fn_s.first + fn_c.first) << "\n";
    cout << "------------------------------------------------------------------\n";

    res[0] = fem_solver.GetNodeCount();
    res[1] = std::sqrt(fn_c.first);
    res[2] = std::sqrt(fn_s.first);
    res[3] = std::sqrt(fn_c.first + fn_s.first);

    return res;
}

int main()
{

    
    // std::vector<int32_t> dividecoefs = {1, 2, 4, 6};

    // std::vector<std::vector<double>> res_test;
    // for(int32_t dc: dividecoefs)
    // {
    //     res_test.push_back(Solve(Test_2, dc, 3000, 1e-10, true, false));
    // }

    // std::ofstream out("Table.md");
    // // Генерация табличных данных
    // out << "## Table\n";
    // out << "|N|Node|f_norm_u_c|f_norm_u_s|f_norm_general|сходимость|\n";
    // out << "|:-:|:-:|:-:|:-:|:-:|:-:|\n";
    // out << "|" << 1 << "|" << int32_t(res_test[0][0]) << "|" << res_test[0][1] << "|" << res_test[0][2] << "|" << res_test[0][3] << "|" << "-" << "|\n";
    // out << "|" << 2 << "|" << int32_t(res_test[1][0]) << "|" << res_test[1][1] << "|" << res_test[1][2] << "|" << res_test[1][3] << "|" << std::log2(res_test[0][3]/res_test[1][3]) << "|\n";
    // out << "|" << 3 << "|" << int32_t(res_test[2][0]) << "|" << res_test[2][1] << "|" << res_test[2][2] << "|" << res_test[2][3] << "|" << std::log2(res_test[1][3]/res_test[2][3]) << "|\n";
    // out << "|" << 4 << "|" << int32_t(res_test[3][0]) << "|" << res_test[3][1] << "|" << res_test[3][2] << "|" << res_test[3][3] << "|" << std::log2(res_test[2][3]/res_test[3][3]) << "|\n";

    // out.close();
 

    FEMSolver fem_solver("Area1.txt", Test_2);
    fem_solver.DivideGrid(6);
    fem_solver.ClearAll(); /* Очистили внутреннее состояние  */
    fem_solver.PreCalc();  /* Подготовили слау */
    fem_solver.Calc(1000, 1e-10, true, false);

    SLAU_SparseMatrix slau_sparse = fem_solver.GetSparseSLAU();
    // SLAU_ProfileMatrix slau_profile;
    // //SparseMatrix2ProfileMatrix(slau_sparse.matrix, slau_profile.Matr);


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
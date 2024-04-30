#include <iostream>

#include "FEMSolver.h"
#include "Test.h"
#include "Integrate_3DStreightQuadPrismatic.h"
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

std::vector<double> Solve_TestTime(ParamDE Test, int32_t dividecoef = 1.0, int32_t maxiter = 15000, double eps = 1e-10, bool LU = true, bool LOS = true)
{
    cout << "\n";
    cout << "------------------------------------------------------------------\n";
    cout << "Solve start\n";

    std::vector<double> res(2);

    FEMSolver fem_solver("Area1.txt", Test);
    BaseGrid3D_Integrate_StreightQuadPrismatic IntegGrid = BaseGrid2BaseGridIntegrate(fem_solver.GetBaseGrid());

    fem_solver.DivideGrid(dividecoef);
    fem_solver.ClearAll(); /* Очистили внутреннее состояние  */
    fem_solver.PreCalc();  /* Подготовили слау */

    double total_time = fem_solver.Calc(maxiter, eps, false, false, LU, LOS);

    /* Функции для проверки сходимости */
    double t = 3.24;
    double fcos_calc = std::cos(Test.omega * t);
    double fsin_calc = std::sin(Test.omega * t);

    Func_3D f_norm_u_c = [&](double x, double y, double z) -> double
    { return std::pow((fem_solver.u_c(x, y, z) - Test.u_c(x, y, z)) * fcos_calc, 2.0); };
    Func_3D f_norm_u_s = [&](double x, double y, double z) -> double
    { return std::pow((fem_solver.u_s(x, y, z) - Test.u_s(x, y, z)) * fsin_calc, 2.0); };

    Integrate_3DStreightQuadPrismatic Integ_c(IntegGrid, f_norm_u_c);
    Integrate_3DStreightQuadPrismatic Integ_s(IntegGrid, f_norm_u_s);

    std::pair<double, double> fn_c = Integ_c.Quad("default", 1e-7, 4);
    std::pair<double, double> fn_s = Integ_s.Quad("default", 1e-7, 4);

    cout << "Finit Element Count = " << fem_solver.GetNodeCount() << "\n";
    cout << "LU mode:" << LU << "  LOS mode: " << LOS << "\n";
    cout << "f_norm_u_c = " << std::scientific << std::setprecision(4) << std::sqrt(fn_c.first) << " ErrorInteg = " << fn_c.second << "\n";
    cout << "f_norm_u_s = " << std::sqrt(fn_s.first) << " ErrorInteg = " << fn_s.second << "\n";
    cout << "General Error = " << std::sqrt(fn_s.first + fn_c.first) << "\n";
    cout << "------------------------------------------------------------------\n";

    res[0] = total_time;
    res[1] = std::sqrt(fn_s.first + fn_c.first);

    return res;
}

int main()
{

    // std::vector<int32_t> dividecoefs = {1, 2, 4, 8};

    // std::vector<std::vector<double>> res_test;
    // for(int32_t dc: dividecoefs)
    // {
    //     res_test.push_back(Solve(Test_1, dc, 100, 1e-10, false, false));
    // }

    // // Генерация табличных данных
    // cout << "Table\n";
    // cout << "|" << 1 << "|" << int32_t(res_test[0][0]) << "|" << res_test[0][1] << "|" << res_test[0][2] << "|" << res_test[0][3] << "|" << "-" << "|\n";
    // cout << "|" << 2 << "|" << int32_t(res_test[1][0]) << "|" << res_test[1][1] << "|" << res_test[1][2] << "|" << res_test[1][3] << "|" << std::log2(res_test[0][3]/res_test[1][3]) << "|\n";
    // cout << "|" << 3 << "|" << int32_t(res_test[2][0]) << "|" << res_test[2][1] << "|" << res_test[2][2] << "|" << res_test[2][3] << "|" << std::log2(res_test[1][3]/res_test[2][3]) << "|\n";
    // cout << "|" << 4 << "|" << int32_t(res_test[3][0]) << "|" << res_test[3][1] << "|" << res_test[3][2] << "|" << res_test[3][3] << "|" << std::log2(res_test[2][3]/res_test[3][3]) << "|\n";

    // std::vector<double> chi = {8.81 * 10e-12, 10e-11, 10e-10, 10e-9};

    // int idx = 4;

    // Test_2.ksi = chi[idx - 1];
    // // std::vector<double> res_LU = Solve_TestTime(Test_2, 6, 15000, 1e-10, true, false);
    // // cout << "|" << idx <<"|" << std::setprecision(2) << Test_2.omega  << "|" << Test_2.lambda  << "|" << Test_2.ksi << "|" << Test_2.sigma << "|" << std::fixed << std::setprecision(2) << std::round(res_LU[0]*100)/100 << "|"
    // //     << std::scientific << res_LU[1] << "|\n";

    std::vector<double> res_LOS = Solve_TestTime(Test_2, 6, 15000, 1e-10, false, true);
    cout << std::fixed << std::setprecision(2) << std::round(res_LOS[0] * 100) / 100 << "|" << std::scientific << res_LOS[1] << "|\n";

    return 0;
}
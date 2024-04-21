#include <iostream>

#include "FEMSolver.h"
#include "Test.h"
#include "Integrate_3DStreightQuadPrismatic.h"


typedef std::function<double(double, double, double)> Func_3D;

BaseGrid3D_Integrate_StreightQuadPrismatic BaseGrid2BaseGridIntegrate(const BaseGrid3DStreightQuadPrismatic &baseGrid)
{
    /* Тупо копируем все содержимое в структуру и все */
    BaseGrid3D_Integrate_StreightQuadPrismatic res;

    res.Nx = baseGrid.Nx;
    res.Ny = baseGrid.Ny;
    res.Nz = baseGrid.Nz;

    res.BaseGridXZ.resize(res.Nz);
    for(int32_t i = 0; i < res.Nz; i++)
        res.BaseGridXZ[i].resize(res.Nx);
    
    for(int32_t i = 0; i < baseGrid.Nz; i++)
    {
        for(int32_t j = 0; j < baseGrid.Nx; j++)
        {
            res.BaseGridXZ[i][j].x = baseGrid.BaseGridXZ[i][j].x;
            res.BaseGridXZ[i][j].z = baseGrid.BaseGridXZ[i][j].z;
        }
    }


    res.BaseGridY.resize(res.Ny);
    res.BaseGridY = baseGrid.BaseGridY;


    res.DivideParam.resize(3);
    res.DivideParam[0].resize(res.Nx-1);
    res.DivideParam[1].resize(res.Nz-1);
    res.DivideParam[2].resize(res.Ny-1);

    for(int32_t i = 0; i < baseGrid.DivideParam[0].size(); i++)
    {
        res.DivideParam[0][i].coef = baseGrid.DivideParam[0][i].coef;
        res.DivideParam[0][i].num = baseGrid.DivideParam[0][i].num;
    }
    

    for(int32_t i = 0; i < baseGrid.DivideParam[1].size(); i++)
    {
        res.DivideParam[1][i].coef = baseGrid.DivideParam[1][i].coef;
        res.DivideParam[1][i].num = baseGrid.DivideParam[1][i].num;
    }

    for(int32_t i = 0; i < baseGrid.DivideParam[2].size(); i++)
    {
        res.DivideParam[2][i].coef = baseGrid.DivideParam[2][i].coef;
        res.DivideParam[2][i].num = baseGrid.DivideParam[2][i].num;
    }

    res.isReadyToUse = true;

    return res;
}



ParamDE Test_1 = Test1();
ParamDE Test_2 = Test2();

void Solve(ParamDE Test,int32_t dividecoef = 1.0 ,int32_t maxiter = 100, double eps = 1e-10 ,bool printInterSlau = true, bool printNodeError = true)
{
    FEMSolver fem_solver("Area1.txt", Test);
    BaseGrid3D_Integrate_StreightQuadPrismatic IntegGrid = BaseGrid2BaseGridIntegrate(fem_solver.GetBaseGrid());

    fem_solver.DivideGrid(dividecoef);
    fem_solver.ClearAll(); /* Очистили внутреннее состояние  */
    fem_solver.PreCalc(); /* Подготовили слау */
    fem_solver.Calc(maxiter, eps, printInterSlau, printNodeError);

    /* Функции для проверки сходимости */
    Func_3D f_norm_u_c = [&](double x, double y, double z) -> double{ return std::pow(fem_solver.u_c(x,y,z)-Test.u_c(x,y,z) , 2.0); };
    Func_3D f_norm_u_s = [&](double x, double y, double z) -> double{ return std::pow(fem_solver.u_s(x,y,z)-Test.u_s(x,y,z) , 2.0); };
    
    Integrate_3DStreightQuadPrismatic Integ_c(IntegGrid, f_norm_u_c);
    Integrate_3DStreightQuadPrismatic Integ_s(IntegGrid, f_norm_u_s);
    
    std::pair<double, double> fn_c = Integ_c.Quad("default", 1e-7, 4);
    std::pair<double, double> fn_s = Integ_s.Quad("default", 1e-7, 4);
    cout << "\n";
    cout << "------------------------------------------------------------------\n";
    cout << "Finit Element Count = " << fem_solver.GetNodeCount() << "\n";
    cout << "f_norm_u_c = " << std::scientific << std::sqrt(fn_c.first) << " ErrorInteg = " << fn_c.second << "\n";
    cout << "f_norm_u_s = " << std::sqrt(fn_s.first) << " ErrorInteg = " << fn_s.second << "\n";
    cout << "General Error = " << std::sqrt(fn_s.first + fn_c.first) << "\n";
    cout << "------------------------------------------------------------------\n";
}

int main()
{

    std::vector<int32_t> dividecoefs = {1,2,4,8};

    for(int32_t dc: dividecoefs)
    {
        Solve(Test_2, dc, 100, 1e-10, false, false);
    }
        
   
   
    return 0;
}
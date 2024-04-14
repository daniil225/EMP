#ifndef FEM_SOLVER_H_
#define FEM_SOLVER_H_

#include "Grid3D_StreightQuadPrismatic.h"
#include "SLAU.h"
#include "Matrix.h"
#include <vector>
#include <tuple>
#include <functional>


struct LocalMatrix
{
private:
    std::vector<std::vector<double>> G1 = {  {1, -1},
                                             {-1, 1}};
    std::vector<std::vector<double>>M1 = {  { 2.0/6.0, 1.0/6.0 },
                                            { 1.0/6.0, 2.0/6.0} };

    std::function<int32_t(int32_t)> mu = [](int32_t i) -> int32_t { return (i-1)%2; };
    std::function<int32_t(int32_t)> vu = [](int32_t i) -> int32_t { return  int32_t((i-1)/2.0) % 2; };
    std::function<int32_t(int32_t)> th = [](int32_t i) -> int32_t { return  int32_t((i-1)/4); };
    double hx, hy, hz;

    /* Метод генерации локальной матрицы  */
    void GenerateLocalMatrix();
    
public:

    std::vector<std::vector<Block>> matrix;                  // Локальная блочная матрица
    std::vector<std::vector<std::pair<int32_t, int32_t>>> Local2GlobalIdx; // Матрица соответсвий локальной и глобальной нумерации


    LocalMatrix(Finit_Element_StreightQuadPrismatic &FE);
};

struct ParamDE
{
    double lambda = 0;
    double omega = 0;
    double sigma =  0;
    double ksi = 0;

    /* Составные части функции решения и правой части */
    /* f(x,y,z) */
    std::function<double (double, double, double) > u_s;
    std::function<double (double, double, double) > u_c;
    std::function<double (double, double, double) > f_s;
    std::function<double (double, double, double) > f_c;

    /* Функция решения и правой части f(x,y,z,t)*/
    std::function<double (double, double, double, double) > u;
    std::function<double (double, double, double, double) > f;


    /* Краевые условия 2-ого рода с указанием границы */
    
};

/* Генерация локальных матриц для глобальной матрицы */
class GenerateLocalMatrix
{
private:
public:
};

class FEMSolver
{
private:
    Grid3D_StreightQuadPrismatic Grid;
    SLAU_BlockMatrix slau_block;
    SLAU_ProfileMatrix slau_profile;
    SLAU_SparseMatrix slau_sparse;

public:

    FEMSolver() = default;

    FEMSolver(const std::string &filename);


    ~FEMSolver() = default;
};

#endif
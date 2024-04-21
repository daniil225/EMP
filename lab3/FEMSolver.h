#ifndef FEM_SOLVER_H_
#define FEM_SOLVER_H_

#include "Grid3D_StreightQuadPrismatic.h"
#include "SLAU.h"
#include "Matrix.h"
#include <vector>
#include <tuple>
#include <functional>



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
    std::vector<std::function<double(double, double, double)>> Theta_s;
    std::vector<std::function<double(double, double, double)>> Theta_c;
};


/* Генерация локальных матриц для глобальной матрицы */
class LocalMatrix
{
private:

    std::vector<std::vector<double>> G1 = {  {1, -1},
                                             {-1, 1}};
    std::vector<std::vector<double>>M1 = {  { 2.0/6.0, 1.0/6.0 },
                                            { 1.0/6.0, 2.0/6.0} };

    std::vector<std::vector<double>> M_2KU = {  {4.0/36.0,2.0/36.0,2.0/36.0,1.0/36.0},
                                                {2.0/36.0,4.0/36.0,1.0/36.0,2.0/36.0},
                                                {2.0/36.0,1.0/36.0,4.0/36.0,2.0/36.0},
                                                {1.0/36.0,2.0/36.0,2.0/36.0,4.0/36.0} };

    std::function<int32_t(int32_t)> mu = [](int32_t i) -> int32_t { return i%2; };  // x
    std::function<int32_t(int32_t)> vu = [](int32_t i) -> int32_t { return  int32_t(i/2.0) % 2; }; // z
    std::function<int32_t(int32_t)> th = [](int32_t i) -> int32_t { return  int32_t(i/4); }; // y
    double hx, hy, hz;

    Finit_Element_StreightQuadPrismatic FE;
    /* Метод генерации локальной матрицы  */
    void GenerateLocalMatrix();

    
public:

    std::vector<std::vector<Block>> matrix;                  // Локальная блочная матрица
    std::vector<std::vector<double>> f;                     //  Локальный вектор правой части 
    std::vector<std::vector<std::pair<int32_t, int32_t>>> Local2GlobalIdx; // Матрица соответсвий локальной и глобальной нумерации
   
    LocalMatrix() = default;
    LocalMatrix(const Finit_Element_StreightQuadPrismatic &FE);

    /* Очищаем все внутреннее состояние */
    void ClearInternalState();

    /* Расчет локальной матрицы + Устанавливает размер заново */
    void CalcLocalMatrix(const ParamDE &param);


    /* Устанавливаем конечный элемент */
    void SetFE(const Finit_Element_StreightQuadPrismatic &FE);

    friend ostream& operator<<(ostream &os, const LocalMatrix &lcmatr);

    ~LocalMatrix() = default;
};





class FEMSolver
{
private:
    Grid3D_StreightQuadPrismatic Grid; // Ну собственно сетка 
    SLAU_BlockMatrix slau_block; // СЛАУ с матрицей в блочном разряженном строчно-столбцовом формате 
    SLAU_ProfileMatrix slau_profile; // СЛАУ с матрицей в профильном формате 
    SLAU_SparseMatrix slau_sparse; // СЛАУ с матрицей в разряженном строчно столбцовом формате 
    ParamDE param; // Параметры ДУ

    std::vector<double> q; // Вектор решения 

    void GenerateSLAU();

public:

    FEMSolver() = default;

    FEMSolver(const std::string &filename, const ParamDE &param, bool startCalc = false);

    /* Подготовка к решению задачи. Включает построение портрета матрицы + генерация СЛАУ */
    void PreCalc();

    /* Расчет решения. Конкретно решение СЛАУ */
    void Calc(int32_t maxiter = 1000, double eps = 1e-15, bool showiterslau = true, bool showcalcfinnode = true);

    /* Дробление сетки */
    void DivideGrid(int32_t coef);

    /* Полная очистка внутреннего состояния. Включает только очистку струтктур для СЛАУ */
    void ClearAll();

    int32_t inline GetNodeCount() const {return Grid.GetGridSize().FEMCount; }


    /* Построенная функция */
    double u(double x, double y, double z, double t);

    double u_c(double x, double y, double z);
    double u_s(double x, double y, double z);

    /* Получить базовую сетку для интегрирования по области */
    inline BaseGrid3DStreightQuadPrismatic GetBaseGrid() { return Grid.GetBaseGrid(); }
    

    ~FEMSolver() = default;
};

#endif
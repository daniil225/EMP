#include "FEMSolver.h"
#include "GenerateMatrixPortrait.h"
#include "SlauSolver.h"
#include <tuple>
#include <iomanip>
#include <chrono>
/* LocalMatrix Class Realisation */
using namespace std;
/* Private */
void LocalMatrix::GenerateLocalMatrix()
{
}

ostream &operator<<(ostream &os, const LocalMatrix &lcmatr)
{
    os << "Matrix: \n";

    for (int32_t i = 0; i < 8; i++)
    {
        for (int32_t j = 0; j < 8; j++)
        {
            os << lcmatr.matrix[i][j].pij << " " << -lcmatr.matrix[i][j].cij << " ";
        }
        os << "\n";
    }

    for (int32_t i = 0; i < 8; i++)
    {
        for (int32_t j = 0; j < 8; j++)
        {
            os << lcmatr.matrix[i][j].cij << " " << lcmatr.matrix[i][j].pij << " ";
        }
        os << "\n";
    }

    os << "End Matrix\n\n";

    os << "f: \n";
    for (int i = 0; i < 8; i++)
    {
        os << lcmatr.f[i][0] << " " << lcmatr.f[i][1] << "\n";
    }
    os << "End f\n\n";

    os << "Local to Global: \n";

    for (int32_t i = 0; i < 8; i++)
    {
        for (int32_t j = 0; j < 8; j++)
        {
            os << "( " << i << ";" << j << ") => (" << lcmatr.Local2GlobalIdx[i][j].first
               << ";" << lcmatr.Local2GlobalIdx[i][j].second << ") ";
        }
        os << "\n";
    }
    os << "Local To Global end\n\n";

    return os;
}

/* Public */
LocalMatrix::LocalMatrix(const Finit_Element_StreightQuadPrismatic &FE) : FE(FE)
{
    /* Локальные матрицы 8*8 */

    matrix.resize(8);
    Local2GlobalIdx.resize(8);
    f.resize(8);
    for (int32_t i = 0; i < 8; i++)
    {
        matrix[i].resize(8);
        Local2GlobalIdx[i].resize(8);
        f.resize(2);
    }

    /* Определаем шаги */
    hx = FE.e[1].x - FE.e[0].x;
    hy = FE.e[4].y - FE.e[0].y;
    hz = FE.e[2].z - FE.e[0].z;
}

void LocalMatrix::ClearInternalState()
{

    if (!matrix.empty() && !f.empty() && !Local2GlobalIdx.empty())
    {
        for (int32_t i = 0; i < 8; i++)
        {
            matrix[i].clear();
            Local2GlobalIdx[i].clear();
            f[i].clear();
        }
        matrix.clear();
        Local2GlobalIdx.clear();
        f.clear();
    }
    /* Локальные матрицы 8*8 */

    matrix.resize(8);
    Local2GlobalIdx.resize(8);
    f.resize(8);
    for (int32_t i = 0; i < 8; i++)
    {
        matrix[i].resize(8);
        Local2GlobalIdx[i].resize(8);
        f[i].resize(2);
    }
}

/* Расчет локальной матрицы */
void LocalMatrix::CalcLocalMatrix(const ParamDE &param)
{
    /* В цикле по Точкам */
    for (int32_t i = 0; i < FE.FinitElementSize; i++)
    {
        /* Локальный вектор для правой части */
        double f_c = 0;
        double f_s = 0;
        /* Тут же учет 2-х КУ */
        for (int32_t j = 0; j < FE.FinitElementSize; j++)
        {
            double x, y, z;
            x = FE.e[j].x;
            y = FE.e[j].y;
            z = FE.e[j].z;

            Block block; // Получившийся блок
            /* Подготовливаем p */

            double p1 = (hy * hz / hx) * G1[mu(i)][mu(j)] * M1[vu(i)][vu(j)] * M1[th(i)][th(j)];
            double p2 = (hx * hz / hy) * M1[mu(i)][mu(j)] * G1[vu(i)][vu(j)] * M1[th(i)][th(j)];
            double p3 = (hx * hy / hz) * M1[mu(i)][mu(j)] * M1[vu(i)][vu(j)] * G1[th(i)][th(j)];
            /* При желании можно добавить 3-и КУ. Это все лишь добавить второй интеграл по поверхности */
            block.pij = param.lambda * (p1 + p2 + p3) - param.ksi * param.omega * param.omega * hx * hy * hz * M1[mu(i)][mu(j)] * M1[vu(i)][vu(j)] * M1[th(i)][th(j)];

            block.cij = param.omega * param.sigma * hx * hy * hz * M1[mu(i)][mu(j)] * M1[vu(i)][vu(j)] * M1[th(i)][th(j)];
            /* Заносим в матрицу */
            matrix[i][j] = block;

            f_c += param.f_c(x, y, z) * hx * hy * hz * M1[mu(i)][mu(j)] * M1[vu(i)][vu(j)] * M1[th(i)][th(j)];
            f_s += param.f_s(x, y, z) * hx * hy * hz * M1[mu(i)][mu(j)] * M1[vu(i)][vu(j)] * M1[th(i)][th(j)];

            /* Генерация матрицы соответствия локальной и глобальной матрицы */
            Local2GlobalIdx[i][j] = {FE.GlobalIdx[i], FE.GlobalIdx[j]};
        }

        f[i][0] = f_s;
        f[i][1] = f_c;
    }

    /* Учет 2-КУ */
    /* В цикле по границе КЭ */
    /* Шуруем по границе */
    for (int32_t idx = 0; idx < FE.BoundCount; idx++)
    {
        /* Если это хотя бы граница и ее тип 2-ой*/
        if (FE.Bound[idx].IsBound == 1 && FE.Bound[idx].BoundType == 2)
        {
            /* Обнуляем строки в матрице там где КУ 2 */
            int32_t formula = FE.Bound[idx].BoundInfo;
            double h = 0.0;
            if (formula == 0 || formula == 1)
                h = hx * hz;
            else if (formula == 2 || formula == 3)
                h = hy * hz; // +
            else if (formula == 4 || formula == 5)
                h = hx * hy;

            // cout << "formula = " << formula << " " << "h = " << h << "\n";
            //  Вставлаем вектор в правую часть
            for (int32_t i = 0; i < 4; i++)
            {
                double f_c = 0.0;
                double f_s = 0.0;
                for (int32_t j = 0; j < 4; j++)
                {
                    double x = FE.e[FE.Bound[idx].LocalIdx[j]].x;
                    double y = FE.e[FE.Bound[idx].LocalIdx[j]].y;
                    double z = FE.e[FE.Bound[idx].LocalIdx[j]].z;

                    f_s += h * param.Theta_s[formula](x, y, z) * M_2KU[i][j];
                    f_c += h * param.Theta_c[formula](x, y, z) * M_2KU[i][j];

                    double fs = param.Theta_s[formula](x, y, z);
                    double fc = param.Theta_c[formula](x, y, z);

                    // cout << param.Theta_s[formula](x,y,z) << " " << x << " " << y << " " << z << "\n";
                }
                // cout << "Add to local IDX = " << FE.Bound[idx].LocalIdx[i] << "\n";
                f[FE.Bound[idx].LocalIdx[i]][0] += f_s;
                f[FE.Bound[idx].LocalIdx[i]][1] += f_c;
                // cout << "\n";
            }

            // std::cout << formula << "Global IDX = " << FE.Bound[idx].GlobalIdx[0] << ", " <<
            //                      FE.Bound[idx].GlobalIdx[1] << ", " << FE.Bound[idx].GlobalIdx[2] << ", " << FE.Bound[idx].GlobalIdx[3]<< "\n";
        }
    }
}

/* Устанавливаем конечный элемент */
void LocalMatrix::SetFE(const Finit_Element_StreightQuadPrismatic &FE)
{
    ClearInternalState(); // Сброс состяния + перевыделение памяти (не самое лучшее рещение) !!!
    this->FE = FE;

    /* Определаем шаги */
    hx = FE.e[1].x - FE.e[0].x;
    hy = FE.e[4].y - FE.e[0].y;
    hz = FE.e[2].z - FE.e[0].z;
}

/* FEMSolver Realisation */

/* Private */
void FEMSolver::GenerateSLAU()
{

    LocalMatrix localMatr;

    for (int32_t i = 0; i < Grid.GetGridSize().FEMCount; i++)
    {
        localMatr.SetFE(Grid.GetElement(i));
        localMatr.CalcLocalMatrix(param);

        /* Вставлаем элементы в матрицу */
        int32_t ii = 0;
        for (auto row : localMatr.Local2GlobalIdx)
        {
            slau_block.b[row[0].first][0] += localMatr.f[ii][0];
            slau_block.b[row[0].first][1] += localMatr.f[ii][1];

            int32_t j = 0;
            for (auto e : row)
            {
                slau_block.matrix.AddBlock(localMatr.matrix[ii][j], e.first, e.second);
                j++;
            }
            ii++;
        }
    }

    /* Учет 1 - КУ */
    /* Ну стратегия простоя и тупая
        идем по КЭ и смотрим если у нас там КУ 1-рода то фигачим их
     */
    for (int32_t elmi = 0; elmi < Grid.GetGridSize().FEMCount; elmi++)
    {
        Finit_Element_StreightQuadPrismatic FE = Grid.GetElement(elmi);
        // FE.PrintElement();
        /* Шуруем по границе */
        for (int32_t i = 0; i < FE.BoundCount; i++)
        {
            /* Если это хотя бы граница и ее тип 1-ый*/
            if (FE.Bound[i].IsBound == 1 && FE.Bound[i].BoundType == 1)
            {
                /* Обнуляем строки в матрице там где КУ 1 */
                int32_t ii = 0;
                for (int32_t idx : FE.Bound[i].GlobalIdx)
                {
                    /*  Собственно по всей строке */
                    for (int32_t j = 0; j < slau_block.N; j++)
                    {
                        Block elem;
                        elem.cij = 0;
                        elem.pij = 0;
                        slau_block.matrix.InsertBlock(elem, idx, j);
                    }
                    /* А теперь на дигональ фигачим 1-ы */
                    Block elem;
                    elem.pij = 1.0;
                    elem.cij = 0.0;
                    slau_block.matrix.InsertBlock(elem, idx, idx);

                    /* Модификация вектора правой части */

                    double x = FE.e[FE.Bound[i].LocalIdx[ii]].x;
                    double y = FE.e[FE.Bound[i].LocalIdx[ii]].y;
                    double z = FE.e[FE.Bound[i].LocalIdx[ii]].z;
                    slau_block.b[idx][0] = param.u_s(x, y, z);
                    slau_block.b[idx][1] = param.u_c(x, y, z);
                    ii++;
                }
            }
        }
    }
}

/* Public  */
FEMSolver::FEMSolver(const std::string &filename, const ParamDE &param_, bool startCalc) : param(param_)
{

    /* Генерация сетки */
    GridStatus Status = Grid.Load(filename);
    if (Status.GetState() == State::OK)
        Status = Grid.GenerateGrid();
    else
        std::cout << Status.GetMsg();

    /* Начать ли вычислять СЛАУ  */
    if (startCalc)
    {
        ClearAll();
        PreCalc();
    }
}

void FEMSolver::ClearAll()
{
    if (slau_block.N > 0)
    {
        slau_block.matrix.di.clear();
        slau_block.matrix.ggl.clear();
        slau_block.matrix.ggu.clear();
        slau_block.matrix.ig.clear();
        slau_block.matrix.jg.clear();
        slau_block.b.clear();

        slau_block.matrix.N = -1;
        slau_block.matrix.size = -1;
        slau_block.size = -1;
        slau_block.N = -1;
    }

    if (slau_sparse.N > 0)
    {
        slau_sparse.matrix.di.clear();
        slau_sparse.matrix.ggl.clear();
        slau_sparse.matrix.ggu.clear();
        slau_sparse.matrix.ig.clear();
        slau_sparse.matrix.jg.clear();
        slau_sparse.b.clear();

        slau_sparse.matrix.N = -1;
        slau_sparse.matrix.size = -1;
        slau_sparse.size = -1;
        slau_sparse.N = -1;
    }

    if (slau_profile.N > 0)
    {
        slau_profile.Matr.di.clear();
        slau_profile.Matr.al.clear();
        slau_profile.Matr.au.clear();
        slau_profile.Matr.ia.clear();
        slau_profile.f.clear();
        slau_profile.Matr.N = -1;
        slau_profile.Matr.size = -1;
        slau_profile.size = -1;
        slau_profile.N = -1;
    }

    if (!q.empty())
        q.clear();
}

void FEMSolver::PreCalc()
{
    // cout << "Generate SLAU start...\n";
    // Grid.PrintGridSlice(0);
    /* Генерация портрета матрицы */
    GeneratePortrait(slau_block, Grid);
    slau_block.N = slau_block.matrix.N;
    slau_block.size = slau_block.matrix.size;
    /* Выделение памяти под СЛАУ */
    slau_block.matrix.ggl.resize(slau_block.size);
    slau_block.matrix.ggu.resize(slau_block.size);
    slau_block.matrix.di.resize(slau_block.N);
    slau_block.b.resize(slau_block.N);
    for (int32_t i = 0; i < slau_block.N; i++)
        slau_block.b[i].resize(2);
    /* Генерация СЛАУ */
    GenerateSLAU();

    // slau_block.matrix.PrintDenseMatrix();
    /* перегоняем в разреженный формат матрицы */
    BlockMatrix2SparseMatrix(slau_block.matrix, slau_sparse.matrix);
    // slau_sparse.matrix.PrintDenseMatrix();

    /* Перегон вектора правой части */
    slau_sparse.N = slau_sparse.matrix.di.size();
    slau_sparse.size = slau_sparse.matrix.ggu.size();

    slau_sparse.matrix.N = slau_sparse.N;
    slau_sparse.matrix.size = slau_sparse.size;

    slau_sparse.b.resize(slau_sparse.N);
    for (int32_t i = 0, j = 0; i < slau_block.N; i++)
    {
        slau_sparse.b[j] = slau_block.b[i][0];
        j++;
        slau_sparse.b[j] = slau_block.b[i][1];
        j++;
    }

    /* Перегоняем в профильный формат матрицы */
    //SparseMatrix2ProfileMatrix(slau_sparse.matrix, slau_profile.Matr);
    // slau_profile.N = slau_profile.Matr.N;
    // slau_profile.size = slau_profile.Matr.size;
    // slau_profile.f.resize(slau_profile.N);
    // slau_profile.f = slau_sparse.b;
    // slau_profile.x.resize(slau_profile.N);

    // cout << "Generate SLAU end\n\n";
    /* А теперь здесь нужно решить СЛАУ */
}

/* Расчет решения */
double FEMSolver::Calc(int32_t maxiter, double eps, bool showiterslau, bool showcalcfinnode, bool LU, bool LOS)
{
    // cout << "SLAU solve start...\n";
    /* Собственно решение СЛАУ */
    /* Для начала перегоним форматы */
    double res = 0.0;
    if (LOS)
    {
        SLAUSolvers::IterationSolvers::SLAU iterslau;
        SLAUSolvers::IterationSolvers::Allocate_Memory_SLAU(iterslau, slau_sparse.N, slau_sparse.size, SLAUSolvers::MatrixType::SPARSE);
        /* Теперь инициализация самой слау ну тут тупо кипирование */
        std::copy(slau_sparse.matrix.di.begin(), slau_sparse.matrix.di.end(), iterslau.matrix.di);
        std::copy(slau_sparse.matrix.ggl.begin(), slau_sparse.matrix.ggl.end(), iterslau.matrix.ggl);
        std::copy(slau_sparse.matrix.ggu.begin(), slau_sparse.matrix.ggu.end(), iterslau.matrix.ggu);
        std::copy(slau_sparse.matrix.ig.begin(), slau_sparse.matrix.ig.end(), iterslau.matrix.ig);
        std::copy(slau_sparse.matrix.jg.begin(), slau_sparse.matrix.jg.end(), iterslau.matrix.jg);
        std::copy(slau_sparse.b.begin(), slau_sparse.b.end(), iterslau.f);
        iterslau.maxiter = maxiter;
        iterslau.eps = eps;

        SLAUSolvers::IterationSolvers::SetSolveMode(iterslau, SLAUSolvers::IterationSolvers::Solvemode::LOS_NOSYMETRIC_LU_FACT);
        
        auto begin = std::chrono::steady_clock::now();
        SLAUSolvers::IterationSolvers::SolveSLAU(iterslau, showiterslau);
        auto end = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        res = elapsed_ms.count()/1000.0;
        q.resize(iterslau.N);

        for (int32_t i = 0; i < iterslau.N; i++)
            q[i] = iterslau.x[i];

        if (showcalcfinnode)
        {
            for (int32_t i = 0, j = 0; i < iterslau.N; i += 2, j++)
            {
                double x, y, z;
                x = Grid[j].x;
                y = Grid[j].y;
                z = Grid[j].z;
                double u_s = param.u_s(x, y, z);
                double u_c = param.u_c(x, y, z);
                cout << "Node num = " << j << fixed << std::setprecision(4) << " u_s - u_s* = " << this->u_s(x, y, z) - u_s << "  u_c - u_c* = " << this->u_c(x, y, z) - u_c << "\n";
            }
        }
    }

    if (LU)
    {
        auto begin = std::chrono::steady_clock::now();
        SolveSlau(slau_profile, slau_profile.x);
        auto end = std::chrono::steady_clock::now();
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
        res = elapsed_ms.count()/1000.0;
        q = slau_profile.x;
    }


    return res;
    /* Решение при помощи  LU решателя */

    // cout << "SLAU solve end\n\n";
    // SLAUSolvers::IterationSolvers::PrintX(iterslau);
}

/* Дробление сетки */
void FEMSolver::DivideGrid(int32_t coef)
{
    /* Деление сетки + перегерация */
    Grid.DivideGrid(coef);
    Grid.ReGenerateGrid();

    // Grid.PrintGridSlice(0);
}

/* Построенная функция */
double FEMSolver::u(double x, double y, double z, double t)
{
}

/* cos компонента решения */
double FEMSolver::u_c(double x, double y, double z)
{
    double res = 0.0;
    /* Генерация Базисных функций для представления решения */
    /* Всего 8 базисных функций по количеству узлов на КЭ */
    Finit_Element_StreightQuadPrismatic FE = Grid.GetElement(x, y, z);

    double x0 = FE.e[0].x;
    double x1 = FE.e[1].x;
    double hx = x1 - x0;

    double y0 = FE.e[0].y;
    double y1 = FE.e[4].y;
    double hy = y1 - y0;

    double z0 = FE.e[0].z;
    double z1 = FE.e[2].z;
    double hz = z1 - z0;

    auto X1 = [&](double x) -> double
    { return (x1 - x) / hx; };
    auto X2 = [&](double x) -> double
    { return (x - x0) / hx; };

    auto Y1 = [&](double y) -> double
    { return (y1 - y) / hy; };
    auto Y2 = [&](double y) -> double
    { return (y - y0) / hy; };

    auto Z1 = [&](double z) -> double
    { return (z1 - z) / hz; };
    auto Z2 = [&](double z) -> double
    { return (z - z0) / hz; };

    auto psi1 = [&](double x, double y, double z)
    { return X1(x) * Z1(z) * Y1(y); };
    auto psi2 = [&](double x, double y, double z)
    { return X2(x) * Z1(z) * Y1(y); };
    auto psi3 = [&](double x, double y, double z)
    { return X1(x) * Z2(z) * Y1(y); };
    auto psi4 = [&](double x, double y, double z)
    { return X2(x) * Z2(z) * Y1(y); };

    auto psi5 = [&](double x, double y, double z)
    { return X1(x) * Z1(z) * Y2(y); };
    auto psi6 = [&](double x, double y, double z)
    { return X2(x) * Z1(z) * Y2(y); };
    auto psi7 = [&](double x, double y, double z)
    { return X1(x) * Z2(z) * Y2(y); };
    auto psi8 = [&](double x, double y, double z)
    { return X2(x) * Z2(z) * Y2(y); };

    function<double(double, double, double)> PSI[8] = {psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8};

    /* само решение */
    // FE.PrintElement();
    for (int32_t i = 0; i < 8; i++)
        res += q[2 * FE.GlobalIdx[i] + 1] * PSI[i](x, y, z);

    return res;
}

/* sin компонента решения */
double FEMSolver::u_s(double x, double y, double z)
{
    double res = 0.0;
    /* Генерация Базисных функций для представления решения */
    /* Всего 8 базисных функций по количеству узлов на КЭ */
    Finit_Element_StreightQuadPrismatic FE = Grid.GetElement(x, y, z);

    double x0 = FE.e[0].x;
    double x1 = FE.e[1].x;
    double hx = x1 - x0;

    double y0 = FE.e[0].y;
    double y1 = FE.e[4].y;
    double hy = y1 - y0;

    double z0 = FE.e[0].z;
    double z1 = FE.e[2].z;
    double hz = z1 - z0;

    auto X1 = [&](double x) -> double
    { return (x1 - x) / hx; };
    auto X2 = [&](double x) -> double
    { return (x - x0) / hx; };

    auto Y1 = [&](double y) -> double
    { return (y1 - y) / hy; };
    auto Y2 = [&](double y) -> double
    { return (y - y0) / hy; };

    auto Z1 = [&](double z) -> double
    { return (z1 - z) / hz; };
    auto Z2 = [&](double z) -> double
    { return (z - z0) / hz; };

    auto psi1 = [&](double x, double y, double z)
    { return X1(x) * Z1(z) * Y1(y); };
    auto psi2 = [&](double x, double y, double z)
    { return X2(x) * Z1(z) * Y1(y); };
    auto psi3 = [&](double x, double y, double z)
    { return X1(x) * Z2(z) * Y1(y); };
    auto psi4 = [&](double x, double y, double z)
    { return X2(x) * Z2(z) * Y1(y); };

    auto psi5 = [&](double x, double y, double z)
    { return X1(x) * Z1(z) * Y2(y); };
    auto psi6 = [&](double x, double y, double z)
    { return X2(x) * Z1(z) * Y2(y); };
    auto psi7 = [&](double x, double y, double z)
    { return X1(x) * Z2(z) * Y2(y); };
    auto psi8 = [&](double x, double y, double z)
    { return X2(x) * Z2(z) * Y2(y); };

    double q_local[8]; // Коэфициенты разложения
    function<double(double, double, double)> PSI[8] = {psi1, psi2, psi3, psi4, psi5, psi6, psi7, psi8};

    /* само решение */
    // std::cout << "\n";
    for (int32_t i = 0; i < 8; i++)
    {
        // cout << "Psi(" << x << "," << y << "," << z << ") = " << PSI[i](x,y,z) << "\n";
        res += q[2 * FE.GlobalIdx[i]] * PSI[i](x, y, z);
    }
    return res;
}

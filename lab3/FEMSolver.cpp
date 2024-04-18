#include "FEMSolver.h"
#include "GenerateMatrixPortrait.h"
#include "SlauSolver.h"
#include <tuple>
#include <iomanip>
/* LocalMatrix Class Realisation */

/* Private */
void LocalMatrix::GenerateLocalMatrix()
{

}

ostream& operator<<(ostream &os, const LocalMatrix &lcmatr)
{
    os << "Matrix: \n";

    for(int32_t i = 0; i < 8; i++)
    {
        for(int32_t j = 0; j < 8; j++)
        {
            os << lcmatr.matrix[i][j].pij << " " << -lcmatr.matrix[i][j].cij << " ";
        }
        os << "\n";
    }


    for(int32_t i = 0; i < 8; i++)
    {
        for(int32_t j = 0; j < 8; j++)
        {
            os << lcmatr.matrix[i][j].cij << " " << lcmatr.matrix[i][j].pij << " ";
        }
        os << "\n";
    }

    os << "End Matrix\n\n";

    os << "f: \n";
    for(int i = 0; i < 8; i++)
    {
        os << lcmatr.f[i][0] << " " << lcmatr.f[i][1] << "\n";
    }
    os << "End f\n\n";

    os << "Local to Global: \n";

    for(int32_t i = 0; i < 8; i++)
    {
        for(int32_t j = 0; j < 8; j++)
        {
            os << "( " << i << ";" << j << ") => (" << lcmatr.Local2GlobalIdx[i][j].first
                <<";" << lcmatr.Local2GlobalIdx[i][j].second << ") ";
        }
        os << "\n";
    }
    os << "Local To Global end\n\n";

    return os;
}

/* Public */
LocalMatrix::LocalMatrix(const Finit_Element_StreightQuadPrismatic &FE):FE(FE)
{
    /* Локальные матрицы 8*8 */

    matrix.resize(8);
    Local2GlobalIdx.resize(8);
    f.resize(8);
    for(int32_t i = 0; i < 8; i++)
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

    if(!matrix.empty() && !f.empty() && !Local2GlobalIdx.empty())
    {
        for(int32_t i = 0; i < 8; i++)
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
    for(int32_t i = 0; i < 8; i++)
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
    for(int32_t i = 0; i < FE.FinitElementSize; i++)
    {
        /* Локальный вектор для правой части */
        double f_c = 0;
        double f_s = 0;
        /* Тут же учет 2-х КУ */
        for(int32_t j = 0; j < FE.FinitElementSize; j++)
        {
            double x, y, z;
            x = FE.e[j].x;
            y = FE.e[j].y;
            z = FE.e[j].z;

            Block block; // Получившийся блок 
            /* Подготовливаем p */
            
            double p1 = (hy*hz/hx)*G1[mu(i)][mu(j)]*M1[vu(i)][vu(j)]*M1[th(i)][th(j)];
            double p2 = (hx*hy/hz)*M1[mu(i)][mu(j)]*G1[vu(i)][vu(j)]*M1[th(i)][th(j)];
            double p3 = (hx*hz/hy)*M1[mu(i)][mu(j)]*M1[vu(i)][vu(j)]*G1[th(i)][th(j)];
            /* При желании можно добавить 3-и КУ. Это все лишь добавить второй интеграл по поверхности */
            block.pij = param.lambda*(p1+p2+p3) - param.ksi*param.omega*param.omega*hx*hy*hz*M1[mu(i)][mu(j)]*M1[vu(i)][vu(j)]*M1[th(i)][th(j)];
            
            block.cij = param.omega*param.sigma*hx*hy*hz*M1[mu(i)][mu(j)]*M1[vu(i)][vu(j)]*M1[th(i)][th(j)];
            /* Заносим в матрицу */
            matrix[i][j] = block;

            /* Вектор правой части. Здесь же и происходит добавление 2-х КУ*/
            f_c += param.f_c(x,y,z)*hx*hy*hz*M1[mu(i)][mu(j)]*M1[vu(i)][vu(j)]*M1[th(i)][th(j)];
            f_s += param.f_s(x,y,z)*hx*hy*hz*M1[mu(i)][mu(j)]*M1[vu(i)][vu(j)]*M1[th(i)][th(j)];

            /* Генерация матрицы соответствия локальной и глобальной матрицы */
            Local2GlobalIdx[i][j] = {FE.GlobalIdx[i],FE.GlobalIdx[j]};
        }

        f[i][0] = f_s;
        f[i][1] = f_c;
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

    for(int32_t i = 0; i < Grid.GetGridSize().FEMCount; i++)
    {
        localMatr.SetFE(Grid.GetElement(i));
        localMatr.CalcLocalMatrix(param);

        /* Вставлаем элементы в матрицу */
        int32_t ii = 0;
        for(auto row: localMatr.Local2GlobalIdx)
        {
            slau_block.b[row[0].first][0] += localMatr.f[ii][0];
            slau_block.b[row[0].first][1] += localMatr.f[ii][1];

            int32_t j = 0;
            for(auto e: row)
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
    for(int32_t elmi = 0; elmi < Grid.GetGridSize().FEMCount; elmi++)
    {
        Finit_Element_StreightQuadPrismatic FE = Grid.GetElement(elmi);
        //FE.PrintElement();
        /* Шуруем по границе */
        for(int32_t i = 0; i < FE.BoundCount; i++)
        {
            /* Если это хотя бы граница и ее тип 1-ый*/
            if(FE.Bound[i].IsBound == 1 && FE.Bound[i].BoundType == 1)
            {
                /* Обнуляем строки в матрице там где КУ 1 */
                int32_t ii = 0;
                for(int32_t idx: FE.Bound[i].GlobalIdx)
                {
                    /*  Собственно по всей строке */
                    for(int32_t j = 0; j < slau_block.N; j++)
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
                    slau_block.b[idx][0] = param.u_s(x,y,z);
                    slau_block.b[idx][1] = param.u_c(x,y,z);
                    ii++;
                }
            }
        }
    }
}


/* Public  */
FEMSolver::FEMSolver(const std::string &filename, const ParamDE &param_): param(param_)
{

    /* Генерация сетки */
    GridStatus Status = Grid.Load(filename);
    if(Status.GetState() == State::OK) Status = Grid.GenerateGrid();
    else cout << Status.GetMsg();

    Grid.PrintGridSlice(0);
    /* Генерация портрета матрицы */
    GeneratePortrait(slau_block, Grid);
    slau_block.N = slau_block.matrix.N;
    slau_block.size = slau_block.matrix.size;
    /* Выделение памяти под СЛАУ */
    slau_block.matrix.ggl.resize(slau_block.size);
    slau_block.matrix.ggu.resize(slau_block.size);
    slau_block.matrix.di.resize(slau_block.N);
    slau_block.b.resize(slau_block.N);
    for(int32_t i = 0; i < slau_block.N; i++)
        slau_block.b[i].resize(2);
    /* Генерация СЛАУ */
    GenerateSLAU();

    slau_block.matrix.PrintDenseMatrix();
    BlockMatrix2SparseMatrix(slau_block.matrix, slau_sparse.matrix);
    slau_sparse.matrix.PrintDenseMatrix();

    /* Перегон вектора правой части */
    slau_sparse.N = slau_sparse.matrix.di.size();
    slau_sparse.size = slau_sparse.matrix.ggu.size();
    
    slau_sparse.matrix.N = slau_sparse.N;
    slau_sparse.matrix.size = slau_sparse.size;

    slau_sparse.b.resize(slau_sparse.N);
    for(int32_t i = 0, j = 0; i < slau_block.N; i++)
    {
        slau_sparse.b[j] = slau_block.b[i][0];
        j++;
        slau_sparse.b[j] = slau_block.b[i][1];
        j++;
    }
    /* А теперь здесь нужно решить СЛАУ */


}


/* Расчет решения */
void FEMSolver::Calc()
{
    /* Собственно решение СЛАУ */
    /* Для начала перегоним форматы */
    SLAUSolvers::IterationSolvers::SLAU iterslau;
    SLAUSolvers::IterationSolvers::Allocate_Memory_SLAU(iterslau,slau_sparse.N, slau_sparse.size,SLAUSolvers::MatrixType::SPARSE);
    /* Теперь инициализация самой слау ну тут тупо кипирование */
    std::copy(slau_sparse.matrix.di.begin(), slau_sparse.matrix.di.end(), iterslau.matrix.di);
    std::copy(slau_sparse.matrix.ggl.begin(), slau_sparse.matrix.ggl.end(),iterslau.matrix.ggl);
    std::copy(slau_sparse.matrix.ggu.begin(), slau_sparse.matrix.ggu.end(),iterslau.matrix.ggu);
    std::copy(slau_sparse.matrix.ig.begin(), slau_sparse.matrix.ig.end(),iterslau.matrix.ig);
    std::copy(slau_sparse.matrix.jg.begin(), slau_sparse.matrix.jg.end(),iterslau.matrix.jg);
    std::copy(slau_sparse.b.begin(), slau_sparse.b.end(), iterslau.f);
    iterslau.maxiter = 100;
    iterslau.eps = 1e-8;

    SLAUSolvers::IterationSolvers::SetSolveMode(iterslau, SLAUSolvers::IterationSolvers::Solvemode::LOS_NOSYMETRIC_DIAG_FACT);
    SLAUSolvers::IterationSolvers::SolveSLAU(iterslau, true);

    //SLAUSolvers::IterationSolvers::PrintX(iterslau);
    for(int32_t i = 0, j = 0; i < iterslau.N; i+=2, j++)
    {
        double x, y, z;
        x = Grid[j].x;
        y = Grid[j].y;
        z = Grid[j].z;
        double u_s = param.u_s(x,y, z);
        double u_c = param.u_c(x, y, z);
        cout << "Node num = " << j  << fixed << std::setprecision(4) << " u_s - u_s* = " << iterslau.x[i] - u_s << "  u_c - u_c* = " << iterslau.x[i+1] - u_c << "\n";
    }

    
}

/* Дробление сетки */
void FEMSolver::DivideGrid(int32_t coef)
{   
    /* Деление сетки + перегерация */
    Grid.DivideGrid(coef);
    Grid.ReGenerateGrid();

    //Grid.PrintGridSlice(0);
}


/* Построенная функция */
double FEMSolver::u(double x, double y, double z, double t)
{

}


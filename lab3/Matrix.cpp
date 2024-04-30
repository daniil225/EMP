#include "Matrix.h"
#include <iomanip>
#include <fstream>

SparseBlockOfMatrix::SparseBlockOfMatrix(int32_t N, int32_t size)
{
    this->N = N;
    this->size = size;
    ggu.resize(size);
    ggl.resize(size);
    di.resize(N);
    ig.resize(N + 1);
    jg.resize(size);
}

void SparseBlockOfMatrix::AllocMemory(int32_t N, int32_t size)
{
    this->N = N;
    this->size = size;
    ggu.resize(size);
    ggl.resize(size);
    di.resize(N);
    ig.resize(N + 1);
    jg.resize(size);
}

void SparseBlockOfMatrix::PrintDenseMatrix()
{
    std::ofstream out("SparseBlockMatrix.txt");

    out << "Dense Block Matrix\n";
    for (int32_t i = 0; i < N; i++)
    {

        for (int32_t j = 0; j < N; j++)
        {
            Block elem = GetBlockMatrix(*this, i, j);
            out << std::left << std::fixed << std::setprecision(3) << elem.pij << " " << -elem.cij << "   ";
        }
        out << "\n";

        for (int32_t j = 0; j < N; j++)
        {
            Block elem = GetBlockMatrix(*this, i, j);
            out << std::left << std::fixed << std::setprecision(3) << elem.cij << " " << elem.pij << "   ";
        }

        out << "\n\n";
    }
    out.close();
}

void SparseBlockOfMatrix::InsertBlock(Block &elem, int32_t i, int32_t j)
{
    if (i == j)
    {
        di[i] = elem;
    }

    if (i < j)
    {
        int ind;
        int igj0 = ig[j];
        int igj1 = ig[j + 1];
        for (ind = igj0; ind < igj1; ind++)
            if (jg[ind] == i)
                ggu[ind] = elem;
    }
    else
    {
        int ind;
        int igi0 = ig[i];
        int igi1 = ig[i + 1];
        for (ind = igi0; ind < igi1; ind++)
            if (jg[ind] == j)
                ggl[ind] = elem;
    }
}

void SparseBlockOfMatrix::AddBlock(Block &elem, int32_t i, int32_t j)
{
    if (i == j)
    {
        di[i].cij += elem.cij;
        di[i].pij += elem.pij;
    }

    if (i < j)
    {
        int ind;
        int igj0 = ig[j];
        int igj1 = ig[j + 1];
        for (ind = igj0; ind < igj1; ind++)
        {
            if (jg[ind] == i)
            {
                ggu[ind].cij += elem.cij;
                ggu[ind].pij += elem.pij;
            }
        }
    }
    else
    {
        int ind;
        int igi0 = ig[i];
        int igi1 = ig[i + 1];
        for (ind = igi0; ind < igi1; ind++)
        {
            if (jg[ind] == j)
            {
                ggl[ind].cij += elem.cij;
                ggl[ind].pij += elem.pij;
            }
        }
    }
}

bool SparseBlockOfMatrix::IsBlockInMatrix(int32_t i, int32_t j)
{
    if (i == j)
    {
        return true;
    }

    if (i < j)
    {
        int ind;
        int igj0 = ig[j];
        int igj1 = ig[j + 1];
        for (ind = igj0; ind < igj1; ind++)
            if (jg[ind] == i)
                return true;
    }
    else
    {
        int ind;
        int igi0 = ig[i];
        int igi1 = ig[i + 1];
        for (ind = igi0; ind < igi1; ind++)
            if (jg[ind] == j)
                return true;
    }

    return false;
}

SparseMatrix::SparseMatrix(int32_t N, int32_t size)
{
    this->N = N;
    this->size = size;
    ggu.resize(size);
    ggl.resize(size);
    di.resize(N);
    ig.resize(N + 1);
    jg.resize(size);
}

void SparseMatrix::AllocMemory(int32_t N, int32_t size)
{
    this->N = N;
    this->size = size;
    ggu.resize(size);
    ggl.resize(size);
    di.resize(N);
    ig.resize(N + 1);
    jg.resize(size);
}

void SparseMatrix::PrintDenseMatrix()
{
    std::ofstream out("SparseMatrix.txt");
    out << "Dense Sparse Matrix\n";
    for (int32_t i = 0; i < N; i++)
    {

        for (int32_t j = 0; j < N; j++)
        {
            if(j % 2 == 0 && j > 0) out << "   ";
            double elem = GetElementSparsematrix(*this, i, j);
            out << std::left << elem << " ";
        }
        if(i % 2 != 0) out << "\n\n";
        else out << "\n";
    }
    out.close();
}

void SparseMatrix::InsertElem(double elem, int32_t i, int32_t j)
{
    if (i == j)
    {
        di[i] = elem;
    }

    if (i < j)
    {
        int ind;
        int igj0 = ig[j];
        int igj1 = ig[j + 1];
        for (ind = igj0; ind < igj1; ind++)
            if (jg[ind] == i)
                ggu[ind] = elem;
    }
    else
    {
        int ind;
        int igi0 = ig[i];
        int igi1 = ig[i + 1];
        for (ind = igi0; ind < igi1; ind++)
            if (jg[ind] == j)
                ggl[ind] = elem;
    }
}

Block GetBlockMatrix(SparseBlockOfMatrix &Matr, int32_t i, int32_t j)
{
    std::vector<Block> &di = Matr.di;
    std::vector<Block> &ggu = Matr.ggu;
    std::vector<Block> &ggl = Matr.ggl;
    std::vector<int32_t> &ig = Matr.ig;
    std::vector<int32_t> &jg = Matr.jg;
    int N = Matr.N;

    if (i == j)
    {
        return di[i];
    }

    if (i < j)
    {
        int ind;
        int igj0 = ig[j];
        int igj1 = ig[j + 1];
        for (ind = igj0; ind < igj1; ind++)
            if (jg[ind] == i)
                return ggu[ind];
    }
    else
    {
        int ind;
        int igi0 = ig[i];
        int igi1 = ig[i + 1];
        for (ind = igi0; ind < igi1; ind++)
            if (jg[ind] == j)
                return ggl[ind];
    }

    return Block();
}

double GetElementSparsematrix(SparseMatrix &Matr, int32_t i, int32_t j)
{
    std::vector<double> &di = Matr.di;
    std::vector<double> &ggu = Matr.ggu;
    std::vector<double> &ggl = Matr.ggl;
    std::vector<int32_t> &ig = Matr.ig;
    std::vector<int32_t> &jg = Matr.jg;
    int N = Matr.N;

    if (i == j)
    {
        return di[i];
    }

    if (i < j)
    {
        int ind;
        int igj0 = ig[j];
        int igj1 = ig[j + 1];
        for (ind = igj0; ind < igj1; ind++)
            if (jg[ind] == i)
                return ggu[ind];
    }
    else
    {
        int ind;
        int igi0 = ig[i];
        int igi1 = ig[i + 1];
        for (ind = igi0; ind < igi1; ind++)
            if (jg[ind] == j)
                return ggl[ind];
    }

    return 0;
}



ProfileMatrix::ProfileMatrix(int32_t N, int32_t size)
{
    this->N = N;
    this->size = size;
    
    ia.resize(N+1);
    au.resize(size);
    al.resize(size);
    di.resize(N);
}

void ProfileMatrix::AllocateMemory(int32_t N, int32_t size)
{
    this->N = N;
    this->size = size;
    
    ia.resize(N+1);
    au.resize(size);
    al.resize(size);
    di.resize(N);
}

void ProfileMatrix::InsertElem(double elem, int32_t i, int32_t j)
{
    if(i == j)
    {
        di[i] = elem;
    }

    if(i < j)
    {
        int32_t ind = ia[j+1] - (j-i);
        au[ind] = elem;
    }
    else
    {
        int32_t ind = ia[i+1] - (i-j);
        al[ind] = elem;
    }
}


void BlockMatrix2SparseMatrix(SparseBlockOfMatrix &src_Matr, SparseMatrix &dst_Matr)
{
    int32_t size_row1 = 0;
    int32_t size_row2 = 0;
    int32_t idxjg = 0; // Индексация для jg не Блочной 
    int32_t idxjg_ggu = 0; //  Индексация для массива ggu не Блочный 
    int32_t idxdi = 0; // Индекс для диагонали разряженной матрицы не Блочной
    dst_Matr.N = 2 * src_Matr.N;
    dst_Matr.ig.resize(dst_Matr.N + 1);
    dst_Matr.ig[0] = 0;

    for (int32_t i = 0, k = 1; i < src_Matr.N; i++, k += 2)
    {
        size_row1 = 2 * (src_Matr.ig[i + 1] - src_Matr.ig[i]);
        size_row2 = size_row1 + 1;
        dst_Matr.ig[k] = dst_Matr.ig[k - 1] + size_row1;
        dst_Matr.ig[k + 1] = dst_Matr.ig[k] + size_row2;
    }
    
    /* В цикле по всем элементам в Блочной матрице */

    dst_Matr.jg.resize(dst_Matr.ig[dst_Matr.N]);
    dst_Matr.ggl.resize(dst_Matr.ig[dst_Matr.N]);
    dst_Matr.ggu.resize(dst_Matr.ig[dst_Matr.N]);
    dst_Matr.di.resize(dst_Matr.N);

    for (int32_t i = 0; i < src_Matr.N; i++)
    {
        /* Цикл для верхней строки */
        for (int32_t j = 0; j < src_Matr.N; j++)
        {
            /* Проверить что элемент с индексами i и j есть в матрице */
            if(src_Matr.IsBlockInMatrix(i,j))
            {
                Block elem = GetBlockMatrix(src_Matr, i, j);
                /* Вне диагональные элементы */
                if(i > j)
                {
                    dst_Matr.jg[idxjg] = 2*j;
                    dst_Matr.ggl[idxjg] = elem.pij;
                    idxjg++;
                    dst_Matr.jg[idxjg] = 2*j + 1;
                    dst_Matr.ggl[idxjg] = -elem.cij; /// Не забуль минуса
                    idxjg++;
                }
                else if(i == j)
                {
                    /* Диагональ */
                    dst_Matr.di[idxdi] = elem.pij;
                    idxdi++;
                }
            }

            /* Верхная часть матрицы */
            if(src_Matr.IsBlockInMatrix(j, i))
            {
                Block elem = GetBlockMatrix(src_Matr, j, i);

                if(i > j)
                {
                    dst_Matr.ggu[idxjg_ggu] = elem.pij;
                    idxjg_ggu++;
                    dst_Matr.ggu[idxjg_ggu] = elem.cij;
                    idxjg_ggu++;

                    //std::cout << elem.pij << " " << elem.cij << " ";
                }
            }
        }

        //std::cout << "\n";
        /* Цикл для нижней строки */
        for(int32_t j = 0; j < src_Matr.N; j++)
        {
            /* Проверить что элемент с индексами i и j есть в матрице */
            if(src_Matr.IsBlockInMatrix(i,j))
            {
                Block elem = GetBlockMatrix(src_Matr, i, j);
                /* Вне диагональные элементы */
                if(i > j)
                {
                    dst_Matr.jg[idxjg] = 2*j;
                    dst_Matr.ggl[idxjg] = elem.cij;
                    idxjg++;

                    dst_Matr.jg[idxjg] = 2*j + 1;
                    dst_Matr.ggl[idxjg] = elem.pij;
                    idxjg++;
                }
                else if(i == j)
                {
                    dst_Matr.jg[idxjg] = 2*j;
                    dst_Matr.ggl[idxjg] = elem.cij;
                    idxjg++;

                    /* Диагональ */
                    dst_Matr.di[idxdi] = elem.pij;
                    idxdi++;
                }
            }

            if(src_Matr.IsBlockInMatrix(j,i))
            {
                Block elem = GetBlockMatrix(src_Matr, j, i);
                if(i > j)
                {

                    dst_Matr.ggu[idxjg_ggu] = -elem.cij; // Не забудь минус
                    idxjg_ggu++;
                    dst_Matr.ggu[idxjg_ggu] = elem.pij;
                    idxjg_ggu++;
                }
                else if(i == j)
                {
                    dst_Matr.ggu[idxjg_ggu] = -elem.cij; // Не забудь минуса 
                    idxjg_ggu++;
                }
            }
        }

    }
}

void SparseMatrix2ProfileMatrix(SparseMatrix &src_Matr, ProfileMatrix& dst_Matr)
{
    /* В цикле по элементам */
    dst_Matr.N = src_Matr.N;

    dst_Matr.ia.resize(src_Matr.N+1);
    dst_Matr.ia[0] = 0;
    //dst_Matr.ia[1] = 0;

    for(int32_t i = 0; i < src_Matr.N; i++)
    {
        int32_t cnt_elem = 0;
        for(int32_t j = 0; j < i; j++)
        {
            double elem = GetElementSparsematrix(src_Matr, i, j);
            if(std::abs(elem) >= 0)
            {
                cnt_elem = i - j;
                break;
            }
        }
        dst_Matr.ia[i+1] = dst_Matr.ia[i] + cnt_elem;

        
    }

    // Выделение памяти
    dst_Matr.size = dst_Matr.ia[dst_Matr.N];
    dst_Matr.au.resize(dst_Matr.size);
    dst_Matr.al.resize(dst_Matr.size);
    dst_Matr.di.resize(dst_Matr.N);


    for(int32_t i = 0; i < src_Matr.N; i++)
    {
        int32_t cnt_elem = 0;
        for(int32_t j = 0; j < src_Matr.N; j++)
        {
            double elem = GetElementSparsematrix(src_Matr, i, j);

            if(std::abs(elem) > 0)
                dst_Matr.InsertElem(elem, i, j);
        }
    }

}
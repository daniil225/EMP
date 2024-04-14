#include "Matrix.h"
#include <iomanip>

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
    std::cout << "Dense Block Matrix\n";
    for (int32_t i = 0; i < N; i++)
    {

        for (int32_t j = 0; j < N; j++)
        {
            Block elem = GetBlockMatrix(*this, i, j);
            std::cout << std::left << elem.pij << " " << -elem.cij << "   ";
        }
        std::cout << "\n";

        for (int32_t j = 0; j < N; j++)
        {
            Block elem = GetBlockMatrix(*this, i, j);
            std::cout << std::left << elem.cij << " " << elem.pij << "   ";
        }

        std::cout << "\n\n";
    }
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
    std::cout << "Dense Sparse Matrix\n";
    for (int32_t i = 0; i < N; i++)
    {

        for (int32_t j = 0; j < N; j++)
        {
            if(j % 2 == 0 && j > 0) std::cout << "   ";
            double elem = GetElementSparsematrix(*this, i, j);
            std::cout << std::left << elem << " ";
        }
        if(i % 2 != 0) std::cout << "\n\n";
        else std::cout << "\n";
    }
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
#include "GenerateMatrixPortrait.h"
#include <cmath>

void GeneratePortrait(SLAU_BlockMatrix &slau, const Grid3D_StreightQuadPrismatic  &grid)
{
    Grid3D_Size gridParam = grid.GetGridSize();
    int Nx = gridParam.GlobalNx;
    int Nz = gridParam.GlobalNz;
    /* По номеру конечного элемента и локальному номеру на нем возвращает глобальный индекс конечного элемента
     * @param
     * int ielem - номер конечного элемента
     * int j - локальный номер
     */
    auto IndexOfUnknown = [&](int idx, int j) -> int
    {
        int shiftXZ = 0;
        int projidx = idx % (((Nx - 1) * (Nz - 1)));
        if (projidx < Nx - 1)
        {
            shiftXZ = projidx;
        }
        else
        {
            int level = std::floor(projidx / (Nx - 1));
            int start = level * Nx;
            int shift = projidx - (Nx - 1) * level;
            shiftXZ = start + shift;
        }

        int levelY = std::floor((idx) / ((Nx - 1) * (Nz - 1)));
        int shiftY = levelY * Nx * Nz;

        int res = shiftY + shiftXZ;

        switch (j)
        {
        case 0:
            res += 0;
            break;
        case 1:
            res += 1;
            break;
        case 2:
            res += Nx;
            break;
        case 3:
            res += (Nx + 1);
            break;
        case 4:
            res += (Nx * Nz);
            break;
        case 5:
            res += (Nx * Nz + 1);
            break;
        case 6:
            res += (Nx * Nz) + Nx;
            break;
        case 7:
            res += (Nx * Nz) + Nx + 1;
            break;
        default:
            break;
        }

        return res;
    };

    int globalN = gridParam.Dim;
    int regionsNum = gridParam.FEMCount;

    slau.matrix.ig.clear();
    slau.matrix.ig.resize(globalN + 1);
    std::vector<int32_t> &ig = slau.matrix.ig;

    int *list[2]{};
    list[0] = new int[2 * globalN * (globalN - 2)];
    list[1] = new int[2 * globalN * (globalN - 2)];
    int *listbeg = new int[globalN];

    int listsize = -1;
    for (int i = 0; i < globalN; i++)
        listbeg[i] = -1;
    for (int ielem = 0; ielem < regionsNum; ielem++)
    {
        for (int i = 0; i < 8; i++)
        {
            int k = IndexOfUnknown(ielem, i);
            for (int j = i + 1; j < 8; j++)
            {
                int ind1 = k;
                int ind2 = IndexOfUnknown(ielem, j);
                if (ind2 < ind1)
                {
                    ind1 = ind2;
                    ind2 = k;
                }
                int iaddr = listbeg[ind2];
                if (iaddr == -1)
                {
                    listsize++;
                    listbeg[ind2] = listsize;
                    list[0][listsize] = ind1;
                    list[1][listsize] = -1;
                }
                else
                {
                    while (list[0][iaddr] < ind1 && list[1][iaddr] > 0)
                        iaddr = list[1][iaddr];

                    if (list[0][iaddr] > ind1)
                    {
                        listsize++;
                        list[0][listsize] = list[0][iaddr];
                        list[1][listsize] = list[1][iaddr];
                        list[0][iaddr] = ind1;
                        list[1][iaddr] = listsize;
                    }
                    else
                    {
                        if (list[0][iaddr] < ind1)
                        {
                            listsize++;
                            list[1][iaddr] = listsize;
                            list[0][listsize] = ind1;
                            list[1][listsize] = -1;
                        }
                    }
                }
            }
        }
    }

    slau.matrix.jg.clear();
    slau.matrix.jg.resize(listsize + 1);
    std::vector<int32_t>& jg = slau.matrix.jg;

    ig[0] = 0;
    for (int i = 0; i < globalN; i++)
    {
        ig[i + 1] = ig[i];
        int iaddr = listbeg[i];

        while (iaddr != -1)
        {
            jg[ig[i + 1]] = list[0][iaddr];
            ig[i + 1]++;
            iaddr = list[1][iaddr];
        }
    }
    slau.N = globalN;
    slau.size = slau.matrix.ig[globalN];
    slau.matrix.N = globalN;
    slau.matrix.size = slau.size;

    /* Очистка памяти */
    delete[] listbeg;
    delete[] list[0];
    delete[] list[1];
}
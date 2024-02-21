#include "Slau.h"

void LUDecomposition(Slau &slau)
{
    Matrix &Matr = slau.Matr;

    const int n = Matr.di.size();
    const int TriangMat = Matr.ia.size();

    for (int i = 0; i < n; i++)
    {
        int i0 = Matr.ia[i];
        int i1 = Matr.ia[i + 1];

        int j = i - (i1 - i0);
        double sd = 0;
        for (int m = i0; m < i1; m++, j++)
        {
            double sl = 0;
            double su = 0;

            int j0 = Matr.ia[j];
            int j1 = Matr.ia[j + 1];

            int mi = i0;
            int mj = j0;

            int kol_i = m - i0;
            int kol_j = j1 - j0;
            int kol_r = kol_i - kol_j;

            if (kol_r < 0)
                mj -= kol_r;
            else
                mi += kol_r;

            for (; mi < m; mi++, mj++)
            {
                sl += Matr.al[mi] * Matr.au[mj];
                su += Matr.au[mi] * Matr.al[mj];
            }

            Matr.au[m] = Matr.au[m] - su;
            Matr.al[m] = (Matr.al[m] - sl) / Matr.di[j];

            sd += Matr.al[m] * Matr.au[m];
        }
        Matr.di[i] = Matr.di[i] - sd;
    }
}

void GausForward(Slau &slau, vector<double> &y)
{
    Matrix &Matr = slau.Matr;
    int n = Matr.di.size();

    // Решение системы Ly = b
    for (int i = 0; i < n; i++)
    {
        auto &al = Matr.al;
        int i0 = Matr.ia[i];
        int i1 = Matr.ia[i + 1];
        double s = 0;
        int j = i - (i1 - i0);
        for (int k = i0; k < i1; k++, j++)
        {
            s += al[k] * y[j];
        }
        y[i] = slau.f[i] - s;
    }
}

void GausBack(Slau &slau, vector<double> &x)
{
    Matrix &Matr = slau.Matr;
    int n = Matr.di.size();
    //  Решение системы Ux = y
    for (int i = n - 1; i >= 0; i--)
    {
        double xi = x[i] / Matr.di[i];
        auto &au = Matr.au;
        int i0 = Matr.ia[i];
        int i1 = Matr.ia[i + 1];
        // int m = Matr.ai[i+1] - Matr.ai[i];
        int j = i - 1;
        for (int k = i1 - 1; k >= i0; k--, j--)
        {
            x[j] -= au[k] * xi;
        }
        x[i] = xi;
    }
}

void SolveSlau(Slau &slau, vector<double> &x)
{
    LUDecomposition(slau);
    GausForward(slau, x);
    GausBack(slau, x);
}

void TestSlau()
{
    vector<double> q;
    Slau slau;
    slau.Matr.di = { 1, 2.66667, 2.66667, 2.66667, 2.66667, 2.66667, 2.66667, 2.66667, 2.66667, 1 };
	slau.Matr.al = { -0.833333, -0.833333, -0.833333, -0.833333, -0.833333, -0.833333, -0.833333, -0.833333, 0 };
	slau.Matr.au = { 0, -0.833333, -0.833333, -0.833333, -0.833333, -0.833333, -0.833333, -0.833333, -0.833333 };
    slau.Matr.ia = { 0,0,1,2,3,4,5,6,7,8,9 };
	slau.f = { 1, 2.000008, 3.000012, 4.000016, 5.00002, 6.000024, 7.000028, 8.000032, 9.000036, 10 };
    q.resize(10, 0); // выходной вектор 1,2,3,4,5,6,7,8,9,10
	SolveSlau(slau, q);
	//cout << q << endl;
}
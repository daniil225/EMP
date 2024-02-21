#include "Matrix.h"

vector<double> MultAOnq(Matrix &matr, vector<double>& q)
{
	vector1D tmp;
	tmp.resize(matr.di.size());

	if (matr.di.size() >= 2)
		tmp[0] = matr.di[0] * q[0] + matr.au[0] * q[1];

	if (matr.di.size() >= 3)
		for (size_t i = 1; i < matr.di.size() - 1; i++)
			tmp[i] = matr.al[i - 1] * q[i - 1] + matr.di[i] * q[i] + matr.au[i] * q[i + 1];

	int lIndex = matr.di.size() - 1;
	tmp[lIndex] = matr.al[lIndex - 1] * q[lIndex - 1] + matr.di[lIndex] * q[lIndex];
	return tmp;
}

double calcNormE(const vector1D &x) { return sqrt(x*x); }
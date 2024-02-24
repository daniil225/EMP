#include "FEM_1D.h"
#include "head.h"
/* Private */

bool FEM::shouldCalc(int i)
{
    // Выход по макисмальному количеству итераций 
    if(i > maxiter)
    {
        return false;
    }

	/* Выход по измененинию вектора решения */
	if(calcNormE(q-qPrev)/calcNormE(q) < 1e-4)
	{
		return false;
	}	

    // Вывод по невязке
	//cout << "Non-repan: " << calcNormE(MultAOnq(slau.Matr, q) - slau.f)/calcNormE(slau.f) << "\n";
    if( calcNormE(MultAOnq(slau.Matr, q) - slau.f)/calcNormE(slau.f) < eps )
    {	
        return false;
    }

    return true;
}

void FEM::GenerateProfile()
{
    int n = Grid.size();

    slau.Matr.ia[0] = 0;

    for(int i = 1; i < n+1; i++)
        slau.Matr.ia[i] = i-1;

}

void FEM::buildGlobalMatrixA(double _dt) 
{
	int nodesCount = Grid.size();
    int finiteElementsCount = nodesCount-1;
    dt = _dt;
	auto &di = slau.Matr.di;
    auto &au = slau.Matr.au;
    auto &al = slau.Matr.al;
	di.clear();
	au.clear();
	al.clear();
	
	di.resize(nodesCount, 0);
	al.resize(nodesCount - 1, 0);
	au.resize(nodesCount - 1, 0);
	
    for (size_t elemNumber = 0; elemNumber < finiteElementsCount; elemNumber++)
	{
		buildLocalmatrixA(elemNumber);

		//cout << ALocal << endl;
		di[elemNumber] += ALocal[0][0];		au[elemNumber] += ALocal[0][1];
		al[elemNumber] += ALocal[1][0];		di[elemNumber + 1] += ALocal[1][1];
	}

	// Первые краевые условия
	di[0] = 1;
	au[0] = 0;
	di[nodesCount - 1] = 1;
	al[al.size() - 1] = 0;
}

void FEM::buildGlobalVectorb() 
{
    auto &f = slau.f;
    int nodesCount = Grid.size();
    int finiteElementsCount = nodesCount-1;

    f.clear();
	f.resize(nodesCount, 0);

	for (size_t elemNumber = 0; elemNumber < finiteElementsCount; elemNumber++)
	{
		buildLocalVectorf(elemNumber);
		f[elemNumber] += bLocal[0];
		f[elemNumber + 1] += bLocal[1];
	}

	f[0] = u(Grid[0], t);
	f[nodesCount - 1] = u(Grid[nodesCount - 1], t);
}

void FEM::printGlobalMatrixA() 
{

}

void FEM::printGlobalVectorb() 
{

}

void FEM::buildLocalMatrixG(int elemNumber) 
{
    double hx = Grid[elemNumber+1] - Grid[elemNumber];
	lambda0 = lambda(q[elemNumber]);
	lambda1 = lambda(q[elemNumber + 1]);
	double numerator = (lambda0 + lambda1) /(2.0 * hx);
	GLocal[0][0] = GLocal[1][1] = numerator;
	GLocal[0][1] = GLocal[1][0] = -numerator;
}

void FEM::buildLocalMatrixM(int elemNumber) 
{
    double hx = Grid[elemNumber+1] - Grid[elemNumber];
	double numerator = (sigma * hx) / (6 * dt);
	MLocal[0][0] = MLocal[1][1] = 2 * numerator;
	MLocal[0][1] = MLocal[1][0] = numerator;
}

void FEM::buildLocalmatrixA(int elemNumber)
{
	ALocal = GLocal = MLocal = { {0,0}, {0,0} };
	buildLocalMatrixG(elemNumber);
	buildLocalMatrixM(elemNumber);
	for (size_t i = 0; i < 2; i++)
	{
		for (size_t j = 0; j < 2; j++)
		{
			ALocal[i][j] = GLocal[i][j] + MLocal[i][j];
		}
	}
}

void FEM::buildLocalVectorf(int elemNumber)
{
    double hx = Grid[elemNumber+1] - Grid[elemNumber];
	bLocal = { 0, 0 };
	bLocal[0] = hx * (2.0 * f(Grid[elemNumber], t) + f(Grid[elemNumber + 1], t)) / 6.0
		+ sigma *hx * (2.0 * qPrev[elemNumber] + qPrev[elemNumber + 1]) / (6.0 * dt);
	bLocal[1] = hx * (f(Grid[elemNumber], t) + 2.0 * f(Grid[elemNumber + 1], t)) / 6.0
		+ sigma * hx *(qPrev[elemNumber] + 2.0 * qPrev[elemNumber + 1]) / (6.0 * dt);
}

/* Public */
void FEM::init(const function2D &_u, const function2D &_f, const function1D &_lambda, double _sigma, const string &Grid_, const string &TimeGrid_) 
{
    Grid = Grid_1D(Grid_);
    TimeGrid = Grid_1D(TimeGrid_);
	Grid.GenerateGrid();
	TimeGrid.GenerateGrid();

    u = _u;
	f = _f;
	lambda = _lambda;
	sigma = _sigma;

    /* Allocate memory for matrix and right part */
    int n = Grid.size();
    slau.Matr.di.resize(n);
    slau.Matr.au.resize(n-1);
    slau.Matr.al.resize(n-1);
    slau.Matr.ia.resize(n+1);
    slau.f.resize(n);
    /* Генерация профиля матрицы она имеет 3-х диагональную структуру */
    GenerateProfile();

    /* Память под вектора решений */
    q.resize(n);
    qPrev.resize(n);
	Q.resize(TimeGrid.size());

    /* Память под локальные матрицы  */
    GLocal = vector(2, vector<double>(2));
    MLocal = vector(2, vector<double>(2));
    ALocal = vector(2, vector<double>(2));

}

pair<int, double> FEM::solve() 
{
// Задаём начальные условия
int n = Grid.size();
	//q.resize(n, 0);
	//qPrev.resize(n, 0);
	vector<double> qExact(n);
	for (size_t i = 0; i < n; i++)
		qExact[i] = u(Grid[i], TimeGrid[0]);
	qPrev = qExact;

	Q[0] = qPrev;
	//cout << "QTrue: [" << qPrev << "]\n";
	int count = 0;
	// Решаем в каждый момент временной сетки
	double sumNormQ = 0;
	for (size_t i = 1; i < TimeGrid.size(); i++)
	{
		dt = TimeGrid[i] - TimeGrid[i - 1];
		t = TimeGrid[i];
		bool sC = true;
		do {
			buildGlobalMatrixA(dt);
			buildGlobalVectorb();
			SolveSlau(slau, q);
			count++; // Итеарция прошла
			qPrev = q; 
			sC =  shouldCalc(count);
			cout << "time: " << t << " Norma: " << calcNormAtMainNodes(q) << "\n"; 
			//cout << "time: " << t << " q = [" << q << "]\n";
			
		} while (sC);

		Q[i] = qPrev; // Заносим в решение очередной временной слой
		sumNormQ += calcNormAtMainNodes(q);
	}

	sumNormQ /= double(TimeGrid.size() -1);

	return make_pair(count, sumNormQ);
}

double FEM::CalculateU(double x, double t)
{
	double res = 0;
	// Пусть мы берем последний временной слой пока что 
	vector<double>& TmpQ = Q[Q.size()-1];
	
	/* Определим отрезок в котором будем расчитвыать значение функции */
	int32_t NumElement = 0;
	for(; NumElement < Grid.size()-1; NumElement++)
	{
		if(Grid[NumElement] <= x && Grid[NumElement+1] >= x) break;
	}

	double xm = Grid[NumElement];
	double xm1 = Grid[NumElement+1];

	auto Psi1 = [&](double x) { return (xm1 - x)/(xm1-xm); };
	auto Psi2 = [&](double x) { return (x-xm)/(xm1-xm); };

	return TmpQ[NumElement]*Psi1(x) + TmpQ[NumElement+1]*Psi2(x);
}
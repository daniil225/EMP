#ifndef FEM_H_
#define FEM_H_


#include "Grid_1D.h"
#include "Slau.h"
#include <string>
#include <fstream>
#include <fstream>
#include <functional>
#include <cmath>

using namespace std;


typedef function<double (double, double)> function2D;
typedef function<double (double)> function1D;

class FEM 
{
    private:
    Grid_1D Grid;
    Grid_1D TimeGrid;
    Slau slau; // Она же глобальная матрицы и вектор правой части 
    vector<double> q, qPrev;

    vector<vector<double>> Q; // Общее решение по временным слоям 

    double lambda0, lambda1;
    double sigma;
    double t;
    double dt;
    const int maxiter = 10000; // Максимальное количество итераций 
    double eps = 1e-7;
    
    int StepCoef = 1; // Шаговый коэффициент по пространству 
    int TimeCoef = 1; // Шаговый коэффициент по времени 

    function2D f, u;
    function1D lambda;

    void GenerateProfile();
    void buildGlobalMatrixA(double _dt);
	void buildGlobalVectorb();
	void printGlobalMatrixA();
	void printGlobalVectorb();

	void buildLocalMatrixG(int elemNumber);
	void buildLocalMatrixM(int elemNumber);
	void buildLocalmatrixA(int elemNumber);
	void buildLocalVectorf(int elemNumber);
    bool shouldCalc(int i);

    double calcNormAtMainNodes(const vector<double> &x)
    {
		double tmp = 0;
		for (size_t i = 0; i < x.size(); i+=StepCoef)
			tmp += pow((x[i] - u(Grid[i], t)), 2);
		return sqrt(tmp) / Grid.size();
	}


    vector<vector<double>>  GLocal ,MLocal, ALocal;
	vector<double> bLocal;

    public:

    FEM() = default;

    void init(const function2D &_u, const function2D &_f, const function1D &_lambda, double _sigma, const string &Grid, const string &TimeGrid);
	pair<int, double> solve();
	inline int getNodesCount() { return Grid.size(); }

    double CalculateU(double x, double t);

};


#endif
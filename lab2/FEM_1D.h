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
    vector<double> qExact; // Вектор на прошлой итерации по нелинейности 
    vector<vector<double>> Q; // Общее решение по временным слоям 

    double lambda0, lambda1;
    double sigma;
    double t;
    double dt;
    const int maxiter = 1000; // Максимальное количество итераций 
    double eps = 1e-7;
    double delta = 1e-7;
    
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

    double calcNormAtMainNodes()
    {
        double res = 0;
		auto normsub = [&](double x, double t)
        {
            return pow(u(x,t) - CalculateU(x,t), 2.0);
        };

        /* Расчет интеграла. Берем шаг h = 0.002 и пройдемся по отрезку и вычислим интеграл */
        double h = 0.0002;
        double start = Grid[0];
        int32_t N = int32_t((Grid[Grid.size()-1] - Grid[0])/h);

        for(int32_t i = 0; i < N; i++)
        {
            double a = start + i*h;
            double b = a + h;
            double arg = (b+a)/2.0;
            //cout << "[ " << a << "; " << b << "]\n";
            //cout << "hx = " << h  << " arg = " << arg << " normsub() = " << normsub(arg, 1.0) << "\n"; 
            res += h*normsub(arg, 1.0);
        }

        return sqrt(res);
	}


    vector<vector<double>>  GLocal ,MLocal, ALocal;
	vector<double> bLocal;

    public:

    FEM() = default;

    void init(const function2D &_u, const function2D &_f, const function1D &_lambda, double _sigma, const string &Grid, const string &TimeGrid);
	pair<int, double> solve();
	inline int getNodesCount() { return Grid.size(); }
    void DivideGridAndPrepareInternalParametrs(const int32_t coef);

    double CalculateU(double x, double t);

};


#endif
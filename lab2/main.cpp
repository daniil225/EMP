#include "FEM_1D.h"

function1D calcFirstDerivative(const function1D& f) {
	return [f](double x) -> double {
		const double h = 0.00001;
		return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (12 * h);
	};
}


function2D calcRightPart(const function1D& lambda, const function2D& u, double sigma) {
	return [=](double x, double t) -> double {
		using namespace std::placeholders;
		auto duBydt = calcFirstDerivative(std::bind(u, x, _1));
		auto duBydx = calcFirstDerivative(std::bind(u, _1, t));
		auto lambda_grad = [=](double x, double t) -> double {
			return lambda(u(x, t)) * duBydx(x);
			//return lambda(duBydt(t)) * duBydx(x);
		};
		auto div = calcFirstDerivative(std::bind(lambda_grad, _1, t));
		return -div(x) + sigma * duBydt(t);
	};
}


int main()
{

    TestSlau();
    vector <function2D> u(14), f(14);
	u[0] = { [](double x, double t) -> double { return 3 * x + t; } };
	u[1] = { [](double x, double t) -> double { return 2 * x*x; } };
	u[2] = { [](double x, double t) -> double { return x * x*x + t; } };
	u[3] = { [](double x, double t) -> double { return x * x*x*x + t; } };
	u[4] = { [](double x, double t) -> double { return exp(x) + t; } };
	u[5] = { [](double x, double t) -> double { return 3 * x + t; } };
	u[6] = { [](double x, double t) -> double { return 3 * x + t * t; } };
	u[7] = { [](double x, double t) -> double { return 3 * x + t * t*t; } };
	u[8] = { [](double x, double t) -> double { return 3 * x + exp(t); } };
	u[9] = { [](double x, double t) -> double { return 3 * x + sin(t); } };
	u[10] = { [](double x, double t) -> double { return exp(x) + t * t; } };
	u[11] = { [](double x, double t) -> double { return exp(x) + t * t*t; } };
	u[12] = { [](double x, double t) -> double { return exp(x) + exp(t); } };
	u[13] = { [](double x, double t) -> double { return exp(x) + sin(t); } };

	vector <function1D>  lambda(8);
	lambda[0] = { [](double u) -> double {return 1; } };
	lambda[1] = { [](double u) -> double {return u; } };
	lambda[2] = { [](double u) -> double {return u * u; } };
	lambda[3] = { [](double u) -> double {return u * u + 1; } };
	lambda[4] = { [](double u) -> double {return u * u*u; } };
	lambda[5] = { [](double u) -> double {return u * u*u*u; } };
	lambda[6] = { [](double u) -> double {return exp(u); } };
	lambda[7] = { [](double u) -> double {return sin(u); } };

    double sigma = 1;

    f[0] = calcRightPart(lambda[1], u[1], sigma);
	//f[0] = [](double x, double t) -> double { return 1-exp(x); };
	//cout << f[0](0.3, 0.5) << "\n";

    FEM fem;
    fem.init(u[1], f[0], lambda[1], sigma, "Grid.txt", "TimeGrid.txt");
    auto res = fem.solve();

	cout << "UCalc(0.3, 1) = " << fem.CalculateU(0.3, 1.0) << "\n";
	cout << "|UCalc(0.3, 1) - u(0.3, 1)| = " << abs(fem.CalculateU(0.3, 1.0) - u[1](0.3, 1.0)) << "\n";
    cout << res.first << " " << res.second << "\n";
    return 0;
}
#include "FEM_1D.h"
#include <string>
#include <unistd.h> 

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

    vector <function2D> u(14), f(14);
	u[0] = { [](double x, double t) -> double { return 3 * x + t; } };
	u[1] = { [](double x, double t) -> double { return 2 * x*x + t; } };
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
	lambda[7] = { [](double u) -> double {return 2+sin(u); } };

    double sigma = 1;

	std::vector<std::string> lambda_str = {"1", "u", "u*u", "u*u + 1", "u^3", "u^4", "exp(u)", "2 + sin(u)"};
	std::vector<std::string> u_str = {"3x + t", "2x^2 + t", "x^3 + t","x^4 + t", "exp(x) + t", "3x + t", "3x + t^2", "3x + t^3", "3x + exp(t)", "3x + sin(t)", "exp(x) + t^2", "exp(x) + t^3", "exp(x) + exp(t)", "exp(t) + sin(t)"};

	/* Генерация таблички для различных lambda */
	
	for(int i = 0; i < lambda.size(); i++)
	{
		printf("lambda(u) = %s\n", lambda_str[i].c_str());
		printf("|-------------------------------------------------------------------------------------------------|\n");
		printf("|         u(x,t)          |          || u(x,t)* - u(x,t) ||             |        Iteration        |\n");
		printf("|-------------------------------------------------------------------------------------------------|\n");
		for(int j = 0; j < u.size(); j++)
		{
			f[j] = calcRightPart(lambda[i], u[j], sigma);
			FEM fem;
    		fem.init(u[j], f[j], lambda[i], sigma, "Grid.txt", "TimeGrid.txt");
			auto res = fem.solve();
			printf("|%s|%.6f||%d|\n", u_str[j].c_str(), res.second, res.first);
			printf("|-------------------------------------------------------------------------------------------------|\n");
			sleep(1);
		}
		printf("\n\n\n");
	}
	


    f[0] = calcRightPart(lambda[0], u[2], sigma);
	//f[0] = [](double x, double t) -> double { return 1; };

   FEM fem;
   fem.init(u[2], f[0], lambda[0], sigma, "Grid.txt", "TimeGrid.txt");
   auto res = fem.solve();

	//cout << "UCalc(0.3, 1) = " << fem.CalculateU(0.3, 1.0) << "\n";
	//cout << "|UCalc(0.3, 1) - u(0.3, 1)| = " << abs(fem.CalculateU(0.3, 1.0) - u[1](0.3, 1.0)) << "\n";
    cout << "Iterations = " <<  res.first << "\n  Norm = " << res.second << "\n";
    return 0;
}
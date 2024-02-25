#include "FEM_1D.h"
#include <string>
#include <unistd.h> 


function1D calcFirstDerivative(const function1D& f) {
	return [f](double x) -> double {
		const double h = 0.00001;
		return (f(x - 2 * h) - 8 * f(x - h) + 8 * f(x + h) - f(x + 2 * h) ) / (12 * h);
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

/* Печатет табличку для исследования на сходимость */
void PrintTable(vector<string> &param_str, int CountOfDivide, function1D &lambda_, function2D &u_)
{
	/* 0 = lambda_str
	   1 = u_str
	*/
	FEM fem;
	function2D f_ = calcRightPart(lambda_, u_, 1);
    fem.init(u_, f_, lambda_, 1, "Grid.txt", "TimeGrid.txt");

	pair<int, double> res_Prev = fem.solve();
	pair<int, double> res_Curr;
	printf("lambda(u) = %s  u(x,t) = %s", param_str[0].c_str(), param_str[1].c_str());
	printf("|--------------------------------------------------------------------------------------------------|\n");
	printf("|  N  |  Nodes  |   Iteration   |  || u(x,t)* - u(x,t) ||  |  log(|| * ||(h)/|| * ||(h/2))  |\n");
	printf("|--------------------------------------------------------------------------------------------------|\n");
	printf("|%3d  |%6d   |%12d   |%22e    |%30e    |\n", 1, fem.getNodesCount(), res_Prev.first, 0); // Первая строка в таблице

	for(int i = 2; i <= CountOfDivide; i++)
	{
		fem.DivideGridAndPrepareInternalParametrs(2);
		res_Curr = fem.solve();
		printf("|%3d  |%6d   |%12d   |%22e    |%30e    |\n", i, fem.getNodesCount(), res_Curr.first, log2(res_Prev.second/res_Curr.second));
		printf("|--------------------------------------------------------------------------------------------------|\n");
		res_Prev = res_Curr;
	}
	printf("\n\n\n");
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
	

	// for(int i = 0; i < lambda.size(); i++)
	// {
	// 	printf("lambda(u) = %s\n", lambda_str[i].c_str());
	// 	printf("|------------------------------------------------------------------|\n");
	// 	printf("|         u(x,t)          |  || u(x,t)* - u(x,t) ||  |  Iteration  |\n");
	// 	printf("|------------------------------------------------------------------|\n");
	// 	for(int j = 0; j < u.size(); j++)
	// 	{
	// 		f[j] = calcRightPart(lambda[i], u[j], sigma);
	// 		FEM fem;
    // 		fem.init(u[j], f[j], lambda[i], sigma, "Grid.txt", "TimeGrid.txt");
	// 		auto res = fem.solve();
	// 		printf("|%23s  |%20e      |%8d     |\n", u_str[j].c_str(), res.second, res.first);
	// 		printf("|------------------------------------------------------------------|\n");
	// 		//sleep(1);
	// 	}
	// 	printf("\n\n\n");
	// }

	FEM fem;
	f[3] = calcRightPart(lambda[3], u[3], sigma);
	fem.init(u[3], f[3], lambda[3], sigma, "Grid.txt", "TimeGrid.txt");
	auto res_Prev = fem.solve();

	printf("lambda(u) = %s  u(x,t) = %s\n", lambda_str[3].c_str(), u_str[3].c_str());
	printf("|--------------------------------------------------------------------------------------------------|\n");
	printf("|  N  |  Nodes  |   Iteration   |  || u(x,t)* - u(x,t) ||  |  log(|| * ||(h)/|| * ||(h/2))  |\n");
	printf("|--------------------------------------------------------------------------------------------------|\n");
	printf("|%3d  |%6d   |%12d   |%22e    |%30e    |\n", 1, fem.getNodesCount(), res_Prev.first, res_Prev.second, 0); // Первая строка в таблице


    return 0;
}
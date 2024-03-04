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
	lambda[7] = { [](double u) -> double {return 2*u+sin(15*u); } };

    double sigma = 1;

	std::vector<std::string> lambda_str = {"1", "u", "u^2", "u^2 + 1", "u^3", "u^4", "e^u", "2u + sin(15u)"};
	std::vector<std::string> u_str = {"3x + t", "2x^2 + t", "x^3 + t","x^4 + t", "e^x + t", "3x + t", "3x + t^2", "3x + t^3", "3x + e^t", "3x + sin(t)", "e^x + t^2", "e^x + t^3", "e^x + e^t", "e^t + sin(t)"};

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

	int Lambda_num = 1;
	/* Делаем проход по всем функциям в списке и фиксированной lambda */
	// for(int k = 0; k < u.size(); k++)
	// {

	// 	FEM fem;
	// 	f[k] = calcRightPart(lambda[Lambda_num], u[k], sigma);
	// 	fem.init(u[k], f[k], lambda[Lambda_num], sigma, "Grid.txt", "TimeGrid.txt");
	// 	auto res_Prev = fem.NutonSolve(); //  Решение уравнения по методу Ньютона 
	// 	pair<int, double> res_Curr;
	// 	printf("$\\lambda(u) = %s \\space\\space\\space  u(x,t) = %s$\n", lambda_str[Lambda_num].c_str(), u_str[k].c_str());
	// 	printf("\n");
	// 	//printf("|-------------------------------------------------------------------------------------------|\n");
	// 	printf("|  N  |  Nodes  |   Iteration   |$\\|u(x,t)^*-u(x,t)\\|_ {L_2}$|$\\log(\\frac{ \\| u(x,t)^* - u(x,t)_ {h} \\|_ {L_2} }{\\| u(x,t)^* - u(x,t)_ {\\frac{h}{2}} \\|_ {L_ 2}})$|\n");
	// 	printf("|-----|:-------:|:-------------:|:-----------------------------:|:-------------------------:|\n");
	// 	printf("|%2d   |%6d   |%7d        |%18e        |%16s                |\n", 1, fem.getNodesCount(), res_Prev.first, res_Prev.second, "0"); // Первая строка в таблице
	// 	//printf("\n");
	// 	//printf("|-------------------------------------------------------------------------------------------|\n");

	// 	for(int i = 2; i <= 5; i++)
	// 	{
	// 		fem.DivideGridAndPrepareInternalParametrs(2);
	// 		res_Curr = fem.NutonSolve();
	// 		printf("|%2d   |%6d   |%7d        |%18e        |        %13e           |\n", i, fem.getNodesCount(), res_Curr.first, res_Curr.second, log2(res_Prev.second/res_Curr.second) ); // Первая строка в таблице
	// 		//printf("|-------------------------------------------------------------------------------------------|\n");
	// 		res_Prev = res_Curr;
	// 	}

	// 	printf("\n\n\n");
	// }

		for(int k = 0; k < u.size(); k++)
	{

		FEM fem;
		f[k] = calcRightPart(lambda[Lambda_num], u[k], sigma);
		fem.init(u[k], f[k], lambda[Lambda_num], sigma, "Grid.txt", "TimeGrid.txt");
		auto res_Prev_NI = fem.NutonSolve(); //  Решение уравнения по методу Ньютона 
		auto res_Prev_SI = fem.solve();
		pair<int, double> res_Curr_SI;
		pair<int, double> res_Curr_NI;
		printf("$\\lambda(u) = %s \\space\\space\\space  u(x,t) = %s$\n", lambda_str[Lambda_num].c_str(), u_str[k].c_str());
		printf("\n");
		double r = (double)res_Prev_SI.first/(double)res_Prev_NI.first;
		//printf("|-------------------------------------------------------------------------------------------|\n");
		printf("|  N  |  Nodes  | Iteration SI | Iteration NI | SI - NI | Ускорение |\n");
		printf("|-----|:-------:|:------------:|:------------:|:-------:|:---------:|\n");
		printf("|%5d|%9d|%14d|%14d|%9d|%11.3lf|\n", 1, fem.getNodesCount(), res_Prev_SI.first, res_Prev_NI.first,
		res_Prev_SI.first - res_Prev_NI.first , r); // Первая строка в таблице
		//printf("\n");
		//printf("|-------------------------------------------------------------------------------------------|\n");

		for(int i = 2; i <= 5; i++)
		{
			fem.DivideGridAndPrepareInternalParametrs(2);
			res_Prev_NI = fem.NutonSolve();
			res_Prev_SI = fem.solve();
			r = (double)res_Prev_SI.first/(double)res_Prev_NI.first;
			printf("|%5d|%9d|%14d|%14d|%9d|%11.3lf|\n", i, fem.getNodesCount(), res_Prev_SI.first, res_Prev_NI.first,
		res_Prev_SI.first - res_Prev_NI.first , r); // Первая строка в таблице
			//printf("|-------------------------------------------------------------------------------------------|\n");
			
		}

		printf("\n\n\n");
	}

    return 0;
}
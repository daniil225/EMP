#ifndef TEST_H_
#define TEST_H_
#include "FEMSolver.h"
#include <cmath>

std::function<double(double, double, double)> div_grad(std::function<double(double, double, double)>& u);


ParamDE Test1();
ParamDE Test2();
ParamDE Test3();
ParamDE Test4();
ParamDE Test5();
ParamDE Test6();
ParamDE Test7();
ParamDE Test8();
ParamDE Test9();
ParamDE Test10();



#endif




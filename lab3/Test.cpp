#include "Test.h"
#include <cmath>

ParamDE Test1()
{
    ParamDE Test;
    Test.ksi = -2;
    Test.lambda = 2;
    Test.omega = 2;
    Test.sigma = 0;

    Test.f = [](double x, double y, double z, double t) { return -2.0*(x+y+z)*sin(t); };
    Test.u = [](double x, double y, double z, double t) { return (x+y+z)*(sin(t) + cos(t)); };

    Test.f_c = [](double x, double y, double z) { return /*2*(x+y+z)*/  8*(x+2*y+z); };
    Test.f_s = [](double x, double y, double z) { return /*-2*(x+2*y+z)*/  8*(x+y+z); };

    Test.u_c = [](double x, double y, double z) { return x+2*y+z; };
    Test.u_s = [](double x, double y, double z) { return x+y+z; };

    return Test;
}
#include "Test.h"
#include <cmath>


std::function<double(double, double, double)> div_grad(std::function<double(double, double, double)>& u)
{
    std::function<double(double, double, double)> divgrad_u = [&](double x, double y, double z) -> double 
    {
        double h = 0.0001;
        double d2udx2 = (u(x-h, y, z) - 2*u(x,y,z) + u(x+h, y, z))/(h*h);
        double d2udy2 = (u(x, y-h, z) - 2*u(x,y,z) + u(x, y+h, z))/(h*h);
        double d2udz2 = (u(x, y, z-h) - 2*u(x,y,z) + u(x, y, z+h))/(h*h);

        return d2udx2 + d2udy2 + d2udz2;
    };

    return divgrad_u;
}

/* Валидация  */
ParamDE Test1()
{
    ParamDE Test;
    Test.ksi = 0;
    Test.lambda = 2;
    Test.omega = 3;
    Test.sigma = 1;

    Test.u_c = [](double x, double y, double z) { return x*x*x+y*y*y+z*z*z; };
    Test.u_s = [](double x, double y, double z) { return  7*x*x*x+y*y*y+4*z*z*z; };


    Test.f_c = [&](double x, double y, double z) { return -Test.lambda*div_grad(Test.u_c)(x,y,z) + Test.omega*Test.sigma*Test.u_s(x,y,z) - Test.omega*Test.omega*Test.ksi*Test.u_c(x,y,z); };
    Test.f_s = [&](double x, double y, double z) { return -Test.lambda*div_grad(Test.u_s)(x,y,z) -Test.omega*Test.sigma*Test.u_c(x,y,z) - Test.omega*Test.omega*Test.ksi*Test.u_s(x,y,z); };

    /* Генерация 2-ых КУ */
    Test.Theta_s.resize(6);
    Test.Theta_c.resize(6);

    Test.Theta_c[0] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_c(x, y+h, z) - Test.u_c(x, y-h, z))/(2.0*h); };
    Test.Theta_c[1] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_c(x, y+h, z) - Test.u_c(x, y-h, z))/(2.0*h); };
    Test.Theta_c[2] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_c(x+h, y, z) - Test.u_c(x-h, y, z))/(2.0*h); };
    Test.Theta_c[3] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_c(x+h, y, z) - Test.u_c(x-h, y, z))/(2.0*h); };
    Test.Theta_c[4] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_c(x, y, z+h) - Test.u_c(x, y, z-h))/(2.0*h); };
    Test.Theta_c[5] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_c(x, y, z+h) - Test.u_c(x, y, z-h))/(2.0*h); };

    Test.Theta_s[0] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_s(x, y+h, z) - Test.u_s(x, y-h, z))/(2.0*h); };
    Test.Theta_s[1] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_s(x, y+h, z) - Test.u_s(x, y-h, z))/(2.0*h); };
    Test.Theta_s[2] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_s(x+h, y, z) - Test.u_s(x-h, y, z))/(2.0*h); };
    Test.Theta_s[3] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_s(x+h, y, z) - Test.u_s(x-h, y, z))/(2.0*h); };
    Test.Theta_s[4] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_s(x, y, z+h) - Test.u_s(x, y, z-h))/(2.0*h); };
    Test.Theta_s[5] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_s(x, y, z+h) - Test.u_s(x, y, z-h))/(2.0*h); };
    
    
    return Test;
}


ParamDE Test2()
{
    ParamDE Test;
    Test.ksi = 10e-11;
    Test.lambda = 10e3;
    Test.omega = 10e5;
    Test.sigma = 10e4;

    Test.u_c = [](double x, double y, double z) { return std::cos(x+y)+ 4*z*z*std::tan(z); };
    Test.u_s = [](double x, double y, double z) { return std::sin(x+y) - 12*z*z*std::log2(5+z); };


    Test.f_c = [&](double x, double y, double z) { return -Test.lambda*div_grad(Test.u_c)(x,y,z) + Test.omega*Test.sigma*Test.u_s(x,y,z) - Test.omega*Test.omega*Test.ksi*Test.u_c(x,y,z); };
    Test.f_s = [&](double x, double y, double z) { return -Test.lambda*div_grad(Test.u_s)(x,y,z) -Test.omega*Test.sigma*Test.u_c(x,y,z) - Test.omega*Test.omega*Test.ksi*Test.u_s(x,y,z); };

    /* Генерация 2-ых КУ */
    Test.Theta_s.resize(6);
    Test.Theta_c.resize(6);

    Test.Theta_c[0] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_c(x, y+h, z) - Test.u_c(x, y-h, z))/(2.0*h); };
    Test.Theta_c[1] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_c(x, y+h, z) - Test.u_c(x, y-h, z))/(2.0*h); };
    Test.Theta_c[2] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_c(x+h, y, z) - Test.u_c(x-h, y, z))/(2.0*h); };
    Test.Theta_c[3] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_c(x+h, y, z) - Test.u_c(x-h, y, z))/(2.0*h); };
    Test.Theta_c[4] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_c(x, y, z+h) - Test.u_c(x, y, z-h))/(2.0*h); };
    Test.Theta_c[5] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_c(x, y, z+h) - Test.u_c(x, y, z-h))/(2.0*h); };

    Test.Theta_s[0] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_s(x, y+h, z) - Test.u_s(x, y-h, z))/(2.0*h); };
    Test.Theta_s[1] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_s(x, y+h, z) - Test.u_s(x, y-h, z))/(2.0*h); };
    Test.Theta_s[2] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_s(x+h, y, z) - Test.u_s(x-h, y, z))/(2.0*h); };
    Test.Theta_s[3] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_s(x+h, y, z) - Test.u_s(x-h, y, z))/(2.0*h); };
    Test.Theta_s[4] = [&](double x, double y, double z) -> double{ double h = 0.00001; return -Test.lambda*(Test.u_s(x, y, z+h) - Test.u_s(x, y, z-h))/(2.0*h); };
    Test.Theta_s[5] = [&](double x, double y, double z) -> double{ double h = 0.00001; return  Test.lambda*(Test.u_s(x, y, z+h) - Test.u_s(x, y, z-h))/(2.0*h); };
    
    return Test;
    
}
ParamDE Test3()
{

}
ParamDE Test4()
{

}
ParamDE Test5()
{

}
ParamDE Test6()
{

}
ParamDE Test7()
{

}
ParamDE Test8()
{

}
ParamDE Test9()
{

}

ParamDE Test10()
{

}
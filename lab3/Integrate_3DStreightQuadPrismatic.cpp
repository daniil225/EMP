#include "Integrate_3DStreightQuadPrismatic.h"

#include <algorithm>

/* private */
std::pair<double, double> Integrate_3DStreightQuadPrismatic::Method_Rectangle(double eps, int32_t MaxInter)
{
    std::pair<double, double> res = std::make_pair<double, double>(-1, -1);
    double res_h = 0;
    double res_h2 = 0;
    int32_t Iter = 1;
    double curreps = 1;
    /* Строим сетку */
    GridStatus Status = Grid.SetBaseGrid(baseGrid);
    /* Проверяем состояние */
    if (Status.GetState() != State::OK)
    {
        std::cerr << "Error: " << Status.GetMsg();
        return res;
    }

    Status = Grid.GenerateGrid(); // Генерация сетки
    /* Проверяем состояние */
    if (Status.GetState() != State::OK)
    {
        std::cerr << "Error: " << Status.GetMsg();
        return res;
    }

    /* Сам алгоритм обхода  */

    while (Iter <= MaxInter && eps < curreps)
    {
        res_h = res_h2 = 0;

        /* Решение на сетке с шагом h */
        for (int32_t i = 0; i < Grid.GetElementCount(); i++)
        {
            StreightQuadPrismatic Element = Grid.GetElement(i);

            /* Пока что предполагаем, что это произвольный параллелепипед основание в плоскости XZ, нумерация узлов с нуля */
            double x0 = Element.points[0].x;
            double y0 = Element.points[0].y;
            double z0 = Element.points[0].z;

            double x1 = Element.points[1].x;
            double y1 = Element.points[4].y;
            double z1 = Element.points[2].z;

            double x = x0 + (x1 - x0) / 2.0;
            double y = y0 + (y1 - y0) / 2.0;
            double z = z0 + (z1 - z0) / 2.0;
            double J = (x1 - x0) * (y1 - y0) * (z1 - z0);
            res_h += J * f(x, y, z);

        }

        /* Решение на сетке с шагом h/2 */
        Status = Grid.DivideGrid(2);
        /* Проверяем состояние */
        if (Status.GetState() != State::OK)
        {
            std::cerr << "Error: " << Status.GetMsg();
            return res;
        }

        Status = Grid.ReGenerateGrid();

        if (Status.GetState() != State::OK)
        {
            std::cerr << "Error: " << Status.GetMsg();
            return res;
        }

        for (int32_t i = 0; i < Grid.GetElementCount(); i++)
        {
            StreightQuadPrismatic Element = Grid.GetElement(i);

            /* Пока что предполагаем, что это произвольный параллелепипед основание в плоскости XZ, нумерация узлов с нуля */
            double x0 = Element.points[0].x;
            double y0 = Element.points[0].y;
            double z0 = Element.points[0].z;

            double x1 = Element.points[1].x;
            double y1 = Element.points[4].y;
            double z1 = Element.points[2].z;

            double x = x0 + (x1 - x0) / 2.0;
            double y = y0 + (y1 - y0) / 2.0;
            double z = z0 + (z1 - z0) / 2.0;
            double J = (x1 - x0) * (y1 - y0) * (z1 - z0);
            res_h2 += J * f(x, y, z);
        }
        Iter += 1;
        curreps= std::abs(res_h - res_h2) / 3.0;
    }

    res.first = res_h2;
    res.second = curreps;

    return res;
}

/* public */

std::pair<double, double> Integrate_3DStreightQuadPrismatic::Quad(std::string Method, double eps, int32_t MaxIter)
{
    std::pair<double, double> res = std::make_pair<double, double>(-1, -1);

    if (std::find(std::begin(Methods), std::end(Methods), Method) != std::end(Methods))
    {
        if (Method == Methods[0] || Method == Methods[1])
        {
            res = Method_Rectangle(eps, MaxIter);
        }
    }
    else
    {
        throw "Ты что дебил?? \n";
    }

    return res;
}

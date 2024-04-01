#ifndef INTEGRATE_3DSTREIGHT_QUAD_PRISMATIC_
#define INTEGRATE_3DSTREIGHT_QUAD_PRISMATIC_

#include "Grid3D_Integrate_StreightQuadPrismatic.h"
#include <functional>
#include <tuple>
#include <string>

class Integrate_3DStreightQuadPrismatic
{
    private:
        BaseGrid3D_Integrate_StreightQuadPrismatic baseGrid;
        Grid3D_Integrate_StreightQuadPrismatic Grid;
        std::function<double(double, double, double)> f; // Интегрируемая функция 
        const std::string Methods[2] = {"default", "Rectangle"};

        /*
            @brief  Интегрирование по методу прямоугольников. Пока что реализованное для произвольного параллелипипеда. 
                    В будущем можно будет реализовать интегрирование по произвольному призматическому элементу.
            @param double eps = 1e-7 - Желаемая точность расчета  
            @param int32_t MaxIter = 10 максимальное число дроблений расчетной области. Дробление происходит на 2.  
            @return std::pair<double, double> - пара значений. 0 - значение интеграла. 1 - точность полученного решения. 
        */
        std::pair<double, double> Method_Rectangle(double eps = 1e-7, int32_t MaxInter = 10);

    protected:

    public:

    /*
        @brief Конструктор по умолчанияю 
    */
    Integrate_3DStreightQuadPrismatic() = default;



    /*
        @brief Конструктор инициализации
        @param BaseGrid3D_Integrate_StreightQuadPrismatic &baseGrid_ - 
        @param std::function<double(double, double, double)>& f_ - 
        @result Установится сетка и функция для расчета
    */
    Integrate_3DStreightQuadPrismatic(const BaseGrid3D_Integrate_StreightQuadPrismatic &baseGrid_, const std::function<double(double, double, double)>& f_):
        baseGrid(baseGrid_), f(f_) {}


    /*
        @brief Интегрирование по заданной области.  
        @param std::string Method = "default" - Метод расчета интеграла. По умолчанию это метод прямоугольников. Потом может быть добавлю методы Гауса
        @param double eps = 1e-7 - Желаемая точность расчета  
        @param int32_t MaxIter = 10 максимальное число дроблений расчетной области. Дробление происходит на 2.  
        @return std::pair<double, double> - пара значений. 0 - значение интеграла. 1 - точность полученного решения.   
    */
   std::pair<double, double> Quad(std::string Method = "default", double eps = 1e-7, int32_t MaxInter = 10);


    ~Integrate_3DStreightQuadPrismatic() = default;

};

#endif
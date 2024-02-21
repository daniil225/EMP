#ifndef GRID_1D_H
#define GRID_1D_H
#include <vector>
#include <string>
#include <fstream>
using  namespace std;


/* Формат файла */
/* 
    * Nx - Количество опорных узлов по оси x
    * x1 x2 ... xn
    * далее коэффиценты дробления сетки 
    * n1 q1 n2 q2 ... 

 */

class Grid_1D
{
    private:
    const double eps = 1e-10;
    ifstream fin;
    ofstream fout;
    struct DivideParamS
    {
        int num;     // количество интервалов на которое нужно разделить отрезок
        double coef; // Коэффициент растяжения или сжатия
    };
    int32_t Nx;
    vector<double> BaseGridX;
    vector<DivideParamS> DivideParam;
    vector<double> Grid;
    int stepCoef = 1; // Коэффициент с шагом для базовой сетки умножается при дроблениях 

    public:

    Grid_1D() = default;
    Grid_1D(const string &filename);
    void Load(const string &filenmae);

    void GenerateGrid();
    void DivideGrid(const int coef);
    void ReGenerateGrid();

    inline int32_t size() const { return (int32_t)Grid.size(); }
    inline double& operator[](const int32_t idx) { return Grid[(uint64_t)idx]; }


    void PrintGrid() const;

    Grid_1D& operator=(const Grid_1D &Grid_)
    {
        this->BaseGridX = Grid_.BaseGridX;
        this->Nx = Grid_.Nx;
        this->DivideParam = Grid_.DivideParam;
        this->Grid = Grid_.Grid;

        return *this;
    }

    ~Grid_1D() = default;
};



#endif
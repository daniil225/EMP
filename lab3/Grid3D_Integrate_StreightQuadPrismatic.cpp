#include "Grid3D_Integrate_StreightQuadPrismatic.h"
#include <functional>
#include <cmath>
#include <iostream>

using namespace std;

/* Private */

void Grid3D_Integrate_StreightQuadPrismatic::GetTotalNumberOfNodes() noexcept
{
    for (int32_t i = 0; i < baseGrid.Nx - 1; i++)
        GlobalNx += baseGrid.DivideParam[0][static_cast<uint64_t>(i)].num;

    for (int32_t i = 0; i < baseGrid.Nz - 1; i++)
        GlobalNz += baseGrid.DivideParam[1][static_cast<uint64_t>(i)].num;

    for (int32_t i = 0; i < baseGrid.Ny - 1; i++)
        GlobalNy += baseGrid.DivideParam[2][static_cast<uint64_t>(i)].num;

    ElementCount = GlobalNx * GlobalNy * GlobalNz;
    GlobalNx++;
    GlobalNy++;
    GlobalNz++;
    Dim = GlobalNx * GlobalNy * GlobalNz;
}

int32_t Grid3D_Integrate_StreightQuadPrismatic::Getlevel(int32_t i, int32_t axis) const noexcept
{
    int32_t res = 0;
    for (int32_t k = 0; k < i; k++)
        res += baseGrid.DivideParam[static_cast<uint64_t>(axis)][static_cast<uint64_t>(k)].num;
    return res;
}

void Grid3D_Integrate_StreightQuadPrismatic::GenerateBaseGrid() noexcept
{
    /* Вспомогальельная структура для определения параметров разбиения */

    struct SettingForDivide
    {
        double step = 0; // Шаг на отрезке
        double coef = 0; // Коэффициент увеличения шага
        int32_t num = 0; // Количество интервалов идем то num-1 и потом явно вставляем элемент

        /* Копирование и присваивание */
        SettingForDivide() = default;
        SettingForDivide(const SettingForDivide &) = default;
        SettingForDivide &operator=(const SettingForDivide &) = default;
        SettingForDivide(SettingForDivide &&) = default;
        SettingForDivide &operator=(SettingForDivide &&) = default;
    };
    /*******************************************************************/

    /* Вспомогательные функции для генерации сетки опеределены в виде лямбда функций */

    /*
        @param:
            int32_t i - Номер массива от 0 до 2
            int32_t j - Номер элемента в массиве
            double left - левая грани отрезка
            double right - правая граница отрезка
        @return: SettingForDivide -  структура с вычесленными параметрами деления сетки
        @result: Расчитываем шаг для сетки
    */
    std::function<SettingForDivide(int32_t, int32_t, double, double)> CalcSettingForDivide =
        [&](int32_t i, int32_t j, double left, double right) -> SettingForDivide
    {
        SettingForDivide res;
        int32_t num = baseGrid.DivideParam[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].num;
        double coef = baseGrid.DivideParam[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].coef;

        if (std::abs(coef - 1.0) > eps)
        {
            double coefStep = 1.0 + (coef * (std::pow(coef, num - 1) - 1.0)) / (coef - 1.0);

            res.step = (right - left) / coefStep;
        }
        else
        {
            res.step = (right - left) / num;
        }

        //  Убираем погрешность
        if (std::abs(res.step) < eps)
            res.step = 0.0;

        res.num = num;
        res.coef = coef;
        return res;
    };

    /*
        @param:
            SettingForDivide &param - параметр разбиения
            double left - левая граница отрезка
            double right - правая граница отрезка
            vector<double>& Line - генерируемый массив
            int &idx - индекс в массиве на какую позицию ставить элемент
        @return: void
        @result: Генерация разбиения по X или Z( гороизонтальная линия или вертикальная ) с учетом разбиения
    */
    std::function<void(SettingForDivide &, double, double, vector<double> &, int &idx)> GenerateDivide =
        [](SettingForDivide &param, double left, double right, vector<double> &Line, int &idx) -> void
    {
        int32_t num = param.num;
        double coef = param.coef;
        double step = param.step;

        Line[static_cast<uint64_t>(idx)] = left;
        idx++;
        double ak = left;
        for (int32_t k = 0; k < num - 1; k++)
        {
            ak = ak + step * std::pow(coef, k);
            Line[static_cast<uint64_t>(idx)] = ak;
            idx++;
        }
        Line[static_cast<uint64_t>(idx)] = right;
    };

    try
    {
        /* Подготовка переменных для генерации сетки */
        int32_t total = GlobalNx * GlobalNy * GlobalNz;
        Grid.resize(static_cast<uint64_t>(total));

        /* Псевдонимы поддержка локализации */
        int32_t Nx = baseGrid.Nx;
        int32_t Ny = baseGrid.Ny;
        int32_t Nz = baseGrid.Nz;
        vector<vector<BaseGrid3D_Integrate_StreightQuadPrismatic::PointXZS>> &BaseGridXZ = baseGrid.BaseGridXZ;
        vector<double> &BaseGridY = baseGrid.BaseGridY;

        /* Глобальные размеры есть */

        /* разбиения покоординатные */
        vector<double> LineX(static_cast<uint64_t>(GlobalNx)); // Массив элементов в строке по Х
        vector<double> LineZ(static_cast<uint64_t>(GlobalNz)); // Массив элементов в строке по Z
        vector<double> LineY(static_cast<uint64_t>(GlobalNy)); // Массив элементов в строке по Y
                                                               /**********************************************************************/

        /* Блок генерации разбиения */

        /*  Описание основной идеи генерации разбиения
            1) Сгенерируем одну плоскость xz, а потом ее растиражируем по y разбиение по x и z
            2) Расстановка элементов основных линий с учетом их расположения в сетке ( Опорная сетка ):
                Расстановка элементов по x ( Опорных )
                Нужно значения х расставить в соответсвующие строки они соответсвуют разбиению по z
        */

        for (int32_t i = 0; i < Nz; i++)
        {
            int32_t idx = 0;

            for (int32_t j = 0; j < Nx - 1; j++)
            {
                double left = BaseGridXZ[static_cast<uint64_t>(i)][static_cast<uint64_t>(j)].x;
                double right = BaseGridXZ[static_cast<uint64_t>(i)][static_cast<uint64_t>(j + 1)].x;
                SettingForDivide param = CalcSettingForDivide(0, j, left, right);
                GenerateDivide(param, left, right, LineX, idx);
            }
            /* Заносим соответствующие значения x на свои позиции */

            int32_t startIdx = Getlevel(i, 1) * GlobalNx;
            int32_t endIdx = startIdx + GlobalNx;
            for (int32_t k = startIdx, kk = 0; k < endIdx; k++, kk++)
                Grid[static_cast<uint64_t>(k)].x = LineX[static_cast<uint64_t>(kk)];
        }

        /* Расстановка элементов по z ( Опорных ) */
        for (int32_t i = 0; i < Nx; i++)
        {
            int32_t idx = 0;
            for (int32_t j = 0; j < Nz - 1; j++)
            {
                double left = BaseGridXZ[static_cast<uint64_t>(j)][static_cast<uint64_t>(i)].z;
                double right = BaseGridXZ[static_cast<uint64_t>(j + 1)][static_cast<uint64_t>(i)].z;
                SettingForDivide param = CalcSettingForDivide(1, j, left, right);
                GenerateDivide(param, left, right, LineZ, idx);
            }

            /* Процедура расстановки узлов в глобальный массив */
            int32_t startIdx = Getlevel(i, 0); // Стартовый индекс для прохода по массиву
            for (int32_t k = 0; k < GlobalNz; k++)
            {
                // Скачки будут ровно на величину GlobalNx
                Grid[static_cast<uint64_t>(startIdx)].z = LineZ[static_cast<uint64_t>(k)];
                startIdx += GlobalNx;
            }
        }

        /************************************************************************/
        /*
                    Кратко Алгоритм:
                    Работаем с осью Z соответственно индексация будет происходить по этой оси
                    в цикле идем по всем столбцам массива сетки

                    Нужно получить левую и правую границу на каждом интервале
                    Сформировать массив отрезков по  данной координате
                    Занести полученный массив в Глобальную сетку
                */
        /* Цикл по всем горизонтальным линиям */
        for (int32_t i = 0; i < GlobalNx; i++)
        {
            int32_t idx = 0;
            /* Цикл по интеравалам оси Z */
            for (int32_t j = 0; j < Nz - 1; j++)
            {
                int32_t startIdx = i + GlobalNx * Getlevel(j, 1);
                int32_t endIdx = i + GlobalNx * Getlevel(j + 1, 1);
                double left = Grid[static_cast<uint64_t>(startIdx)].x; // Левая граница по х
                double right = Grid[static_cast<uint64_t>(endIdx)].x;  // Правая граница по х

                // Разбиение интервала подчиняется разбиению по оси z
                SettingForDivide param = CalcSettingForDivide(1, j, left, right);
                GenerateDivide(param, left, right, LineZ, idx);
            }
            /* Занесение результата в Итоговый массив */

            int32_t startIdx = i; // Стартовая позиция
            for (int32_t k = 0; k < GlobalNz; k++)
            {
                Grid[static_cast<uint64_t>(startIdx)].x = LineZ[static_cast<uint64_t>(k)];
                startIdx += GlobalNx;
            }
        }

        /*
            Генерация вспомогательных горизонтальных линий
            Цикл по всем горизонтальным линиям
        */
        for (int32_t i = 0; i < GlobalNz; i++)
        {
            int32_t idx = 0;
            for (int32_t j = 0; j < Nx - 1; j++)
            {
                int32_t startIdx = Getlevel(j, 0) + i * GlobalNx;
                int32_t endIdx = Getlevel(j + 1, 0) + i * GlobalNx;
                double left = Grid[static_cast<uint64_t>(startIdx)].z;
                double right = Grid[static_cast<uint64_t>(endIdx)].z;
                // Разбиение интервала подчиняется разбиению по оси x
                SettingForDivide param = CalcSettingForDivide(0, j, left, right);
                GenerateDivide(param, left, right, LineX, idx);
            }
            /* Занесение результатов в Глобальную сетку */
            int32_t startIdx = i * GlobalNx;
            for (int32_t k = 0; k < GlobalNx; k++)
            {
                Grid[static_cast<uint64_t>(startIdx)].z = LineX[static_cast<uint64_t>(k)];
                startIdx++;
            }
        }
        /************************************************************************/

        /*
            Разбиение по y теражируем сечение с учетом сечения плоскостью
            Создадим разбиение элементов массива Y
        */

        int32_t idxY = 0;
        for (int32_t i = 0; i < Ny - 1; i++)
        {
            double left = BaseGridY[static_cast<uint64_t>(i)];
            double right = BaseGridY[static_cast<uint64_t>(i + 1)];
            SettingForDivide param = CalcSettingForDivide(2, i, left, right);
            GenerateDivide(param, left, right, LineY, idxY);
        }

        // Вставляем элемент Y0 в базовом сечении
        int32_t idx = 0;
        for (int32_t j = 0; j < GlobalNz; j++)
        {
            for (int32_t k = 0; k < GlobalNx; k++)
            {
                // // Инициализация битовых полей без этой команды будет Ошибка !!! ( В актуальной версии исправлен этот недочет )
                // PointInfo::ClearInfo(Grid[idx].info);
                Grid[static_cast<uint64_t>(idx)].y = LineY[0];
                idx++;
            }
        }

        // Тиражирование сетки по всем узлам
        // i отвечает за индексацию в массиве по Y
        // Idx за индексацию по элементам

        for (int32_t i = 1; i < GlobalNy; i++)
        {
            int32_t ListIdx = 0;
            for (int32_t j = 0; j < GlobalNz; j++)
            {
                for (int32_t k = 0; k < GlobalNx; k++)
                {
                    // Инициализация битовых полей без этой команды будет Ошибка !!! ( В актуальной версии исправлен этот недочет )
                    // PointInfo::ClearInfo(Grid[idx].info);
                    Grid[static_cast<uint64_t>(idx)].y = LineY[static_cast<uint64_t>(i)];
                    Grid[static_cast<uint64_t>(idx)].x = Grid[static_cast<uint64_t>(ListIdx)].x;
                    Grid[static_cast<uint64_t>(idx)].z = Grid[static_cast<uint64_t>(ListIdx)].z;
                    ListIdx++;
                    idx++;
                }
            }
        }
        /************************************************************************/
    }
    catch (const std::bad_alloc &e)
    {
        Status.SetStatus(State::MEMORY_ALLOC_ERROR, e.what());
        // exit(State::MEMORY_ALLOC_ERROR);
    }
    catch (...)
    {
        Status.SetStatus(State::UNKNOWN_ERROR, "Unknown error in Grid3D_StreightQuadPrismatic::GenerateBaseGrid\n");
    }
}

/* Public */

GridStatus Grid3D_Integrate_StreightQuadPrismatic::GenerateGrid() noexcept
{
    // Расчет общего количества узлов в сетке
    // Для этого пройдемся по массиву разбиения каждого отрезка и вычислим общее число узлов
    GetTotalNumberOfNodes();

    // Генерация базовой сетки
    if (Status.GetState() == State::OK)
        GenerateBaseGrid();
    else
        return Status;

    return Status;
}

GridStatus Grid3D_Integrate_StreightQuadPrismatic::DivideGrid(const int coef) noexcept
{
    /* Для этого нужно сделать преобразование базовой сетки на заданный коэффициент и пересчитать всю сетку*/

    for (uint64_t i = 0; i < static_cast<uint64_t>(baseGrid.Nx - 1); i++)
    {
        baseGrid.DivideParam[0][i].num *= coef;
        baseGrid.DivideParam[0][i].coef = pow(baseGrid.DivideParam[0][i].coef, 1.0 / (static_cast<double>(coef)));
    }

    for (uint64_t i = 0; i < static_cast<uint64_t>(baseGrid.Nz - 1); i++)
    {
        baseGrid.DivideParam[1][i].num *= coef;
        baseGrid.DivideParam[1][i].coef = pow(baseGrid.DivideParam[1][i].coef, 1.0 / (static_cast<double>(coef)));
    }

    for (uint64_t i = 0; i < static_cast<uint64_t>(baseGrid.Ny - 1); i++)
    {
        baseGrid.DivideParam[2][i].num *= coef;
        baseGrid.DivideParam[2][i].coef = pow(baseGrid.DivideParam[2][i].coef, 1.0 / (static_cast<double>(coef)));
    }

    return Status;
}

GridStatus Grid3D_Integrate_StreightQuadPrismatic::ReGenerateGrid() noexcept
{
    /* Предварительно выставляем по дефолту все параметры класса */
    ElementCount = Dim = GlobalNx = GlobalNy = GlobalNz = 0;

    Grid.clear();

    /* Перегенерация сетки */
    Status = GenerateGrid();

    return Status;
}

StreightQuadPrismatic Grid3D_Integrate_StreightQuadPrismatic::GetElement(const int32_t idx) const noexcept
{
    int32_t Nx = GlobalNx;
    int32_t Ny = GlobalNy;
    int32_t Nz = GlobalNz;

    /*
        @param int32_t idx - Индекс конечного элемента
        @return int32_t\\
        @return  Вернет по номеру конечного элемента значение его праой нижней границе (локальный номер 1)
    */
    auto K = [Nx, Ny, Nz](int32_t _idx) -> int32_t
    {
        int32_t shiftXZ = 0;
        int32_t projidx = _idx % (((Nx - 1) * (Nz - 1)));
        if (projidx < Nx - 1)
        {
            shiftXZ = projidx;
        }
        else
        {
            int32_t level = static_cast<int32_t>(std::floor(projidx / (Nx - 1)));
            int32_t start = level * Nx;
            int32_t shift = projidx - (Nx - 1) * level;
            shiftXZ = start + shift;
        }

        int32_t levelY = static_cast<int32_t>(std::floor((_idx) / ((Nx - 1) * (Nz - 1))));
        int32_t shiftY = levelY * Nx * Nz;
        return shiftY + shiftXZ;
    };

    StreightQuadPrismatic Element;
    /* Расчет глобальных индексов */
    Element.GlobalIdx[0] = K(idx);
    Element.GlobalIdx[1] = Element.GlobalIdx[0] + 1;
    Element.GlobalIdx[2] = Element.GlobalIdx[0] + Nx;
    Element.GlobalIdx[3] = Element.GlobalIdx[1] + Nx;

    Element.GlobalIdx[4] = Element.GlobalIdx[0] + Nx * Nz;
    Element.GlobalIdx[5] = Element.GlobalIdx[4] + 1;
    Element.GlobalIdx[6] = Element.GlobalIdx[4] + Nx;
    Element.GlobalIdx[7] = Element.GlobalIdx[5] + Nx;

    for (int32_t i = 0; i < Element.ElementSize; i++)
    {
        Element.points[static_cast<uint64_t>(i)] = Grid[static_cast<uint64_t>(Element.GlobalIdx[static_cast<uint64_t>(i)])];
    }

    return Element;
}

GridStatus Grid3D_Integrate_StreightQuadPrismatic::SetBaseGrid(const BaseGrid3D_Integrate_StreightQuadPrismatic &baseGrid_) noexcept
{
    // Если структура не валидная то выбрасываем ошибку
    if (!baseGrid_.isReadyToUse)
        Status.SetStatus(State::LOAD_ERROR, "Incorrect Base Grid. Call from Grid3D_StreightQuadPrismatic(const BaseGrid3DStreightQuadPrismatic &baseGrid_)");
    else
        baseGrid = baseGrid_;
    return Status;
}

void Grid3D_Integrate_StreightQuadPrismatic::PrintGridSlice(int32_t level) const
{
uint64_t idx = static_cast<uint64_t>(level * GlobalNx * GlobalNz);
    cout << "Start idx: " << idx << "  End idx: " << idx + static_cast<uint64_t>(GlobalNx * GlobalNz - 1) << " Step Row: " << GlobalNx << "\n";
    cout << "Format = (x,z)\n";
    cout << "y = " << Grid[idx].y << "\n";
    for (int32_t i = 0; i < GlobalNz; i++)
    {
        for (int32_t j = 0; j < GlobalNx; j++)
        {
            // cout << fixed << std::setprecision(2) << "(" << grid.Grid[idx].x << ";" << grid.Grid[idx].z << ";" << grid.Grid[idx].y << ") ";
            printf("(%.2f;%.2f)", Grid[idx].x, Grid[idx].z);
            idx++;
        }
        cout << "\n";
    }
}
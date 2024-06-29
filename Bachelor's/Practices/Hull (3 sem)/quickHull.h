#include "lib.h"

struct Quick {
private:

    // Хранит результат в неотсортированном виде, т.е. точки выпуклой оболочки
    unordered_set<Pair, pair_hash> hull;

    // Функция находит расстояние от прямой p1p2 до точки p
    double LineDist(const Pair& p1, const Pair& p2, const Pair& p) {

        return fabs(LineL(p1, p2, p));
    }

    // Функция поиска точек искомой оболочки
    void QuickHull(const vector<Pair>& point, const Pair& p1, const Pair& p2, const int& side) {

        // текущий индикатор
        int currentI = -1;
        // максимальное расстояние
        double maxDist = 0;
        {
            // Расстояние от прямой p1p2 до точки point[i]
            double dist = 0;

            // Найдем точку на максимальном расстояни от прямой p1p2, находящуюся на стороне side
            for (int i = 0; i < point.size(); i++) {

                // Расстояние от прямой p1p2 до точки point[i]
                dist = LineDist(p1, p2, point[i]);

                if (detSide(p1, p2, point[i]) == side && dist > maxDist) {

                    currentI = i;
                    maxDist = dist;
                }
            }
        }

        // Если не одной точки не найдено, добавляем p1 и p2 в оболочку hull
        if (currentI == -1) {
            //cout << p1.first << " " << p1.second << endl;
            //cout << p2.first << " " << p2.second << endl;

            hull.insert(p1);
            hull.insert(p2);

            return;
        }

        // Рекурсия двух частей, разделенных по point[currentI]

        QuickHull(point, point[currentI], p1, -detSide(point[currentI], p1, p2));

        QuickHull(point, point[currentI], p2, -detSide(point[currentI], p2, p1));

        //cout << "-------------------\n";
    }

    // Функция вычисления точек оболочки
    void CalcHull(const vector<Pair>& point) {

        if (point.size() < 3) {

            cout << "Convex hull not possible, n < 3\n";
            return;
        }

        int min = 0, max = 0;

        for (int i = 1; i < point.size(); ++i) {

            //Найдем самую левую точку min и самую правую точку max
            if (point[i].first < point[min].first)

                min = i;

            if (point[i].first > point[max].first)

                max = i;

            // Если таких несколько, выберем среди таких нижнюю и верхнюю соответственно
            if ((fabs(point[i].first - point[min].first) < accuracy) && fabs(point[i].second - point[min].second) < accuracy)
                min = i;

            if ((fabs(point[i].first - point[max].first) < accuracy) && fabs(point[i].second - point[max].second) > accuracy)
                max = i;

        }

        // Найдем точки оболочки на одной стороне прямой point[min] point[max] с помощью рекурсии
        QuickHull(point, point[min], point[max], 1);

        // Найдем точки оболочки на противоположной стороне прямой point[min] point[max] с помощью рекурсии
        QuickHull(point, point[min], point[max], -1);
    }

    void exe() {

        //std::cout << "Quick:";
        vector<Pair> arr = { };
        ReadFromFile(arr, inFileQuick);
        CalcHull(arr);
        arr.clear();
        SortAndPrintHull(hull, outFileQuick);
    }


public:

    void operator()() {
        this->exe();
    }

};
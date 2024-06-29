#include "lib.h"

struct NonEffective {
private:

    // Хранит результат в неотсортированном виде, т.е. точки выпуклой оболочки
    vector<Pair> hull;

    // Функция вычисления точек оболочки
    void  calcHull(vector<Pair>& point) {


        if (point.size() < 3) {

            cout << "Convex hull not possible, n < 3\n";
            return;
        }

        int min = 0, max = 0;

        for (int i = 1; i < point.size(); ++i) {

            //Найдем самую правую точку max по x
            if (point[i].first > point[max].first)
                max = i;

            if ((fabs(point[i].first - point[max].first) < accuracy) && (point[i].second - point[max].second) > accuracy)
                max = i;
        }

        // добавим первый элемент в исходную оболочку
        iter_swap(point.begin(), point.begin() + max);
        hull.push_back(point[0]);

        // предыдущий элемент оболочки
        Pair prevPoint = point[0];
        for (auto it = point.begin() + 1; it < point.end(); ++it) {

            // Возможный следующий элемент в оболочке
            const Pair nextPoint = *it;

            // находятся ли точки на одной стороне прямой
            bool SameSide = true;
            auto iterator = point.begin();
            int side = 0;

            // найдем первую точку не лежащую на прямой
            for (; iterator < point.end(); ++iterator) {

                side = detSide(prevPoint, nextPoint, *iterator);
                if (side != 0) {
                    ++iterator;
                    break;
                }
            }

            // сравним оставшиеся точки с полученной (будут ли они по одну сторону относительно прямой prevPoint nextPoint)
            for (; iterator < point.end(); ++iterator) {

                int currPointSide = detSide(prevPoint, nextPoint, *iterator);
                if ((currPointSide != side) && (currPointSide != 0)) {
                    SameSide = false;
                    break;
                }
            }

            if (SameSide) {
                // если в hull уже есть такой элемент, то продолжаем цикл for
                if (find_if(hull.begin(),
                    hull.end(),
                    [nextPoint](const Pair& pointOfConvexHull) {
                        return ((nextPoint.first == pointOfConvexHull.first) && (nextPoint.second == pointOfConvexHull.second));
                    }) != hull.end()) {
                    continue;
                }

                hull.push_back(nextPoint);
                prevPoint = nextPoint;
                it = point.begin();
            }
        }
    }

    void exe() {

        //std::cout << "None effective:";
        vector<Pair> arr = { };
        ReadFromFile(arr, inFileNonEff);
        calcHull(arr);
        arr.clear();
        WriteToFile<vector<Pair>>(hull, outFileNonEff);
        //SortAndPrintHull(hull, outFileNonEff);
    }
public:

    void operator ()() {
        this->exe();
    }

};
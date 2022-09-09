#include <iostream> 
#include <unordered_set>
#include <set>
#include <fstream> 
#include <vector>
#include <algorithm>
#include <iomanip> 
using namespace std;

// Пары координат (x; y) 
#define Pair pair<double, double>
// Точность сравнения double
#define accuracy 10e-8
// Входной файл для Quick
#define inFileQuick "problem3.txt"
// Входной файл для NonEff
#define inFileNonEff "problem3.txt"
// Выходной файл Quick
#define outFileQuick "out_Quick.txt"
// Выходной файл Noneff
#define outFileNonEff "out_NonEff.txt"
// Количество значащих знаков после запятой
#define Precision 17

// hash for unorderd_set
struct pair_hash
{
	template <class T1, class T2>
	std::size_t operator () (std::pair<T1, T2> const& pair) const
	{
		std::size_t h1 = std::hash<T1>()(pair.first);
		std::size_t h2 = std::hash<T2>()(pair.second);

		return h1 ^ h2;
	}
};

double LineL(const Pair& p1, const Pair& p2, const Pair& p) {
	return ((p.second - p1.second) * (p2.first - p1.first) -
		(p2.second - p1.second) * (p.first - p1.first));
}

// Функция находит по какю сторону относительно прямой p1p2 лежит точка p
int detSide(const Pair& p1, const Pair& p2, const Pair& p)
{
	// Данное значение определяет по какую сторону от прямой p1p2 лежит точка p
	double side = LineL(p1, p2, p);

	if (side > 0)
		return 1;

	if (side < 0)
		return -1;

	return 0;
}

// Функия для записи координат точек из unordered_set<Pair> pr в текстовый файл (будет неотсортированный вариант)
template <typename T>
void WriteToFile(const T& OutputPoints, const string& path) {

	// Создадим ofstream для записи в файл
	ofstream fin;
	fin.open(path);
	fin << OutputPoints.size();
	fin << setprecision(Precision);

	for (auto it : OutputPoints)
		fin << "\n" << it.first << " " << it.second;

	fin.close();
}

// Функия для считывания координат точек из текстого файла в vector<Pair> pr
void ReadFromFile(vector<Pair>& InputPoints, const string& path) {
	// Создадим ifstream для чтения файла
	ifstream fout;

	fout.open(path);
	size_t n = 0;
	fout >> n;
	InputPoints.resize(n + 1);
	// cout << n << " " << pr.size();
	n = 0;
	// Записываем данные в vector<Pair> pr
	while (fout) {
		fout >> InputPoints[n].first >> InputPoints[n].second;
		++n;
	}

	fout.close();
}

// Функция для вывода результата в консоль (будет неотсортированный вариант)
template <typename T>
void PrintHull(const T& hull) {

	cout << "\nThe points in Convex Hull are: " << hull.size();
	cout << setprecision(Precision);

	for (auto it : hull)
	{
		cout << "\n(" << it.first << ", " << it.second << ") ";
	}
}

// Рассчет центра hull
Pair Calc_Center(const unordered_set<Pair, pair_hash>& hull) {

	double xc = 0.0;
	double yc = 0.0;

	for (auto it : hull) {

		xc += it.first;
		yc += it.second;
	}

	return { xc / hull.size(), yc / hull.size() };
}

void SortAndPrintHull(const unordered_set<Pair, pair_hash>& hull, const string& path) {

	// центр hull
	Pair c = Calc_Center(hull);

	// лямбда функция сравнения
	auto cmp = [c](const Pair a, const Pair b) {
		return atan2(a.first - c.first, a.second - c.second) < atan2(b.first - c.first, b.second - c.second);
	};

	// hull по часовой стрелке
	set<Pair, decltype(cmp)> ordHull(cmp);

	for (auto it : hull) {

		ordHull.insert(it);
	}

	WriteToFile<set<Pair, decltype(cmp)>>(ordHull, path);
	//PrintHull<set<Pair, decltype(cmp)>>(ordHull);
}
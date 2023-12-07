#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

const double EPS = 1e-16; // Машинный ноль

const string My_File_Pass = "C:\\Users\\pinev\\source\\repos\\ConsoleApplication2\\ConsoleApplication2\\";

const double M_PI = 3.14;
size_t NP = 0; // Количество узлов
size_t NE = 0; // Количество элементов
size_t NC = 0; // Количество контуров

vector<vector<size_t>> Elem; // Массив номеров элементов
vector<vector<double>> Points; // Массив координат узлов
vector<double> Square; // Массив площадей элементов
vector<double> FValues; // Массив значений функции в узлах
vector<vector<double>> DFEValues; // Массив значений производных функции на элементах
vector<vector<double>> DFNValues; // Массив значений производных функции в узлах
vector<double> IFEValues; // Массив значений интегралов функции на элементах
vector<double> IFNValues; // Массив значений интегралов функции в узлах

// Заданная функция
double f(const double x, const double y)
{
    return 10.0 + 5.0 * y * sin(2 * M_PI * x);
}

// Площадь треугольника по трем точкам
double S(const vector<double>& p0, const vector<double>& p1, const vector<double>& p2)
{
    return 0.5 * fabs((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]));
}

// Процедура чтения из файла
void My_Read_File(const string file)
{
    ifstream F(My_File_Pass + file + ".txt");

    F >> NE >> NP >> NC;

    Elem.resize(NE + 1);
    Points.resize(NP + 1);
    Square.resize(NE + 1);

    size_t pass = 0;

    for (size_t i = 1; i <= NE; ++i)
    {
        Elem[i].resize(3); // В каждом элементе 3 узла
        F >> pass >> pass >> Elem[i][0] >> Elem[i][1] >> Elem[i][2];
    }

    for (size_t i = 1; i <= NP; ++i)
    {
        Points[i].resize(3); // Каждая точка из R3
        F >> pass >> Points[i][0] >> Points[i][1] >> Points[i][2];
    }

    for (size_t i = 1; i <= NE; ++i)
        Square[i] = S(Points[Elem[i][0]], Points[Elem[i][1]], Points[Elem[i][2]]);

    F.close();
}

// Процедура вычисления значений функции в узлах и запись в файл
void My_Write_File_FValues(const string file)
{
    ofstream F(My_File_Pass + file + ".txt");

    FValues.resize(NP + 1);

    F << NP << "\n";

    double buf = 0.0;

    for (size_t i = 1; i <= NP; ++i)
    {
        buf = f(Points[i][0], Points[i][1]);
        FValues[i] = buf;

        F << i << " " << buf << "\n";
    }

    F.close();
}

// Процедура вычисления производных функции на элементах и запись в файл
void My_Write_File_DFEValues(const string fileX, const string fileY)
{
    ofstream FX(My_File_Pass + fileX + ".txt");
    ofstream FY(My_File_Pass + fileY + ".txt");

    DFEValues.resize(NE + 1);

    FX << NE << "\n";
    FY << NE << "\n";

    double buf = 0.0;
    size_t p0, p1, p2;

    for (size_t i = 1; i <= NE; ++i)
    {
        DFEValues[i].resize(2);

        p0 = Elem[i][0];
        p1 = Elem[i][1];
        p2 = Elem[i][2];

        buf = (Points[p1][1] - Points[p2][1]) * FValues[p0];
        buf += (Points[p2][1] - Points[p0][1]) * FValues[p1];
        buf += (Points[p0][1] - Points[p1][1]) * FValues[p2];

        buf /= 2.0 * Square[i];

        DFEValues[i][0] = buf;

        FX << i << " " << buf << "\n";

        buf = (Points[p2][0] - Points[p1][0]) * FValues[p0];
        buf += (Points[p0][0] - Points[p2][0]) * FValues[p1];
        buf += (Points[p1][0] - Points[p0][0]) * FValues[p2];

        buf /= 2.0 * Square[i];

        DFEValues[i][1] = buf;

        FY << i << " " << buf << "\n";
    }

    FX.close();
    FY.close();
}

// Процедура вычисления интегралов функции на элементах и запись в файл
void My_Write_File_IFEValues(const string file)
{
    ofstream F(My_File_Pass + file + ".txt");

    IFEValues.resize(NE + 1);

    F << NE << "\n";

    double buf = 0.0;

    for (size_t i = 1; i <= NE; ++i)
    {
        buf = (FValues[Elem[i][0]] + FValues[Elem[i][1]] + FValues[Elem[i][2]]) * Square[i] / 3.0;

        IFEValues[i] = buf;

        F << i << " " << buf << "\n";
    }

    F.close();
}

// Процедура вычисления производных функции в узлах и запись в файл
void My_Write_File_DFNValues(const string fileX, const string fileY)
{
    ofstream FX(My_File_Pass + fileX + ".txt");
    ofstream FY(My_File_Pass + fileY + ".txt");

    DFNValues.resize(NP + 1);

    vector<double> a(NP + 1, 0.0);
    auto b0 = a;
    auto b1 = a;

    double bufS = 1.0;
    double buf = 0.0;
    size_t n0, n1, n2;

    for (size_t i = 1; i <= NE; ++i)
    {
        n0 = Elem[i][0];
        n1 = Elem[i][1];
        n2 = Elem[i][2];

        bufS = Square[i];

        a[n0] += bufS;
        a[n1] += bufS;
        a[n2] += bufS;

        buf = bufS * DFEValues[i][0];

        b0[n0] += buf;
        b0[n1] += buf;
        b0[n2] += buf;

        buf = bufS * DFEValues[i][1];

        b1[n0] += buf;
        b1[n1] += buf;
        b1[n2] += buf;
    }

    FX << NP << "\n";
    FY << NP << "\n";

    for (size_t i = 1; i <= NP; ++i)
    {
        DFNValues[i].resize(2);

        buf = b0[i] / a[i];

        DFNValues[i][0] = buf;

        FX << i << " " << buf << "\n";

        buf = b1[i] / a[i];

        DFNValues[i][1] = buf;

        FY << i << " " << buf << "\n";
    }

    FX.close();
    FY.close();
}

// Процедура вычисления интегралов функции в узлах и запись в файл
void My_Write_File_IFNValues(const string file)
{
    ofstream F(My_File_Pass + file + ".txt");

    IFNValues.resize(NP + 1);

    vector<double> a(NP + 1, 0.0);
    auto b = a;

    double buf = 0.0;
    size_t n0, n1, n2;

    for (size_t i = 1; i <= NE; ++i)
    {
        n0 = Elem[i][0];
        n1 = Elem[i][1];
        n2 = Elem[i][2];

        buf = Square[i];

        a[n0] += buf;
        a[n1] += buf;
        a[n2] += buf;

        buf *= IFEValues[i];

        b[n0] += buf;
        b[n1] += buf;
        b[n2] += buf;
    }

    F << NP << "\n";

    for (size_t i = 1; i <= NP; ++i)
    {
        buf = b[i] / a[i];

        IFNValues[i] = buf;

        F << i << " " << buf << "\n";
    }

    F.close();
}

// Процедура вычисления функции, производных и интеграла в точке для проверки
void My_Check_Point(const vector<double> p, const string file)
{
    ofstream F(My_File_Pass + file + ".txt");

    F << p[0] << " " << p[1] << "\n";

    F.close();

    double t1, t2, t3;

    size_t a, b, c;

    size_t num = 0;

    for (size_t i = 1; i <= NE; ++i)
    {
        a = Elem[i][0];
        b = Elem[i][1];
        c = Elem[i][2];

        t1 = (Points[a][0] - p[0]) * (Points[b][1] - Points[a][1]) - (Points[a][1] - p[1]) * (Points[b][0] - Points[a][0]);
        t2 = (Points[b][0] - p[0]) * (Points[c][1] - Points[b][1]) - (Points[b][1] - p[1]) * (Points[c][0] - Points[b][0]);
        t3 = (Points[c][0] - p[0]) * (Points[a][1] - Points[c][1]) - (Points[c][1] - p[1]) * (Points[a][0] - Points[c][0]);

        if (t1 * t2 > -EPS && t1 * t3 > -EPS)
            num = i;
    }

    a = Elem[num][0];
    b = Elem[num][1];
    c = Elem[num][2];

    double bufS = S(Points[a], Points[b], Points[c]);

    double Na = S(p, Points[b], Points[c]);
    double Nb = S(p, Points[a], Points[c]);
    double Nc = S(p, Points[a], Points[b]);

    double buf = 0.0;

    cout << "Proverka (5, 5): " << "\n";
    //cout << "pCheck = {" << p[0] << ", " << p[1] << "} \n";

    buf = (Na * FValues[a] + Nb * FValues[b] + Nc * FValues[c]) / bufS;

    cout << "FValue = " << buf << " \n";

    t1 = (Na * DFNValues[a][0] + Nb * DFNValues[b][0] + Nc * DFNValues[c][0]) / bufS;
    t2 = (Na * DFNValues[a][1] + Nb * DFNValues[b][1] + Nc * DFNValues[c][1]) / bufS;

    cout << "Deriv = {" << t1 << ", " << t2 << "} \n";

    buf = (Na * IFNValues[a] + Nb * IFNValues[b] + Nc * IFNValues[c]) / bufS;

    cout << "Integral = " << buf << " \n";
}

int main()
{
    My_Read_File("output_4_0_2");

    // Точка для проверки
    const vector<double> pCheck = { 5.0, 5.0, 0.0 };

    My_Write_File_FValues("FValues");
    My_Write_File_DFEValues("DFEValuesX", "DFEValuesY");
    My_Write_File_IFEValues("IFEValues");
    My_Write_File_DFNValues("DFNValuesX", "DFNValuesY");
    My_Write_File_IFNValues("IFNValues");

    My_Check_Point(pCheck, "pCheck");

    return 0;
}

#define _USE_MATH_DEFINES
#include <cmath>
#include "CSR_Matrix.h"
#include "Solver.h"
#include <iostream>

const double p = 1;
const double s = 1;

const unsigned int NX = 71;
const unsigned int NY = 71;

double func(double t, double x, double y) { return 0; }

double func_0(double x, double y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
}

int main() {
    Heat_Equation_2D<double, NX, NY> eq;
    eq.coef = 1;
    eq.func = func;
    Grid2D<double, NX, NY> grid(0, p, 0, s);
    std::fill(grid.points.begin(), grid.points.end(), 0);
    for (int i = 1; i < NY - 1; ++i) {
        for (int j = 1; j < NX - 1; ++j) {
            grid(i, j) =
                func_0(grid.X_LEFT + j * grid.DX, grid.Y_BOTTOM + i * grid.DY);
        }
    }

    std::ofstream file;
    file.open("data.txt");
    Solve<double, NX, NY>(eq, grid, 0.5, 0.005, file);
}
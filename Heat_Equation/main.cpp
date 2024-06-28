#include "Crank_Nicolson.h"
#include "Tridiagonal_Matrix.h"
#include <iostream>

double func(double t, double x) { return t * (x + 1); }
double left_func(double t) { return t * t; }
double right_func(double t) { return t * t; }

const int NX = 85;

int main() {
    std::ofstream file;
    file.open("data100.txt");
    Heat_Equation eq;
    eq.coef = 1;
    eq.func = func;
    eq.left_boundary_conditions = { 0, 1 };
    eq.left_func = left_func;
    eq.right_boundary_conditions = { 1, 0 };
    eq.right_func = right_func;

    Grid1D<double, NX> grid(0, 1);
    for (int i = 0; i < NX; ++i) {
        grid[i] = 0;
    }

    Solve<double, NX>(eq, grid, 10, 0.01, file);
}
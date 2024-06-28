#define _USE_MATH_DEFINES
#include <cmath>
#include "Convection_Equation_Solver.h"

const double H = 0.5;
const double X_LEFT = 0;
const double X_RIGHT = 20;
const unsigned int NX = (X_RIGHT - X_LEFT) / H + 1;

double left_corner(const double m1, const double m2, const double CFL) {
    return 0;
}


int main() {
    Grid1D<NX> grid(X_LEFT, X_RIGHT);
    for (std::size_t i = 0; i < NX; ++i) {
        grid.points[i] =
            std::sin(4 * M_PI * (X_LEFT + H * i) / (X_RIGHT - X_LEFT));
    }

    std::ofstream file;
    Convection_Equation_Solver_Left_Corner<NX> solver1(1, grid);
    Convection_Equation_Solver_Left_Corner<NX> solver2(1.01, grid);
    Convection_Equation_Solver_Left_Corner<NX> solver3(0.6, grid);
    Convection_Equation_Solver_Lax_Wendroff<NX> solver4(1, grid);
    Convection_Equation_Solver_Lax_Wendroff<NX> solver5(1.1, grid);
    Convection_Equation_Solver_Lax_Wendroff<NX> solver6(0.6, grid);

    file.open("data_lc_CFL1.txt");
    solver1.Solve(0, 20, file);
    file.close();

    file.open("data_lc_CFL101.txt");
    solver2.Solve(0, 20, file);
    file.close();

    file.open("data_lc_CFL06.txt");
    solver3.Solve(0, 20, file);
    file.close();

    file.open("data_lw_CFL1.txt");
    solver4.Solve(0, 20, file);
    file.close();

    file.open("data_lw_CFL101.txt");
    solver5.Solve(0, 70, file);
    file.close();

    file.open("data_lw_CFL06.txt");
    solver6.Solve(0, 20, file);
    file.close();
}
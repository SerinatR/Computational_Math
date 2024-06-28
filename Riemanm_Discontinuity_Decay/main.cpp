#include "Solver.h"
#include <fstream>
#include <iostream>

int main() {
    Grid1D<Cell, 100> grid(-10, 10);
    for (std::size_t x = 0; x < 50; ++x) {
        // grid[x] = Cell({13, 0, 115385});
        grid[x].density = 13;
        grid[x].velocity = -30;
        grid[x].energy = 115385;
        // grid[99 - x] = Cell({1.3, 0, 115385});
        grid[99 - x].density = 13;
        grid[99 - x].velocity = 30;
        grid[99 - x].energy = 115385;
    }
    std::ofstream file;
    file.open("data_shockwaves.txt");
    Solve_Riemann_Discontinuity_Decay<100>(grid, 5. / 3, 0, 0.02, 1e-5, 0.01, file);
}
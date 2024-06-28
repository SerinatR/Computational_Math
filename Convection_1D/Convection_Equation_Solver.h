#include <array>
#include <fstream>

template <unsigned int NX> struct Grid1D {
    const double X_LEFT;
    const double X_RIGHT;
    const double DX;
    std::array<double, NX> points;

    Grid1D(const double x_left, const double x_right) : X_LEFT(x_left), X_RIGHT(x_right), DX((X_RIGHT - X_LEFT) / (NX - 1)) {}
};

template <unsigned int NX> class Convection_Equation_Solver {
protected:
    const double CFL;
    Grid1D<NX> grid;

    virtual void Update_Grid() = 0;

    void Write_Grid_To_File(std::ofstream& file) const {
        for (std::size_t i = 0; i < NX; ++i) {
            file << grid.X_LEFT + grid.DX * i << " " << grid.points[i]
                << std::endl;
        }
    }

public:
    Convection_Equation_Solver(const double cfl, const double x_left, const double x_right) : CFL(cfl), grid(x_left, x_right) {}

    Convection_Equation_Solver(double cfl, const Grid1D<NX>& grid) : CFL(cfl), grid(grid) {}

    void Solve(const double t_start, const double t_end, std::ofstream& file) {
        Write_Grid_To_File(file);
        file << std::endl << std::endl;
        for (double t = t_start; t < t_end; t += CFL * grid.DX) {
            Update_Grid();
            Write_Grid_To_File(file);
            file << std::endl << std::endl;
        }
    }
};

template <unsigned int NX>
class Convection_Equation_Solver_Left_Corner : public Convection_Equation_Solver<NX> {
    using Convection_Equation_Solver<NX>::CFL;
    using Convection_Equation_Solver<NX>::grid;

public:
    Convection_Equation_Solver_Left_Corner(const double cfl, const double x_left, const double x_right)
        : Convection_Equation_Solver<NX>(cfl, x_left, x_right) {}

    Convection_Equation_Solver_Left_Corner(double cfl, const Grid1D<NX>& grid)
        : Convection_Equation_Solver<NX>(cfl, grid) {}

    void Update_Grid() override {
        std::array<double, NX> points_new;
        for (std::size_t x = 1; x < NX; ++x) {
            points_new[x] =
                grid.points[x] - CFL * (grid.points[x] - grid.points[x - 1]);
        }
        points_new[0] = points_new[NX - 1];
        grid.points = points_new;
    }
};

template <unsigned int NX>
class Convection_Equation_Solver_Lax_Wendroff : public Convection_Equation_Solver<NX> {
    using Convection_Equation_Solver<NX>::CFL;
    using Convection_Equation_Solver<NX>::grid;

public:
    Convection_Equation_Solver_Lax_Wendroff(const double cfl, const double x_left, const double x_right)
        : Convection_Equation_Solver<NX>(cfl, x_left, x_right) {}

    Convection_Equation_Solver_Lax_Wendroff(double cfl, const Grid1D<NX>& grid)
        : Convection_Equation_Solver<NX>(cfl, grid) {}

    void Update_Grid() override {
        std::array<double, NX> points_new;
        for (std::size_t x = 1; x < NX - 1; ++x) {
            points_new[x] =
                grid.points[x] -
                CFL * (grid.points[x + 1] - grid.points[x - 1]) / 2 +
                CFL * CFL *
                (grid.points[x + 1] - 2 * grid.points[x] +
                    grid.points[x - 1]) /
                2;
        }
        points_new[NX - 1] = grid.points[NX - 1] -
            CFL / 2 * (grid.points[1] - grid.points[NX - 2]) +
            CFL * CFL / 2 *
            (grid.points[1] - 2 * grid.points[NX - 1] +
                grid.points[NX - 2]);
        points_new[0] = points_new[NX - 1];
        grid.points = points_new;
    }
};

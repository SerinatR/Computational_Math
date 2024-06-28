#include <array>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <iostream>

struct Cell {
    double density;
    double velocity;
    double energy;

    Cell() : density(0), energy(0), velocity(0) {}

    Cell(const Eigen::Matrix<double, 3, 1>& vec)
        : density(vec(0)), velocity(vec(1) / vec(0)), energy(vec(2) / vec(0)) {}

    operator Eigen::Vector3d() const {
        return Eigen::Vector3d(density, velocity * density, energy * density);
    }
};

template <typename T, unsigned int NX> struct Grid1D {
    const double X_LEFT;
    const double X_RIGHT;
    const double DX;
    std::array<T, NX> points;

    Grid1D(const double x_left, const double x_right)
        : X_LEFT(x_left), X_RIGHT(x_right), DX((X_RIGHT - X_LEFT) / (NX - 1)) {}

    T& operator[](std::size_t index) { return points[index]; }

    const T& operator[](std::size_t index) const { return points[index]; }

    friend std::ostream& operator<<(std::ostream& os,
        const Grid1D<Cell, NX>& grid) {
        for (std::size_t i = 0; i < NX; ++i) {
            double x = grid.X_LEFT + i * grid.DX;
            os << x << " " << grid[i].density << " " << grid[i].velocity << " "
                << grid[i].energy << std::endl;
        }
        return os;
    }
};

template <unsigned int NX>
double Get_Max_Abs_Eigenvalue(const Grid1D<Cell, NX>& grid, const double c) {
    double max_abs_eigenvalue = 0;
    for (const auto& element : grid.points) {
        double abs_velocity_plus = std::abs(element.velocity + c);
        double abs_velocity_minus = std::abs(element.velocity - c);

        max_abs_eigenvalue = std::max(max_abs_eigenvalue, abs_velocity_plus);
        max_abs_eigenvalue = std::max(max_abs_eigenvalue, abs_velocity_minus);
    }
    return max_abs_eigenvalue;
}

template <unsigned int NX>
void Solve_Riemann_Discontinuity_Decay(Grid1D<Cell, NX> grid, const double adiabatic_index,
    const double t_start, const double t_end,
    double t_increment, const double cfl_max,
    std::ostream& file) {
    double t = t_start;
    double speed_of_sound;
    double t_incr = t_increment;
    Grid1D<Cell, NX> grid_next = grid;
    Eigen::Matrix<double, 3, 3> omega_T{{0, 0, adiabatic_index - 1},
        { 0, 0, adiabatic_index - 1 },
        { 0, 0, adiabatic_index - 1 }};
    Eigen::Matrix<double, 3, 3> lambda = Eigen::Matrix<double, 3, 3>::Zero();
    Eigen::Matrix<double, 3, 3> A;

    file << "# Time: " << t << std::endl;
    file << grid;
    file << std::endl << std::endl;

    while (t < t_end) {
        double max_abs_eigenvalue =
            get_max_abs_eigenvalue<NX>(grid, speed_of_sound);
        if (t_incr > grid.DX / max_abs_eigenvalue) {
            t_incr = grid.DX * cfl_max / max_abs_eigenvalue;
        }
        for (std::size_t x = 1; x < NX - 1; ++x) {
            speed_of_sound = std::sqrt(adiabatic_index * (adiabatic_index - 1) *
                grid[x].energy);

            omega_T(0, 0) = -grid[x].velocity * speed_of_sound;
            omega_T(0, 1) = speed_of_sound;
            omega_T(1, 0) = -speed_of_sound * speed_of_sound;
            omega_T(2, 0) = grid[x].velocity * speed_of_sound;
            omega_T(2, 1) = -speed_of_sound;

            lambda(0, 0) = grid[x].velocity + speed_of_sound;
            lambda(1, 1) = grid[x].velocity;
            lambda(2, 2) = grid[x].velocity - speed_of_sound;

            A = omega_T.inverse() * lambda * omega_T;

            grid_next[x] = Cell(Eigen::Vector3d(grid[x]) -
                t_incr * A *
                (Eigen::Vector3d(grid[x + 1]) -
                    Eigen::Vector3d(grid[x - 1])) /
                (2 * grid.DX) +
                t_incr *
                (omega_T.inverse() *
                    lambda.array().abs().matrix() * omega_T) *
                (Eigen::Vector3d(grid[x + 1]) -
                    2 * Eigen::Vector3d(grid[x]) +
                    Eigen::Vector3d(grid[x - 1])) /
                (2 * grid.DX));
        }
        grid_next[0] = grid_next[1];
        grid_next[NX - 1] = grid_next[NX - 2];
        grid.points = grid_next.points;
        file << "# Time: " << t << std::endl;
        file << grid;
        file << std::endl << std::endl;
        t += t_incr;
        t_incr = std::max(t_increment, t_incr);
    }
}
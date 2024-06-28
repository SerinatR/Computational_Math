#include <array>
#include <vector>

#ifndef TRIDIAGONAL_MATRIX
#define TRIDIAGONAL_MATRIX

template <typename T, unsigned int N> class Tridiagonal_Matrix {
private:
    std::array<std::array<T, 3>, N> m_data;

public:
    Tridiagonal_Matrix(const std::array<std::array<T, 3>, N>& data)
        : m_data{ data } {}

    T operator()(const int i, const int j) const {
        if (j == i - 1) {
            return m_data[i][0];
        }
        else if (j == i) {
            return m_data[i][1];
        }
        else if (j == i + 1) {
            return m_data[i][2];
        }
        else {
            return 0;
        }
    }

    std::array<T, N> Solve(const std::array<T, N>& column) const {
        int n = N - 1;
        std::vector<T> p_vector{0}, q_vector{ 0 };
        p_vector.reserve(N);
        q_vector.reserve(N);
        for (int i = 0; i < n + 1; ++i) {
            p_vector.push_back(
                -1 * (*this)(i, i + 1) /
                ((*this)(i, i - 1) * p_vector[i] + (*this)(i, i)));
            q_vector.push_back(
                (column[i] - (*this)(i, i - 1) * q_vector[i]) /
                ((*this)(i, i - 1) * p_vector[i] + (*this)(i, i)));
        }
        std::array<T, N> solution;
        solution[n] = (column[n] - (*this)(n, n - 1) * q_vector[n]) /
            ((*this)(n, n - 1) * p_vector[n] + (*this)(n, n));
        for (int i = n - 1; i >= 0; --i) {
            solution[i] = p_vector[i + 1] * solution[i + 1] + q_vector[i + 1];
        }
        return solution;
    }
};

#endif

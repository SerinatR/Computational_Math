#include "Vec_Math.h"
#include <iomanip>
#include <iostream>
#include <map>
#include <vector>

#ifndef CSR_MATRIX
#define CSR_MATRIX

template <typename T> class CSR_Matrix {
private:
    std::vector<T> m_values;
    std::vector<int> m_column_indexes, m_row_indexation;

public:
    CSR_Matrix(const std::map<std::pair<int, int>, T>& data) {
        m_values.reserve(data.size());
        m_column_indexes.reserve(data.size());
        m_row_indexation.resize(data.rbegin()->first.first + 2, 0);
        int last_row = -1;
        int entry_count = 0;
        for (const auto& entry : data) {
            int current_row = entry.first.first;
            int current_col = entry.first.second;
            const T& value = entry.second;
            m_values.push_back(value);
            m_column_indexes.push_back(current_col);
            if (current_row != last_row) {
                for (int i = last_row + 1; i <= current_row; ++i) {
                    m_row_indexation[i] = entry_count;
                }
                last_row = current_row;
            }
            ++entry_count;
        }
        m_row_indexation[last_row + 1] = entry_count;
        for (int i = last_row + 2; i < m_row_indexation.size(); ++i) {
            m_row_indexation[i] = entry_count;
        }
    }
    const std::vector<T>& get_values() const { return m_values; }
    const std::vector<int>& get_column_indexes() const {
        return m_column_indexes;
    }
    const std::vector<int>& get_row_indexation() const {
        return m_row_indexation;
    }
    T operator()(int i, int j) const {
        if (i + 2 > m_row_indexation.size()) {
            return 0;
        }
        for (int k = m_row_indexation[i]; k < m_row_indexation[i + 1]; ++k) {
            if (m_column_indexes[k] == j) {
                return m_values[k];
            }
        }
        return 0;
    }
    std::vector<T> operator*(const std::vector<T>& column) const {
        std::vector<T> res;
        res.reserve(m_row_indexation.size() - 1);
        for (int k = 0; k < m_row_indexation.size() - 1; ++k) {
            T temp = 0;
            for (int i = m_row_indexation[k]; i < m_row_indexation[k + 1];
                ++i) {
                temp += m_values[i] * column[m_column_indexes[i]];
            }
            res.push_back(temp);
        }
        return res;
    }

    void Print() const {
        int num_rows = m_row_indexation.size() - 1;

        for (int row = 0; row < num_rows; ++row) {
            int start = m_row_indexation[row];
            int end = m_row_indexation[row + 1];
            int col_index = 0;
            for (int idx = start; idx < end; ++idx) {
                while (col_index < m_column_indexes[idx]) {
                    std::cout << std::setw(5) << 0;
                    col_index++;
                }
                std::cout << std::setw(5) << m_values[idx];
                col_index++;
            }
            while (col_index < num_rows) {
                std::cout << std::setw(5) << 0;
                col_index++;
            }
            std::cout << std::endl;
        }
    }
};

template <typename T>
std::vector<T>
Conjugate_Gradients(const CSR_Matrix<T>& A, const std::vector<T>& b,
    const std::vector<T>& x0, const double tolerance) {
    std::vector<double> x = x0;
    std::vector<double> residual = A * x - b;
    std::vector<double> d = residual;
    std::vector<double> residual_next;
    double alpha, beta;
    while (length(residual) > tolerance) {
        alpha = dot_product(residual, residual) / dot_product(d, A * d);
        x -= alpha * d;
        residual_next = residual - alpha * (A * d);
        if (length(d) == 0) {
            break;
        }
        else {
            beta = dot_product(residual_next, residual_next) /
                dot_product(residual, residual);
            d = residual_next + beta * d;
        }
        residual = residual_next;
    }
    return x;
}

#endif
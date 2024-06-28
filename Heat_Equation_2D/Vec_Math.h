#include <cassert>
#include <cmath>
#include <vector>

#ifndef VECTOR_OPERATIONS
#define VECTOR_OPERATIONS

template <typename T>
std::vector<T>& operator+=(std::vector<T>& vec1, const std::vector<T>& vec2) {
    assert(vec1.size() == vec2.size() && "Vectors must be equal length");
    for (std::size_t i = 0; i < vec1.size(); ++i) {
        vec1[i] += vec2[i];
    }
    return vec1;
}

template <typename T>
std::vector<T> operator+(const std::vector<T>& vec1,
    const std::vector<T>& vec2) {
    std::vector<T> vec_res(vec1);
    vec_res += vec2;
    return vec_res;
}

template <typename T>
std::vector<T>& operator*=(std::vector<T>& vec, const T& num) {
    for (std::size_t i = 0; i < vec.size(); ++i) {
        vec[i] *= num;
    }
    return vec;
}

template <typename T>
std::vector<T> operator*(const std::vector<T>& vec, const T& num) {
    std::vector<T> vec_res(vec);
    vec_res *= num;
    return vec_res;
}

template <typename T>
std::vector<T>& operator-=(std::vector<T>& vec1, const std::vector<T>& vec2) {
    assert(vec1.size() == vec2.size() && "Vectors must be equal length");
    for (std::size_t i = 0; i < vec1.size(); ++i) {
        vec1[i] -= vec2[i];
    }
    return vec1;
}

template <typename T>
std::vector<T> operator-(const std::vector<T>& vec1,
    const std::vector<T>& vec2) {
    std::vector<T> vec_res(vec1);
    vec_res -= vec2;
    return vec_res;
}

template <typename T>
std::vector<T> operator*(const T& num, const std::vector<T>& vec) {
    return vec * num;
}

template <typename T>
T dot_product(const std::vector<T>& vec1, const std::vector<T>& vec2) {
    assert(vec1.size() == vec2.size() && "Vectors must be equal length");
    T res = 0;
    for (std::size_t i = 0; i < vec1.size(); ++i) {
        res += vec1[i] * vec2[i];
    }
    return res;
}

double length(const std::vector<double>& vec) {
    return std::sqrt(dot_product(vec, vec));
}

#endif
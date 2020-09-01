#pragma once

#include "base_constant_function.hpp"

#include <iostream>
#include <iomanip>

#include <boost/type_traits/enable_if.hpp>
#include <boost/multi_array.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/random.hpp>
#include <boost/any.hpp>

///2-х мерная матрица.
///
template<size_t N, size_t M, typename T>
struct matrix_2 ;
template<size_t N, size_t M, typename T>
std::ostream& operator << (std::ostream& os, const matrix_2<N, M, T>& A);

template<size_t N, size_t M, typename T>
struct matrix_2 : boost::multi_array<T, 2> {
    static constexpr size_t dimension = 2;
    typedef boost::multi_array<T, 2> type;
    typedef T value_type;
    type MTRX;
    matrix_2() : MTRX({boost::extents[N][M]}) {}
    void fill(const double x) {
        for (size_t i = 0; i < N; ++i)
            for (size_t j = 0; j < M; ++j) MTRX[i][j] = x;
    }
    size_t size(size_t i) const {
        BOOST_ASSERT(i < 2);
        return MTRX.shape()[i];
    }
    constexpr T& operator () (size_t i, size_t j) {
        return MTRX[i][j];
    }
    constexpr T const& operator () (size_t i, size_t j) const {
        return MTRX[i][j];
    }
    constexpr T& at (size_t i, size_t j) {
        return MTRX[i][j];
    }
    matrix_2<N, M, T>& operator = (matrix_2<N, M, T> const& other) {
        MTRX = other.MTRX;
        return *this;
    }
    friend std::ostream& operator << (std::ostream& os, const matrix_2<N, M, T>& A){
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                os << A.MTRX[i][j] << "\t";
            }
            os << std::endl;
        }
        return os;
    }
};

///З-х мерная матрица.
///
template<size_t N, size_t M, size_t K, typename T>
struct matrix_3 ;
template<size_t N, size_t M, size_t K, typename T>
std::ostream& operator << (std::ostream& os, const matrix_3<N, M, K, T>& A);

template<size_t N, size_t M, size_t K, typename T>
struct matrix_3 : boost::multi_array<T, 3> {
    static constexpr size_t dimension = 3;
private:
    typedef boost::multi_array<T, 3> type;
    typedef T value_type;
    type MTRX;
public:
    matrix_3() : MTRX({boost::extents[N][M][K]}) {}

    size_t size(size_t i) const {
        BOOST_ASSERT(i < 3);
        return MTRX.shape()[i];
    }
    constexpr T& operator () (size_t i, size_t j, size_t k) {
        return MTRX[i][j][k];
    }
    constexpr T const& operator () (size_t i, size_t j, size_t k) const {
        return MTRX[i][j][k];
    }
    constexpr T& at (size_t i, size_t j, size_t k) {
        return MTRX[i][j][k];
    }
    matrix_3<N, M, K, T>& operator = (matrix_3<N, M, K, T> const& other) {
        MTRX = other.MTRX;
        return *this;
    }
    friend std::ostream& operator << (std::ostream& os, const matrix_3<N, M, K, T>& A){
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                for (size_t k = 0; k < K; ++k) {
                    os << A.MTRX[i][j][k] << "\t";
                }
                os << std::endl;
            }
            os << std::endl;
        }
        return os;
    }
};

template<typename T> struct is_int : boost::mpl::false_ {};
template<> struct is_int<int> : boost::mpl::true_ {};
template<> struct is_int<unsigned> : boost::mpl::true_ {};


template<typename T, typename Matrix, typename = boost::enable_if_t<(Matrix::dimension > 0)>>
void gen_rand_matrix(Matrix& A, T min, T max) {
    std::time_t now = std::time(0);
    boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
        if constexpr(Matrix::dimension == 3 && is_int<T>::value) {
            boost::random::uniform_int_distribution<> dist{min, max};
            for(size_t i = 0; i < A.size(0); ++i)
                for(size_t j = 0; j < A.size(1); ++j)
                    for(size_t k = 0; k < A.size(2); ++k)
                        A(i, j, k) = dist(gen);
        }
        if constexpr(Matrix::dimension == 3 && !is_int<T>::value) {
            boost::random::uniform_real_distribution<> dist{min, max};
            for(size_t i = 0; i < A.size(0); ++i)
                for(size_t j = 0; j < A.size(1); ++j)
                    for(size_t k = 0; k < A.size(2); ++k)
                        A(i, j, k) = dist(gen);
        }
        if constexpr(Matrix::dimension == 2 && is_int<T>::value) {
            boost::random::uniform_int_distribution<> dist{min, max};
            for(size_t i = 0; i < A.size(0); ++i)
                for(size_t j = 0; j < A.size(1); ++j)
                        A(i, j) = dist(gen);
        }
        if constexpr(Matrix::dimension == 2 && !is_int<T>::value) {
            boost::random::uniform_real_distribution<> dist{min, max};
            for(size_t i = 0; i < A.size(0); ++i)
                for(size_t j = 0; j < A.size(1); ++j)
                        A(i, j) = dist(gen);
        }
}

///4-х мерная матрица.
///
template<size_t N, size_t M, size_t K, size_t L, typename T>
struct matrix_4 ;
template<size_t N, size_t M, size_t K, size_t L, typename T>
std::ostream& operator << (std::ostream& os, const matrix_4<N, M, K, L, T>& A);

template<size_t N, size_t M, size_t K, size_t L, typename T>
struct matrix_4 : boost::multi_array<T, 4> {
    static constexpr size_t dimension = 3;
private:
    typedef boost::multi_array<T, 4> type;
    typedef T value_type;
    type MTRX;
public:
    matrix_4() : MTRX({boost::extents[N][M][K][L]}) {}

    size_t size(size_t i) const {
        BOOST_ASSERT(i < 4);
        return MTRX.shape()[i];
    }
    constexpr T& operator () (size_t i, size_t j, size_t k, size_t l) {
        return MTRX[i][j][k][l];
    }
    constexpr T const& operator () (size_t i, size_t j, size_t k, size_t l) const {
        return MTRX[i][j][k][l];
    }
    constexpr T& at (size_t i, size_t j, size_t k, size_t l) {
        return MTRX[i][j][k][l];
    }
    matrix_4<N, M, K, L, T>& operator = (matrix_4<N, M, K, L, T> const& other) {
        MTRX = other.MTRX;
        return *this;
    }
    friend std::ostream& operator << (std::ostream& os, const matrix_4<N, M, K, L, T>& A){
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < M; ++j) {
                for (size_t k = 0; k < K; ++k) {
                    for (size_t l = 0; l < L; ++l) {
                        os << A.MTRX[i][j][k][l] << "\t";
                    }
                    os << std::endl;
                }
                os << std::endl;
            }
            os << std::endl;
        }
        return os;
    }
};
/*!
 *  Умножение матриц методом метапрограммирования
 */

/*!
 * struct abstract_sum
 */
template<class Closure>
struct abstract_sum_closures {
    typedef typename Closure::value_type value_type;
    abstract_sum_closures(Closure &closure) :  closure(closure), result(value_type()){}

    template<unsigned I>
    void apply(){
        result += closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_sums(Closure &closure) {
    abstract_sum_closures<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}

template <size_t I, size_t J, typename Matrix>
struct matrix_prod_closure
{
    typedef typename Matrix::value_type value_type;
private:
    const Matrix& A, B;
public:
    matrix_prod_closure(const Matrix& A, const Matrix& B) : A(A), B(B) {}
    template<size_t K>
    value_type value() const {
        return A(I, K) * B(K, J);
    }
};

template <size_t I, size_t J, typename Matrix>
struct result_matrx {
    typedef typename Matrix::value_type value_type;
    const Matrix& A, B;
    result_matrx(const Matrix& A, const Matrix& B): A(A), B(B) {}
    template <size_t K>
    value_type value() const {
        matrix_prod_closure<I, J, Matrix> closure(A, B);
        return abstract_sums<K + 1>(closure);
    }
};
template <size_t I, size_t J, typename Matrix>
typename Matrix::value_type calc_result_matrx(const Matrix& A, const Matrix& B) {
    matrix_prod_closure<I, J, Matrix> closure(A, B);
    return abstract_sums<BASE_FUNCTION_COUNT>(closure);
}

template <size_t i, typename Matrix>
struct A_m
{
private:
    const Matrix& A, B;
    Matrix& C;

public:
    A_m(Matrix& C, const Matrix& A, const Matrix& B): C(C), A(A), B(B){}

    template <size_t j>
    void apply() {
        C(i, j) = calc_result_matrx<i, j, Matrix>(A, B);
    }
};
template <typename Matrix>
struct CALC_A
{
private:
    const Matrix& A, B;
    Matrix& C;
public:
    CALC_A(Matrix& C, const Matrix& A, const Matrix& B): C(C), A(A), B(B){}

    template <size_t i>
    void apply() {
        A_m<i, Matrix> closure(C, A, B);
        meta_loop<i + 1>(closure);
    }
};
template <typename Matrix>
void calc_multy_matrix(Matrix& C, const Matrix& A, const Matrix& B) {
    CALC_A<Matrix> closure(C, A, B);
    meta_loop<BASE_FUNCTION_COUNT>(closure);
}

#pragma once

#include <boost/mpl/bool.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/type_traits/enable_if.hpp>
#include <boost/multi_array.hpp>

#include <array>

#define BASE_FUNCTION_COUNT 10

static constexpr std::array<double, 5> omega{0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885};
static constexpr std::array<double, BASE_FUNCTION_COUNT> X_i{-0.97390, -0.86506,
                                                             -0.67940, -0.43339, -0.14887, 0.14887, 0.43339, 0.67940, 0.86506, 0.97390};
static constexpr std::array<double, BASE_FUNCTION_COUNT> Y_i{-0.97390, -0.86506,
                                                             -0.67940, -0.43339, -0.14887, 0.14887, 0.43339, 0.67940, 0.86506, 0.97390};
static constexpr std::array<double, BASE_FUNCTION_COUNT> Z_i{-0.97390, -0.86506,
                                                             -0.67940, -0.43339, -0.14887, 0.14887, 0.43339, 0.67940, 0.86506, 0.97390};

///Factorial

template<size_t N>
struct factorial {
    static constexpr size_t value = N*factorial<N - 1>::value;
};
template<>
struct factorial<0> {
    static constexpr size_t value = 1;
};

///
///Вычисление степени
///

template<int N, typename T>
typename boost::enable_if_t<(N < 0), T> pow_(const T& x) {
    return T(1)/pow_<-N>(x);
}
template<int N, typename T>
typename boost::enable_if_t<(N == 0), T> pow_(const T& x) {
    return T(1);
}
template<int N, typename T>
typename boost::enable_if_t<(N > 0) && (N%2 == 0), T> pow_(const T& x) {
    T p = pow_<N / 2>(x);
    return p*p;
}
template<int N, typename T>
typename boost::enable_if_t<(N > 0) && (N%2 == 1), T> pow_(const T& x) {
    return pow_<N - 1>(x)*x;
}


/*!
 * struct meta_loop
 */
template <size_t N, size_t M, size_t I, size_t J, class Closure>
typename boost::enable_if_t<(I == N) && (J == M)> is_meta_loop2(Closure& closure) {}

template <size_t N, size_t M, size_t I, size_t J, class Closure>
typename boost::enable_if_t<(I < N) && (J < M)> is_meta_loop2(Closure& closure) {
    closure.template apply<I, J>();
    is_meta_loop2<N, M, I + 1, J + 1>(closure);
}

template <size_t N, size_t M, class Closure>
void meta_loop2(Closure& closure) {
    is_meta_loop2<N, M, 0, 0>(closure);
}
/*!
 * struct meta_loop
 */
template <size_t N, size_t I, class Closure>
typename boost::enable_if_t<(I == N)> is_meta_loop(Closure& closure) {}

template <size_t N, size_t I, class Closure>
typename boost::enable_if_t<(I < N)> is_meta_loop(Closure& closure) {
    closure.template apply<I>();
    is_meta_loop<N, I + 1>(closure);
}
template <size_t N, class Closure>
void meta_loop(Closure& closure) {
    is_meta_loop<N, 0>(closure);
}
template <size_t N, class Closure>
void meta_loopUV(Closure& closure) {
    is_meta_loop<N, 1>(closure);
}

/////Calculate Binom

template<size_t N, size_t K>
struct BC {
    static constexpr size_t value = factorial<N>::value / factorial<K>::value / factorial<N - K>::value;
};

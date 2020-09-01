#pragma once

#include "base_constant_function.hpp"
#include "matrix.hpp"

#include <iostream>
#include <algorithm>
#include <functional>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>

namespace mpl = boost::mpl;

template<typename T>
using abc_arr = std::array<T, 4>;

template<typename T>
struct  math_object_base {
    T& self() {
        return static_cast<T&>(*this);
    }
    const T& self() const {
        return static_cast<const T&>(*this);
    }
};

template<size_t ALPHA, size_t BETA, size_t GAMMA>
struct base_function {
    static constexpr size_t alpha = ALPHA;
    static constexpr size_t beta  = BETA;
    static constexpr size_t gamma = GAMMA;
};

using bf_t = mpl::vector10< base_function<0, 0, 0>,
                            base_function<1, 0, 0>,
                            base_function<0, 1, 0>,
                            base_function<0, 0, 1>,
                            base_function<2, 0, 0>,
                            base_function<0, 2, 0>,
                            base_function<0, 0, 2>,
                            base_function<1, 1, 0>,
                            base_function<1, 0, 1>,
                            base_function<0, 1, 1>>;


struct print_bf_t {
    template<class T>
    void operator() (T) const {
        std::cout << typeid (T).name() << std::endl;
    }
};

template<size_t N>
struct get_bf_type {
    typedef typename mpl::at_c<bf_t, N>::type type;
};

template<typename T>
struct ABC_ {
    typedef T value_type;
    ABC_(const abc_arr<T>& a_, const abc_arr<T>& b_, const abc_arr<T>& c_) : a(a_), b(b_), c(c_) {}

    template<size_t I>
    typename std::enable_if_t<I == 1, const abc_arr<T>>
    get() const { return a; }

    template<size_t I>
    typename std::enable_if_t<I == 2, const abc_arr<T>>
    get() const { return b; }

    template<size_t I>
    typename std::enable_if_t<I == 3, const abc_arr<T>>
    get() const { return c; }

private:
    const abc_arr<T> a, b, c;
};

//// monomial3d

template<size_t alpha, size_t beta, size_t gamma>
struct monomial3d {
    static constexpr size_t value =
            factorial<alpha + beta + gamma + 3>::value / factorial<alpha>::value / factorial<beta>::value /factorial<gamma>::value;
};

///abstract_sum
template<size_t N, size_t K, typename T, class action_class, typename A>
typename std::enable_if_t<(K > 0), T>
abstract_sum(const A& args) {
    return abstract_sum<N, K - 1, T, action_class, A>(args) + action_class::template value<N, K>(args);
}

template<size_t N, size_t K, typename T, class action_class, typename A>
typename std::enable_if_t<(K == 0), T>
abstract_sum(const A& args) {
    return action_class::template value<N, K>(args);
}

///Integrate, BC-binom_kooeff
template<size_t step, class next_step, size_t PX, size_t PY, size_t PZ>
struct tetra_integrate_polinomial_sum3 {
    template<size_t J, size_t K, typename T>
    static T value(const ABC_<T>& args) {
        return BC<J, K>::value * pow_<K>(args.template get<step>()[0])*
                pow_<J - K>(args.template get<step>()[1]) * next_step::template value<PX + K, PY + J - K, PZ>(args);
    }
};
template<size_t step, class next_step, size_t PX, size_t PY, size_t PZ>
struct tetra_integrate_polinomial_sum2 {
    template<size_t I, size_t J, typename T>
    static T value(const ABC_<T>& args) {
        return BC<I, J>::value * pow_<I - J>(args.template get<step>()[2])*
                abstract_sum<J, J, T, tetra_integrate_polinomial_sum3<step, next_step, PX, PY, PZ + I - J>>(args);
    }
};
template<size_t step, class next_step, size_t PX, size_t PY, size_t PZ>
struct tetra_integrate_polinomial_sum1 {
    template<size_t N, size_t I, typename T>
    static T value(const ABC_<T>& args) {
        return BC<N, I>::value * pow_<N - I>(args.template get<step>()[3])*
                abstract_sum<I, I, T, tetra_integrate_polinomial_sum2<step, next_step, PX, PY, PZ>>(args);
    }
};

///==============================================
///
struct tetra_integrate_monomial {
    template<size_t PX, size_t PY, size_t PZ, typename T>
    static T value(const ABC_<T>&) {
        return T(1) / monomial3d<PX, PY, PZ>::value;
    }
};

template<size_t gamma>
struct tetra_integrate_polinomial_step_3 {
    template<size_t PX, size_t PY, size_t PZ, typename T>
    static T value(const ABC_<T>& args) {
        return abstract_sum<gamma, gamma, T, tetra_integrate_polinomial_sum1<3, tetra_integrate_monomial, PX, PY, PZ>>(args);
    }
};

template<size_t beta, size_t gamma>
struct tetra_integrate_polinomial_step_2 {
    template<size_t PX, size_t PY, size_t PZ, typename T>
    static T value(const ABC_<T>& args) {
        return abstract_sum<beta, beta, T, tetra_integrate_polinomial_sum1<2, tetra_integrate_polinomial_step_3<gamma>, PX, PY, PZ>>(args);
    }
};

template<size_t alpha, size_t beta, size_t gamma, typename T>
T tetra_integrate_polinomial_3(const abc_arr<T>& a, const abc_arr<T>& b, const abc_arr<T>& c) {
    return abstract_sum<alpha, alpha, T,
            tetra_integrate_polinomial_sum1<1, tetra_integrate_polinomial_step_2<beta, gamma>, 0, 0, 0>>(ABC_<T>(a, b, c));
}

///Матрица масс

template<typename T>
struct d_vector {
    d_vector(T x, T y, T z) : val_x(x), val_y(y), val_z(z) {}

    T value_x() const {
        return val_x;
    }
    T value_y() const {
        return val_y;
    }
    T value_z() const {
        return val_z;
    }
private:
    T val_x, val_y, val_z;
};

template<size_t i, class Matrix>
struct calc_A {
    typedef typename get_bf_type<i>::type type_i;
    typedef typename Matrix::value_type data_type;
    calc_A(Matrix& A_, data_type J_, const abc_arr<data_type>& a_,
           const abc_arr<data_type>& b_, const abc_arr<data_type>& c_,
           const d_vector<data_type>& dx_) : A(A_), J(J_), a(a_), b(b_), c(c_), dx(dx_)  {}

    template<size_t j>
    void apply() {
        typedef typename get_bf_type<j>::type type_j;
        constexpr size_t alpha = type_i::alpha + type_j::alpha;
        constexpr size_t beta  = type_i::beta + type_j::beta;
        constexpr size_t gamma = type_i::gamma + type_j::gamma;
        A(i, j) = J * tetra_integrate_polinomial_3<alpha, beta, gamma>(a, b, c) /
                (pow_<alpha>(dx.value_x())*pow_<beta>(dx.value_y())*pow_<gamma>(dx.value_z()));
    }
private:
    Matrix& A;
    const data_type J;
    const abc_arr<data_type>& a, b, c;
    const d_vector<data_type>& dx;
};

template<class Matrix>
struct calc_A_i {
    typedef typename Matrix::value_type data_type;
    calc_A_i(Matrix& A_, data_type J_, const abc_arr<data_type>& a_,
           const abc_arr<data_type>& b_, const abc_arr<data_type>& c_,
           const d_vector<data_type>& dx_) : A(A_), J(J_), a(a_), b(b_), c(c_), dx(dx_)  {}

    template<size_t i>
    void apply() {
        calc_A<i, Matrix> closure(A, J, a, b, c, dx);
        meta_loop<i + 1>(closure);
    }
private:
    Matrix& A;
    const data_type J;
    const abc_arr<data_type>& a, b, c;
    const d_vector<data_type>& dx;
};

template<typename T, size_t N, typename = std::enable_if<N == BASE_FUNCTION_COUNT, T>>
void calc_Matrix(matrix_2<N, N, T>& A, T J,
                 const abc_arr<T>& a, const abc_arr<T>& b, const abc_arr<T>& c, const d_vector<T>& dx) {
    calc_A_i<matrix_2<N, N, T>> closure(A, J, a, b, c, dx);
    meta_loop<N>(closure);
}

///Вычисление Якобиана

template<typename T>
T Jacobian(const abc_arr<T>& a, const abc_arr<T>& b, const abc_arr<T>& c) {
    return  a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) + a[2]*(b[0]*c[1] - b[1]*c[0]);
}
template<typename T>
T Jacobian1(const abc_arr<T>& a, const abc_arr<T>& b, const abc_arr<T>& c) {
    return  a[0]*(b[1]*c[2] - b[2]*c[1]) - a[1]*(b[0]*c[2] - b[2]*c[0]) + a[2]*(b[0]*c[1] - b[1]*c[0]);
}

///Вычисление Базисных функций
///


////////////////=============================================================================
struct Grid1 {
    inline Grid1() : a({2.5, 2.3, 1.5, .03}), b({1.5, 1.3, 2.5, 1.3}), c({3.5, .3, .5, 2.4}), dx(.1, .2, .1) {
        centr[0] = (a[0] + a[1] + a[2] + a[3])/4;
        centr[1] = (b[0] + b[1] + b[2] + b[3])/4;
        centr[2] = (c[0] + c[1] + c[2] + c[3])/4;
    }
protected:
    abc_arr<double> a;
    abc_arr<double> b;
    abc_arr<double> c;
    d_vector<double> dx;
    std::array<double, 3> centr;
};

struct phi1_ : Grid1 {
    typedef Grid1 super;
    phi1_(double x_, double y_, double z_) : super() {
        xx = (x_*super::a[0] + y_*super::a[1] + z_*super::a[2] - super::a[3] - super::centr[0])/(super::dx.value_x());
        yy = (x_*super::b[0] + y_*super::b[1] + z_*super::b[2] - super::b[3] - super::centr[1])/(super::dx.value_y());
        zz = (x_*super::c[0] + y_*super::c[1] + z_*super::c[2] - super::c[3] - super::centr[2])/(super::dx.value_z());
    }
protected:
    double xx, yy, zz;
};

///Phi
///
template<size_t N>
struct calc_phi1 : phi1_ {
    typedef phi1_ phi_type;
    calc_phi1(std::array<double, N>& value_phi_, double x_, double y_, double z_) :
                                        value_phi(value_phi_), phi_type(x_, y_, z_) {}

    template<size_t i>
    void apply() {
        typedef typename get_bf_type<i>::type type_i;
        value_phi[i] = pow_<type_i::alpha>(phi_type::xx) * pow_<type_i::beta>(phi_type::yy) * pow_<type_i::gamma>(phi_type::zz);
    }
private:
    std::array<double, BASE_FUNCTION_COUNT>& value_phi;
};
template<size_t N, typename = std::enable_if<N == BASE_FUNCTION_COUNT>>
void calc_phi_f1(std::array<double, N>& value_phi_, double x_, double y_, double z_) {
    calc_phi1<N> closure(value_phi_, x_, y_, z_);
    meta_loop<N>(closure);
}

////Вычислим производные

template<size_t alpha, size_t beta, size_t gamma, typename T = double>
typename std::enable_if_t<(alpha == 0), T> deriv_x(T xx, T yy, T zz, T dx) {
    return T(0);
}
template<size_t alpha, size_t beta, size_t gamma, typename T = double>
typename std::enable_if_t<(alpha != 0), T> deriv_x(T xx, T yy, T zz, T dx) {
    return alpha * pow_<alpha - 1>(xx) / dx * pow_<beta>(yy) * pow_<gamma>(zz);
}

template<size_t alpha, size_t beta, size_t gamma, typename T = double>
typename std::enable_if_t<(beta == 0), T> deriv_y(T xx, T yy, T zz, T dx) {
    return T(0);
}
template<size_t alpha, size_t beta, size_t gamma, typename T = double>
typename std::enable_if_t<(beta != 0), T> deriv_y(T xx, T yy, T zz, T dx) {
    return beta * pow_<alpha>(xx) / dx * pow_<beta - 1>(yy) * pow_<gamma>(zz);
}

template<size_t alpha, size_t beta, size_t gamma, typename T = double>
typename std::enable_if_t<(gamma == 0), T> deriv_z(T xx, T yy, T zz, T dx) {
    return T(0);
}
template<size_t alpha, size_t beta, size_t gamma, typename T = double>
typename std::enable_if_t<(gamma != 0), T> deriv_z(T xx, T yy, T zz, T dx) {
    return gamma * pow_<alpha>(xx) / dx * pow_<beta>(yy) * pow_<gamma - 1>(zz);
}


///Вычисление интегральных операторов

template<class IE>
struct integrate_expression :  math_object_base<IE> {};

struct integrate_data_type {
    double rho = 1., E = 0, u = 0.1, v = 0.1, w = 0.1, e = 0, p = 1.2;
};

struct _rho_t : integrate_expression<_rho_t> { //Начальные данные
    double operator() (const integrate_data_type& id) const {
        return id.rho;
    }
};
extern _rho_t _rho;
struct _E_t : integrate_expression<_E_t> { //Начальные данные
    double operator() (const integrate_data_type& id) const {
        return id.E;
    }
};
extern _E_t _E;
struct _u_t : integrate_expression<_u_t> { //Начальные данные
    double operator() (const integrate_data_type& id) const {
        return id.u;
    }
};
extern _u_t _u;
struct _v_t : integrate_expression<_v_t> { //Начальные данные
    double operator() (const integrate_data_type& id) const {
        return id.v;
    }
};
extern _v_t _v;
struct _w_t : integrate_expression<_w_t> { //Начальные данные
    double operator() (const integrate_data_type& id) const {
        return id.w;
    }
};
extern _w_t _w;
struct _e_t : integrate_expression<_e_t> { //Начальные данные
    double operator() (const integrate_data_type& id) const {
        return id.e;
    }
};
extern _e_t _e;
struct _p_t : integrate_expression<_p_t> { //Начальные данные
    double operator() (const integrate_data_type& id) const {
        return id.p;
    }
};
extern _p_t _p;
////===============================================

template<class IE1, class OP, class IE2>
struct integrate_expression_binary_op : integrate_expression<integrate_expression_binary_op<IE1, OP, IE2> > {
     integrate_expression_binary_op(const integrate_expression<IE1> &ie1, OP op, const integrate_expression<IE2> &ie2) :
         ie1(ie1.self()), op(op), ie2(ie2.self()){}
     double operator()(const integrate_data_type &id) const {
         return op(ie1(id), ie2(id));
     }
private:
     const IE1 ie1;
     OP op;
     const IE2 ie2;
};

template<class IE1, class IE2>
integrate_expression_binary_op<IE1, std::plus<double>, IE2>
operator+(const integrate_expression<IE1> &ie1, const integrate_expression<IE2> &ie2) {
    return integrate_expression_binary_op<IE1, std::plus<double>, IE2>(ie1, std::plus<double>(), ie2);
}
template<class IE1, class IE2>
integrate_expression_binary_op<IE1, std::minus<double>, IE2>
operator-(const integrate_expression<IE1> &ie1, const integrate_expression<IE2> &ie2) {
    return integrate_expression_binary_op<IE1, std::minus<double>, IE2>(ie1, std::minus<double>(), ie2);
}
template<class IE1, class IE2>
integrate_expression_binary_op<IE1, std::multiplies<double>, IE2>
operator*(const integrate_expression<IE1> &ie1, const integrate_expression<IE2> &ie2) {
    return integrate_expression_binary_op<IE1, std::multiplies<double>, IE2>(ie1, std::multiplies<double>(), ie2);
}
template<class IE1, class IE2>
integrate_expression_binary_op<IE1, std::divides<double>, IE2>
operator/(const integrate_expression<IE1> &ie1, const integrate_expression<IE2> &ie2) {
    return integrate_expression_binary_op<IE1, std::divides<double>, IE2>(ie1, std::divides<double>(), ie2);
}

///Abstract summ integrate tetra=============================================================

template<class Closure>
struct abstract_sum_closure {
    typedef typename Closure::value_type value_type;
    abstract_sum_closure(Closure &closure) :  closure(closure), result(value_type()){}

    template<unsigned I>
    void apply(){
        result += closure.template value<I>();
    }
    Closure &closure;
    value_type result;
};

template<unsigned N, class Closure>
typename Closure::value_type abstract_sum1(Closure &closure) {
    abstract_sum_closure<Closure> my_closure(closure);
    meta_loop<N>(my_closure);
    return my_closure.result;
}
///time variable

const double dxs = .01;
static integrate_data_type s_;

template<class IE1, class IE2, class IE3>
struct integrate_tetra_t_base {
    typedef double value_type;

    integrate_tetra_t_base(const integrate_expression<IE1> &ie1_,
                           const integrate_expression<IE2> &ie2_,
                           const integrate_expression<IE3> &ie3_
                           ) : ie1(ie1_.self()), ie2(ie2_.self()), ie3(ie3_.self())  {}
protected:
    const IE1 ie1;
    const IE2 ie2;
    const IE3 ie3;
};


template<size_t J, size_t PX, size_t PY, size_t PZ, class IE1, class IE2, class IE3>
struct tetra_integrate_t_sum3 : integrate_tetra_t_base<IE1, IE2, IE3> {
    typedef integrate_tetra_t_base<IE1, IE2, IE3> super;
    typedef typename super::value_type value_type;
    tetra_integrate_t_sum3(const integrate_expression<IE1> &ie1_,
                           const integrate_expression<IE2> &ie2_,
                           const integrate_expression<IE3> &ie3_
                           ) : super(ie1_, ie2_, ie3_) {}

    template<unsigned K>
    value_type value() const {
        return (super::ie1(s_))*deriv_x<PX + K, PY + J - K, PZ>(X_i[K]/4 + 1/4, Y_i[K]/2 + 1/2, Z_i[K]/2 + 1/2, dxs)
                + (super::ie2(s_))*deriv_y<PX + K, PY + J - K, PZ>(X_i[K]/4 + 1/4, Y_i[K]/2 + 1/2, Z_i[K]/2 + 1/2, dxs)
                + (super::ie3(s_))*deriv_z<PX + K, PY + J - K, PZ>(X_i[K]/4 + 1/4, Y_i[K]/2 + 1/2, Z_i[K]/2 + 1/2, dxs);
    }
};
template<size_t I, size_t PX, size_t PY, size_t PZ, class IE1, class IE2, class IE3>
struct tetra_integrate_t_sum2 : integrate_tetra_t_base<IE1, IE2, IE3> {
    typedef integrate_tetra_t_base<IE1, IE2, IE3> super;
    typedef typename super::value_type value_type;
    tetra_integrate_t_sum2(const integrate_expression<IE1> &ie1_,
                           const integrate_expression<IE2> &ie2_,
                           const integrate_expression<IE3> &ie3_ )
        : super(ie1_, ie2_, ie3_) {}

    template<unsigned J>
    value_type value() const {
        const tetra_integrate_t_sum3<J, PX, PY, PZ + I - J, IE1, IE2, IE3> closure(super::ie1, super::ie2, super::ie3);
        return  abstract_sum1<J + 1>(closure);
    }
private:
    double x, y, z;
};

template<size_t PX, size_t PY, size_t PZ, class IE1, class IE2, class IE3>
struct tetra_integrate_t_sum1 : integrate_tetra_t_base<IE1, IE2, IE3> {
    typedef integrate_tetra_t_base<IE1, IE2, IE3> super;
    typedef typename super::value_type value_type;
    tetra_integrate_t_sum1(const integrate_expression<IE1> &ie1_,
                           const integrate_expression<IE2> &ie2_,
                           const integrate_expression<IE3> &ie3_ )
        : super(ie1_, ie2_, ie3_){}

    template<unsigned I>
    value_type value() const {
        const tetra_integrate_t_sum2<I, PX, PY, PZ, IE1, IE2, IE3> closure(super::ie1, super::ie2, super::ie3);
        return abstract_sum1<I + 1>(closure);
    }
};


template<size_t i, class IE1, class IE2, class IE3, class Matrix>
struct calc_tetra_integrate {
    typedef typename get_bf_type<i>::type type_i;
    typedef typename Matrix::value_type data_type;
    calc_tetra_integrate(const integrate_expression<IE1> &ie1_,
                         const integrate_expression<IE2> &ie2_,
                         const integrate_expression<IE3> &ie3_, Matrix& A_ ) :
                         ie1(ie1_.self()), ie2(ie2_.self()), ie3(ie3_.self()),  A(A_)  {}

    template<size_t j>
    void apply() {
        typedef typename get_bf_type<j>::type type_j;
        constexpr size_t alpha = type_i::alpha + type_j::alpha;
        constexpr size_t beta  = type_i::beta  + type_j::beta;
        constexpr size_t gamma = type_i::gamma + type_j::gamma;
        const tetra_integrate_t_sum1<alpha, beta, gamma, IE1, IE2, IE3> closure(ie1, ie2, ie3);
        A(i, j) = abstract_sum1<j + 1>(closure);
    }
private:
    Matrix& A;
    const IE1 ie1;
    const IE2 ie2;
    const IE3 ie3;
};

template<class IE1, class IE2, class IE3, class Matrix>
struct calc_tetra_integrate_i {
    typedef typename Matrix::value_type data_type;
    calc_tetra_integrate_i(const integrate_expression<IE1> &ie1_,
                           const integrate_expression<IE2> &ie2_,
                           const integrate_expression<IE3> &ie3_, Matrix& A_) :
                           ie1(ie1_.self()), ie2(ie2_.self()), ie3(ie3_.self()), A(A_)  {}

    template<size_t i>
    void apply() {
        calc_tetra_integrate<i, IE1, IE2, IE3, Matrix> closure(ie1, ie2, ie3, A);
        meta_loop<i + 1>(closure);
    }
private:
    Matrix& A;
    const IE1 ie1;
    const IE2 ie2;
    const IE3 ie3;
};

template<class IE1, class IE2, class IE3>
void calc_tetra_integrate_Matrix(const integrate_expression<IE1> &ie1_,
                                 const integrate_expression<IE2> &ie2_,
                                 const integrate_expression<IE3> &ie3_, matrix_2<BASE_FUNCTION_COUNT, BASE_FUNCTION_COUNT, double>& A) {
    calc_tetra_integrate_i<IE1, IE2, IE3, matrix_2<BASE_FUNCTION_COUNT, BASE_FUNCTION_COUNT, double>> closure(ie1_, ie2_, ie3_, A);
    meta_loop<BASE_FUNCTION_COUNT>(closure);
}




////////============================================================================================








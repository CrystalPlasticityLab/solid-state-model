#include "array.h"
#include "vector.h"
#include "quat.h"
#include "matrix.h"
#include "tensor.h"

template <> bool tens::is_not_small_value(long double value) { return (abs(value) > 1e-13) ? true : false;}
template <> bool tens::is_not_small_value(double      value) { return (abs(value) > 1e-13) ? true : false;}
template <> bool tens::is_not_small_value(float       value) { return (abs(value) > 1e-7f) ? true : false;}
template <> bool tens::is_not_small_value(int         value) { return (abs(value) == 0) ? false : true; }

template <> bool tens::is_small_value(long double value) { return (abs(value) < 1e-13) ? true : false;}
template <> bool tens::is_small_value(double      value) { return (abs(value) < 1e-13) ? true : false;}
template <> bool tens::is_small_value(float       value) { return (abs(value) < 1e-7f) ? true : false;}
template <> bool tens::is_small_value(int         value) { return (abs(value) == 0) ? true : false; }

template tens::matrix<double, 3> tens::generate_rand_ort();

template<typename T, std::size_t N>
tens::matrix<T, N> tens::generate_rand_ort() {
	tens::array<T, 4> qvalue(0.0);
	for (size_t row = 0; row < 4; row++) qvalue[row] = static_cast<T>(unidistr(gen));
	quat<T> q(qvalue);
	return static_cast<matrix<T, 3>>(q.get_ort_matrix());
}

template tens::matrix<double, 3> tens::generate_rand_ort();
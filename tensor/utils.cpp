#pragma once
#include "quat.h"

template <> bool tens::is_not_small_value(long double value) { return (abs(value) > 1e-13) ? true : false;}
template <> bool tens::is_not_small_value(double      value) { return (abs(value) > 1e-13) ? true : false;}
template <> bool tens::is_not_small_value(float       value) { return (abs(value) > 1e-7f) ? true : false;}
template <> bool tens::is_not_small_value(int         value) { return (abs(value) == 0) ? false : true; }

template <> bool tens::is_small_value(long double value) { return (abs(value) < 1e-13) ? true : false;}
template <> bool tens::is_small_value(double      value) { return (abs(value) < 1e-13) ? true : false;}
template <> bool tens::is_small_value(float       value) { return (abs(value) < 1e-7f) ? true : false;}
template <> bool tens::is_small_value(int         value) { return (abs(value) == 0) ? true : false; }


tens::container<double, 3, 2> tens::generate_rand_ort() {
	container<double, 4, 1> arr_rand;
	arr_rand.fill_rand();
	quat<double> q(arr_rand);
	return get_ort_matrix<double>(q);
}

tens::container<double, 3, 2> tens::generate_indent_ort() {
	quat<double> q;
	return get_ort_matrix<double>(q);
}
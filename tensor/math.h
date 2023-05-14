#pragma once
#include "error.h"

namespace math {

	template<typename T> bool is_not_small_value(T value);
	template<typename T> bool is_small_value(T value);

	namespace dim3 {
		template<typename T>
		inline void mat_scal_mat_transp(const T* lhs, const T* rhs, T* nhs) {
			nhs[0] = lhs[0] * rhs[0] + lhs[4] * rhs[4] + lhs[5] * rhs[5];
			nhs[1] = lhs[1] * rhs[1] + lhs[3] * rhs[3] + lhs[8] * rhs[8];
			nhs[2] = lhs[2] * rhs[2] + lhs[6] * rhs[6] + lhs[7] * rhs[7];
			nhs[3] = lhs[3] * rhs[2] + lhs[1] * rhs[6] + lhs[8] * rhs[7];
			nhs[4] = lhs[4] * rhs[2] + lhs[5] * rhs[6] + lhs[0] * rhs[7];
			nhs[5] = lhs[5] * rhs[1] + lhs[4] * rhs[3] + lhs[0] * rhs[8];
			nhs[6] = lhs[6] * rhs[1] + lhs[2] * rhs[3] + lhs[7] * rhs[8];
			nhs[7] = lhs[7] * rhs[0] + lhs[2] * rhs[4] + lhs[6] * rhs[5];
			nhs[8] = lhs[8] * rhs[0] + lhs[3] * rhs[4] + lhs[1] * rhs[5];
		}

		template<typename T>
		inline void mat_scal_mat(const T* lhs, const T* rhs, T* nhs) {
			nhs[0] = lhs[0] * rhs[0] + lhs[4] * rhs[7] + lhs[5] * rhs[8];
			nhs[1] = lhs[1] * rhs[1] + lhs[8] * rhs[5] + lhs[3] * rhs[6];
			nhs[2] = lhs[2] * rhs[2] + lhs[6] * rhs[3] + lhs[7] * rhs[4];
			nhs[3] = lhs[3] * rhs[2] + lhs[1] * rhs[3] + lhs[8] * rhs[4];
			nhs[4] = lhs[4] * rhs[2] + lhs[5] * rhs[3] + lhs[0] * rhs[4];
			nhs[5] = lhs[5] * rhs[1] + lhs[0] * rhs[5] + lhs[4] * rhs[6];
			nhs[6] = lhs[6] * rhs[1] + lhs[7] * rhs[5] + lhs[2] * rhs[6];
			nhs[7] = lhs[7] * rhs[0] + lhs[2] * rhs[7] + lhs[6] * rhs[8];
			nhs[8] = lhs[8] * rhs[0] + lhs[3] * rhs[7] + lhs[1] * rhs[8];
		}

		template<typename T>
		inline T mat_conv_transp(const T* lhs, const T* rhs) {
			T res(0);
			res += lhs[0] * rhs[0];
			res += lhs[1] * rhs[1];
			res += lhs[2] * rhs[2];
			res += lhs[3] * rhs[6];
			res += lhs[4] * rhs[7];
			res += lhs[5] * rhs[8];
			res += lhs[6] * rhs[3];
			res += lhs[7] * rhs[4];
			res += lhs[8] * rhs[5];
			return res;
		}

		template<typename T>
		inline void mat_scal_vect(const T* m, const T* a, T* nhs) {
			nhs[0] = a[0] * m[0] + a[2] * m[4] + a[1] * m[5];
			nhs[1] = a[1] * m[1] + a[2] * m[3] + a[0] * m[8];
			nhs[2] = a[2] * m[2] + a[1] * m[6] + a[0] * m[7];
		}
		
		template<typename T>
		inline void vect_scal_mat(const T* a, const T* m, T* nhs) {
			nhs[0] = a[0] * m[0] + a[2] * m[7] + a[1] * m[8];
			nhs[1] = a[1] * m[1] + a[0] * m[5] + a[2] * m[6];
			nhs[2] = a[2] * m[2] + a[1] * m[3] + a[0] * m[4];
		}

		template<typename T>
		inline T vect_scal_vect(const T* a, const T* b) {
			return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
		}

		template<typename T>
		inline T det_mat(const T* m){
			return m[0] * m[1] * m[2] - m[0] * m[3] * m[6] - m[1] * m[4] * m[7] + m[3] * m[5] * m[7] - m[2] * m[5] * m[8] + m[4] * m[6] * m[8];
		}

		template<typename T>
		inline void inv_mat(const T* m, T* inv_matr) {
			T div = det_mat(m);
			if (is_small_value(div)) {
				throw ErrorMath::DivisionByZero();
			}
			div = T(1) / div;
			inv_matr[0] = (m[1] * m[2] - m[3] * m[6])*div; inv_matr[5] = (m[4] * m[6] - m[2] * m[5])*div; inv_matr[4] = (m[3] * m[5] - m[1] * m[4])*div;
			inv_matr[8] = (m[3] * m[7] - m[2] * m[8])*div; inv_matr[1] = (m[0] * m[2] - m[4] * m[7])*div; inv_matr[3] = (m[4] * m[8] - m[0] * m[3])*div;
			inv_matr[7] = (m[6] * m[8] - m[1] * m[7])*div; inv_matr[6] = (m[5] * m[7] - m[0] * m[6])*div; inv_matr[2] = (m[0] * m[1] - m[5] * m[8])*div;
		}
	}
}

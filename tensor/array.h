#pragma once
#include "container.h"
#include "matrix.h"

namespace tens {



	template<typename T, std::size_t N> class array;
	template<typename T, std::size_t N> std::ostream& operator<< (std::ostream& o, const array<T, N>& v);

	template<typename T, size_t N>
	class array : public container<T, N>
	{
		void _copy(const std::array<T,N>& arr){
			for (size_t i = 0; i < N; i++){
				(*this)[i] =  arr[i];
			}
		}
	protected:
		array() : container<T, N>(1) {
			container<T, N>::fill_value(T(0));
		};
	public:
		~array() { };
		explicit array(ARRAY_TYPE AT);
		explicit array(const T& val) : container<T, N>(1) { container<T, N>::fill_value(T(0)); };
		explicit array(const std::array<T, N>& arr) : container<T, N>(1) { this->_copy(arr); };
		explicit array(const array& a) : container<T, N>(static_cast<const container<T, N>&>(a)) {};
		array(array&& a)noexcept : container<T, N>(static_cast<container<T, N>&&>(a)) {};

		static array generate_rand() {
			array<T, N> a;
			a.fill_rand();
			return a;
		};

		friend std::ostream& operator<< <>(std::ostream& out, const array& a);

		inline array& operator= (const array& rhs)	{ // copy assign operator
			container<T, N>::operator = (static_cast<const container<T, N>&>(rhs));
			return *this;
		}

		inline array& operator= (array&& rhs) noexcept {  // move assign operator
			container<T, N>::operator = (static_cast<container<T, N>&&>(rhs));
			return *this;
		}
		
		inline array  operator- () const;
		inline array  operator+ (const array& rhs) const;
		inline array  operator- (const array& rhs) const;
		inline     T  operator* (const array& v) const;
		inline array  operator* (const T& mult) const;
		inline array  operator/ (const T& div) const;
		inline array& operator= (const T& val);
		inline bool   operator==(const array& val) const;
		inline array& operator/=(const T& div);
		inline array& operator*=(const T& mult);
		inline array& operator+=(const array& rhs);
		inline array& operator-=(const array& rhs);

		inline matrix<T, N>  outer_product(const array& rhs) const;

		friend array<T, N> array_vector_product(const array<T, N>& lhs, const array<T, N>& rhs) {
			static_assert(N == 3 && " vector_product is implemented only for 3-dimensional vectors.");
			array<T, N> comp;
			comp[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
			comp[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
			comp[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
			return comp;
		}

		virtual inline T   norm() const;
		inline array   get_normalize() const;
		inline virtual void normalize();
		inline static  void normalize(array&);

		friend array<T, N > operator* (const array<T, N >& v, const matrix<T, N >& m) {
			array <T, N> res((T)0);
			for (size_t row = 0; row < N; row++)
				for (size_t col = 0; col < N; col++)
					res[row] += v[col] * m[col][row];

			return res;
		}

		friend array<T, N > operator* (const matrix<T, N >& m, const array<T, N >& v) {
			array <T, N> res((T)0);
			for (size_t row = 0; row < N; row++)
				for (size_t col = 0; col < N; col++)
					res[row] += v[col] * m[row][col];

			return res;
		}
	};

	template<typename T, std::size_t N>
	array<T, N>::array(ARRAY_TYPE AT) :container<T,N>(1) {
		switch (AT)
		{
		case ARRAY_TYPE::ZERO:		
			container<T, N>::fill_value((T)0);
			return;
		case ARRAY_TYPE::RANDOM: 
			container<T, N>::fill_rand();
			return;
		case ARRAY_TYPE::RANDOMUNIT:
			container<T, N>::fill_rand();
			normalize();
			return;
		default:                     	
			container<T, N>::fill_value((T)0);
			return;
		}
		return;
	}

	template<typename T, size_t N>
	inline array<T, N> array<T,N>::operator- () const {
		array<T, N> res;
		for (size_t row = 0; row < N; row++)
			res[row] = -(*this)[row];
		return res;
	};

	template<typename T, size_t N>
	inline array<T, N> array<T, N>::operator+ (const array& rhs) const {
		array<T, N> res;
		for (size_t row = 0; row < N; row++)
			res[row] = (*this)[row] + rhs[row];
		return res;
	};

	template<typename T, size_t N>
	inline array<T, N> array<T, N>::operator- (const array& rhs) const {
		array<T, N> res;
		for (size_t row = 0; row < N; row++)
			res[row] = (*this)[row] - rhs[row];
		return res;
	};

	template<typename T, size_t N>
	inline array<T, N> array<T, N>::operator * (const  T& mult) const {
		array<T, N> res;
		for (size_t row = 0; row < N; row++)
			res[row] = (*this)[row] * mult;
		return res;
	};

	template<typename T, size_t N>
	inline array<T, N> array<T, N>::operator / (const  T& div) const {
		array<T, N> res;
		if (is_not_small_value(div)) {
			T mult = (T)1 / div;
			for (size_t row = 0; row < N; row++)
				res[row] = (*this)[row] * mult;
		}
		return res;
	};

	template<typename T, size_t N>
	inline array<T, N>& array<T, N>::operator /= (const  T& div) {
		array<T, N> res;
		if (abs(div) > (T)0) {
			T mult = (T)1 / div;
			for (size_t row = 0; row < N; row++)
				(*this)[row] *= mult;
		}
		return *this;
	};

	template<typename T, size_t N>
	inline array<T, N>& array<T, N>::operator = (const T& val) {
		fill(val);
		return (*this);
	}

	template<typename T, size_t N>
	inline bool array<T, N>::operator == (const array<T, N>& val) const {
		const array<T, N> diff = val - *this;
		if (is_small_value(diff.norm())) {
			return true;
		}
		return false;
	}

	template<typename T, size_t N>
	inline array<T, N>& array<T, N>::operator *=(const T& mult) {
		for (size_t row = 0; row < N; row++)
			(*this)[row] *= mult;
		return (*this);
	};

	template<typename T, size_t N>
	inline array<T, N>& array<T, N>::operator +=(const array& rhs) {
		for (size_t row = 0; row < N; row++)
			(*this)[row] += rhs[row];
		return *this;
	}

	template<typename T, size_t N>
	inline array<T, N>& array<T, N>::operator -=(const array& rhs) {
		for (size_t row = 0; row < N; row++)
			(*this)[row] -= rhs[row];
		return *this;
	}

	template<typename T, size_t N>
	std::ostream& operator<<(std::ostream& out, const array<T,N>& a){
		for (size_t row = 0; row < N; row++)
			out << a[row] << " ";
		out << "\n";
		return out;
	}

	template<typename T, size_t N>
	inline matrix<T, N> array<T, N >::outer_product(const array<T, N >& rhs) const {
		const array<T, N>& lhs = *this;
		std::array< std::array<T, N>,N> res;
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				res[row][col] = lhs[row] * rhs[col];
		return static_cast<matrix<T, N>> (res);
	}

	template<typename T, size_t N>
	inline T array<T, N >::operator* (const array<T, N >& v) const {
		T res = (T)0;
		for (size_t row = 0; row < N; row++)
			res += (*this)[row] * v[row];

		return res;
	}

	template<typename T, size_t N>
	inline T    array<T, N >::norm() const {
		const array<T, N >& v = (*this);
		return sqrt(v * v);
	}

	template<typename T, size_t N>
	inline void array<T, N>::normalize(array<T, N>& v) {
		T length = v.norm();
		if (is_small_value(length)){
			throw ErrorMath::DivisionByZero();
		}
		v *= (T)1 / length;
	}


	template<typename T, size_t N>
	inline array<T, N > array<T, N>::get_normalize() const {
		array<T, N> res(*this);
		normalize(res);
		return res;
	}

	template<typename T, size_t N>
	inline void array<T, N >::normalize() {
		normalize(*this);
	}
};
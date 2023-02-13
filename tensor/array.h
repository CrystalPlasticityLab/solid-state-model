#pragma once
#include "container.h"
#include "matrix.h"

namespace tens {

	enum class ARRAYTTYPE{
		ZERO,
		RANDOM,
		RANDOMUNIT
	};

	template<typename T, std::size_t N> class array;
	template<typename T, std::size_t N> std::ostream& operator<< (std::ostream& o, const array<T, N>& v);

	template<typename T, size_t N>
	class array : public std::unique_ptr<T[]>
	{
		size_t  _dim = 0;
		void _fill(const T& val){
			for (size_t i = 0; i < N; i++){
				(*this)[i] =  val;
			}
		}
		void _fill_random(){
			for (size_t i = 0; i < N; i++) {
				(*this)[i] = static_cast<T>(unidistr(gen));
			}
		}
		void _copy(const T* const values){
			for (size_t i = 0; i < N; i++){
				(*this)[i] =  values[i];
			}
		}
		void _copy(const std::array<T,N>& arr){
			for (size_t i = 0; i < N; i++){
				(*this)[i] =  arr[i];
			}
		}
		void _copy(const array& arr){
			for (size_t i = 0; i < N; i++){
				(*this)[i] =  arr[i];
			}
		}
	protected:
		array() : std::unique_ptr<T[]>(new T[N]), _dim(N) { 
			this->_fill(T());
		};
	public:
		~array() { };
		explicit array(ARRAYTTYPE AT);
		explicit array(const T& val)               : std::unique_ptr<T[]>(new T[N]), _dim(N) { this->_fill(val); }; // public ctor by const value
		explicit array(const std::array<T,N>& arr) : std::unique_ptr<T[]>(new T[N]), _dim(N) { this->_copy(arr); }; // public ctor by astd::array
		explicit array(const array&  a)            : std::unique_ptr<T[]>(new T[a._dim]), _dim(a._dim)  { this->_copy(a.get()); }; // copy ctor
		array(array&& v)noexcept : std::unique_ptr<T[]>(std::move(v))         { this->_dim = v._dim; }; // move ctor

		friend std::ostream& operator<< <>(std::ostream& out, const array& a);

		inline array& operator= (const array& rhs)	{ // copy assign operator
			_copy(rhs);
			this->_dim = rhs._dim;
			return *this;
		}

		inline array& operator= (array&& rhs) noexcept {  // move assign operator
			this->_dim = rhs._dim;
			static_cast<std::unique_ptr<T[]>&>(*this) = std::move(rhs);
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
	array<T, N>::array(ARRAYTTYPE AT) : std::unique_ptr<T[]>(new T[N]), _dim(N) {
		switch (AT)
		{
		case ARRAYTTYPE::ZERO:		
			this->_fill((T)0);
			return;
		case ARRAYTTYPE::RANDOM: 
			_fill_random();
			return;
		case ARRAYTTYPE::RANDOMUNIT:
			_fill_random();
			normalize();
			return;
		default:                     	
			this->_fill((T)0);
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
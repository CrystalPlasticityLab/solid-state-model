#pragma once
#include <random>
#include <algorithm>
#include "container.h"

extern const size_t DIM;
extern std::random_device rd;  // Will be used to obtain a seed for the random number engine
extern std::mt19937 gen;       // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> unidistr;

namespace tens {
	enum class MATRIXINITTYPE{
		ZERO,
		INDENT,
	};

	enum class TRANSPOSE
	{
		TRUE,
		FALSE
	};

	template<typename T, std::size_t N> class matrix;
	template<typename T, std::size_t N> std::ostream& operator<< (std::ostream& o, const matrix<T,N>& m);
	template<typename T, std::size_t N> matrix<T, N> generate_rand_ort();
	template<typename T, std::size_t N> matrix<T, N> transpose();

	template<typename T, std::size_t N>
	class matrix : public std::unique_ptr<T[]>
	{
		typedef  std::array<std::array<T, N>, N>  array;

		size_t  _dim;
		std::array<T*, N> _rows;
		void set_zero();
		
		void _copy(const std::array<std::array<T, N>, N>& arr){
			this->_dim = N;
			_set_rows();
			for (size_t row = 0; row < N; row++){
				for (size_t col = 0; col < N; col++){
					_rows[row][col] =  arr[row][col];
				}
			}
		}
		
		void _copy(const matrix& m){
			this->_dim = m._dim;
			memcpy(this->get(), m.get(), m._dim*m._dim*sizeof(T));
			_set_rows();
		}

		void _move(matrix& _m){
			matrix&& m = std::move(_m);
			this->_dim = m._dim;
			static_cast<std::unique_ptr<T[]>&>(*this) = std::move(m);
			std::swap(m._rows, this->_rows);
		}

		void _set_rows(){
			const T* ptr = this->get();
			for (size_t row = 0; row < _dim; row++){
				_rows[row] = ptr + row*_dim;
			}
		}

		matrix() : std::unique_ptr<T[]>(new T[N*N]), _dim(N) { _set_rows(); }; // default ctor
	protected:
	public:
		~matrix() { };
		explicit matrix(MATRIXINITTYPE IT) ;
		explicit matrix(const array& arr)  : std::unique_ptr<T[]>(new T[N*N]), _dim(N)  { _copy(arr); };
		matrix(const matrix& m)   : std::unique_ptr<T[]>(new T[m._dim*m._dim]) { _copy(m); }; // copy ctor
		matrix(matrix&& m)noexcept         : std::unique_ptr<T[]>(), _dim(0)   { _move(m);};  // move ctor

		inline T* operator [](size_t i)   { return _rows[i]; };
		
		inline matrix& operator= (const matrix& rhs) {
			this->_copy(rhs); 
			return *this;
		}

		inline matrix& operator= (matrix&& rhs) noexcept {
			this->_move(rhs);
			return *this;
		}

		inline  matrix  operator- () const;
		inline  bool    operator==(const matrix& m) const;
		inline  matrix  operator+ (const matrix& m) const;
		inline  matrix  operator- (const matrix& m) const;
		inline  matrix  operator* (const matrix& m) const;
		inline  matrix  operator* (const T& val) const;
		inline  matrix& operator+=(const matrix& m);
		inline  matrix& operator-=(const matrix& m);
		inline  matrix& operator*=(const matrix& m);
		inline  matrix& operator*=(const T& val);
		inline  matrix& operator/=(const T& val);

		virtual inline T   norm() const;
		bool        check_ort() const;
		inline matrix transpose();

		friend void transpose(matrix& m) {
			for (size_t row = 0; row < N; row++)
				for (size_t col = row + 1; col < N; col++)
					std::swap(m[row][col], m[col][row]);
		}

		inline matrix scal(TRANSPOSE left, const matrix& rhs, TRANSPOSE right) const;
		inline matrix transform(TRANSPOSE left, const matrix& op, TRANSPOSE right) const;
		inline  T   convolution(const matrix<T, N>& rhs) const;
		
		friend std::ostream& operator<< <>(std::ostream& out, const matrix& a);
	};

	template<typename T, std::size_t N>
	matrix<T, N> generate_rand_ort();

	template<typename T, std::size_t N>
	matrix<T, N> generate_rand() {
		std::array<std::array<T, N>, N> a;
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				a[row][col] = static_cast<T>(unidistr(gen));
		return matrix<T, N>(a);
	};

	
	template<typename T, std::size_t N>
	void matrix<T, N>::set_zero() {
		std::fill(this->get(), this->get() +  _dim*_dim, (T)0);
	}

	template<typename T, std::size_t N>
	matrix<T, N>::matrix(MATRIXINITTYPE IT) : std::unique_ptr<T[]>(new T[N*N]), _dim(N) {
		_set_rows();
		switch (IT)
		{
		case MATRIXINITTYPE::ZERO  :		
			set_zero();
			return;
		case MATRIXINITTYPE::INDENT:
			set_zero();
			for (size_t row = 0; row < N; row++) 
				(*this)[row][row] = (T)1; 
			return;
		default:                     		
			set_zero();
			return;
		}
		return;
	}

	template<typename T, std::size_t N>
	matrix<T, N> matrix<T, N>::transpose() {
		matrix<T, N> res(*this);
		for (size_t row = 0; row < N; row++)
			for (size_t col = row + 1; col < N; col++)
				std::swap(res[row][col], res[col][row]);
		return res;
	}


	template<typename T, size_t N>
	inline T    matrix<T, N >::norm() const {
		return sqrt(this->convolution(*this));
	}

	template<typename T, size_t N>
	bool matrix<T, N >::check_ort() const {
		matrix<T, N> m = *this * this->transpose();
		T diag = 0;
		T nondiag = 0;
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				(row == col) ? diag += m[row][row] : nondiag += m[row][col];
		return (is_small_value<T>(abs(diag - (T)N) + abs(nondiag)) ? true : false);
	}

	// this^?l . rhs^?r
	template<typename T, std::size_t N>
	matrix<T, N> matrix<T, N>::scal(TRANSPOSE left, const matrix& rhs, TRANSPOSE right) const {
		if (left == right)	{
			if (left == TRANSPOSE::FALSE)	return (*this) * rhs;
			else		                	return (rhs * (*this)).transpose();
		} else	{
			if (left == TRANSPOSE::FALSE)	return (*this) * rhs.transpose();
			else                            return transpose() * rhs;
		}
	}

	//  op^?l . this . op^?r
	template<typename T, std::size_t N>
	matrix<T, N> matrix<T, N>::transform(TRANSPOSE left, const matrix& op, TRANSPOSE right) const {
		matrix<T, N> opt = op.transpose();
		if (left == right)	{
			if (left == TRANSPOSE::FALSE)	return op * (*this) * op;
			else			                return opt * (*this) * opt;
		} else	{
			if (left == TRANSPOSE::FALSE)	return op * (*this) * opt;
			else                			return opt * (*this) * op;
		}
	}
	template<typename T, size_t N>
	inline matrix<T, N> matrix<T,N>::operator- () const {
		matrix<T, N> res;
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				res[row][col] = -(*this)[row][col];
		return res;
	};
 
	template<typename T, std::size_t N>
	bool matrix<T, N>::operator==(const matrix<T, N>& rhs) const {
		const matrix<T, N> diff = rhs - *this;
		if (is_small_value(diff.convolution(diff))) {
			return true;
		}
		return false;
	}

	template<typename T, std::size_t N>
	matrix<T, N> matrix<T, N>::operator + (const matrix<T, N>& rhs) const {
		matrix<T, N> nhs;
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				nhs[row][col] = (*this)[row][col] + rhs[row][col];
		return nhs;
	}

	template<typename T, std::size_t N>
	matrix<T, N> matrix<T, N>::operator - (const matrix<T, N>& rhs) const {
		matrix<T, N> nhs;
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				nhs[row][col] = (*this)[row][col] - rhs[row][col];
		return nhs;
	}

	template<typename T, std::size_t N>
	matrix<T, N>& matrix<T, N>::operator += (const matrix<T, N>& rhs) {
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				(*this)[row][col] += rhs[row][col];
		return *this;
	}

	template<typename T, std::size_t N>
	matrix<T, N>& matrix<T, N>::operator -= (const matrix<T, N>& rhs) {
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				(*this)[row][col] -= rhs[row][col];
		return *this;
	}

	template<typename T, std::size_t N>
	matrix<T, N> matrix<T, N>::operator * (const T& val) const {
		matrix<T, N> nhs;
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				nhs[row][col] = (*this)[row][col] * val;
		return nhs;
	}

	template<typename T, std::size_t N>
	matrix<T, N> matrix<T, N>::operator * (const matrix<T, N>& rhs) const {
		matrix<T, N> nhs(MATRIXINITTYPE::ZERO);
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				for (size_t i = 0; i < N; i++)
					nhs[row][col] += (*this)[row][i] * rhs[i][col];
		return nhs;
	}

	template<typename T, std::size_t N>
	T   matrix<T, N>::convolution(const matrix<T, N>& rhs) const {
		T res = (T)0;
		const matrix<T, N>& lhs = *this;
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				res += lhs[row][col] * rhs[col][row];
		return res;
	}

	template<typename T, std::size_t N>
	matrix<T, N>& matrix<T, N>::operator *=(const T& val)  {
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				(*this)[row][col] *= val;
		return *this;
	}

	template<typename T, std::size_t N>
	matrix<T, N>& matrix<T, N>::operator /=(const T& val)  {
		const double inv = 1.0/val;
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				(*this)[row][col] *= inv;
		return *this;
	}

	template<typename T, std::size_t N>
	matrix<T, N>& matrix<T, N>::operator *= (const matrix<T, N>& rhs) {
		matrix<T, N> nhs(MATRIXINITTYPE::ZERO);
		for (size_t row = 0; row < N; row++)
			for (size_t col = 0; col < N; col++)
				for (size_t i = 0; i < N; i++)
					nhs[row][col] += (*this)[row][i] * rhs[i][col];

		(*this) = nhs;
		return *this;
	}

	template<typename T, std::size_t N>
	std::ostream& operator<<(std::ostream& out, const matrix<T, N>& a) {
		for (size_t row = 0; row < N; row++){
			for (size_t col = 0; col < N; col++){
					out << a[row][col] << " ";
			}
				out << "\n";
		}
		out << "\n";
		return out;
	};
}
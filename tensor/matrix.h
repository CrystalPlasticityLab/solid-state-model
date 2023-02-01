#pragma once
#include <random>
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
	class matrix : public container <std::array<std::array<T, N>, N>>
	{
		typedef  std::array<std::array<T, N>, N>  _array;
		typedef  container <_array> _cont;

		void set_zero();
		matrix() : _cont() {};
	protected:
	public:
		~matrix() { };
		explicit matrix(const _array& _arr)       : _cont(_arr) { };
		explicit matrix(const matrix& m)     : _cont(static_cast<const _cont&>(m)) { }; // copy constructor
		matrix(matrix&& m)noexcept : _cont(static_cast<_cont&&>(m)) { };      // move constructor
		explicit matrix(MATRIXINITTYPE IT) ;

		inline       std::array<T, N>& operator [](size_t i)       { return (*this->_Elem)[i]; };
		inline const std::array<T, N>& operator [](size_t i) const { return (*this->_Elem)[i]; };

		friend std::ostream& operator<< <>(std::ostream& out, const matrix& a);
		inline matrix& operator= (const matrix& rhs) {
			container<_array>::operator = (static_cast<const container<_array>&>(rhs));
			return *this;
		}
		inline matrix& operator= (matrix&& rhs) noexcept {
			container<_array>::operator = (static_cast<container<_array>&&>(rhs));
			return *this;
		}

		inline  matrix  operator- () const;
		inline  bool         operator==(const matrix& m) const;
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
		inline matrix transpose() const;
		inline matrix scal(TRANSPOSE left, const matrix& rhs, TRANSPOSE right) const;
		inline matrix transform(TRANSPOSE left, const matrix& op, TRANSPOSE right) const;
		inline  T   convolution(const matrix<T, N>& rhs) const;
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
		for (auto& row : *this->_Elem) row.fill((T)0);
	}

	template<typename T, std::size_t N>
	matrix<T, N>::matrix(MATRIXINITTYPE IT) : _cont() {

		switch (IT)
		{
		case MATRIXINITTYPE::ZERO  :		set_zero();
			return;
		case MATRIXINITTYPE::INDENT:
			set_zero();
			for (size_t row = 0; row < N; row++) 
				(*this)[row][row] = (T)1; 
			return;
		default:                     		set_zero();
			return;
		}
		return;
	}

	template<typename T, std::size_t N>
	matrix<T, N> matrix<T, N>::transpose() const {

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
		if ( is_small_value(diff.convolution(diff))) {
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
		return (*this);
	}

	template<typename T, std::size_t N>
	std::ostream& operator<<(std::ostream& out, const matrix<T, N>& a) {
		for (const auto& row : *a._Elem)
		{
			for (const auto& col : row)
				out << col << " ";
			out << "\n";
		}
		out << "\n";
		return out;
	};
}
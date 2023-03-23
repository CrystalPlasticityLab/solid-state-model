#pragma once
#include <memory>
#include <array>
#include <random>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <utility>
#include "error.h"
#include "math.h"

extern std::mt19937 gen;       // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> unidistr;

namespace tens {

	enum class FILL_TYPE {
		ZERO,
		RANDOM,
		RANDOMUNIT,
		INDENT
	};

	template<typename T>
	class container;
	template<typename T> std::ostream& operator<< (std::ostream& o, const container<T>& cont);

	template<typename T, size_t N>
	container<T> Matrix(const std::array<std::array<T, N>, N>& matrix) {
		std::unique_ptr<T[]> arr = std::unique_ptr<T[]>(new T[N*N]);
		if (N == 3) {// ------- for DIM == 3 : {00, 11, 22, 12, 02, 01, 21, 20, 10} 
			arr[0] = matrix[0][0];
			arr[1] = matrix[1][1];
			arr[2] = matrix[2][2];
			arr[3] = matrix[1][2];
			arr[4] = matrix[0][2];
			arr[5] = matrix[0][1];
			arr[6] = matrix[2][1];
			arr[7] = matrix[2][0];
			arr[8] = matrix[1][0];
			return container<T>(3, 2, std::move(arr));
		}
		throw new NoImplemetationYet();
	}
	
	template<typename T, size_t N>
		container<T> Matrix(FILL_TYPE type) {
		return container<T>(N, 2, type);
	}

	template<typename T, size_t N>
	container<T> Array(FILL_TYPE type) {
		return container<T>(N, 1, type);
	}

	template<typename T, size_t N>
	container<T> Array(const std::array<T, N>& array) {
		std::unique_ptr<T[]> arr = std::unique_ptr<T[]>(new T[N]);
		memcpy(arr.get(), &array, N * sizeof(T));
		return container<T>(N, 1, std::move(arr));
	}

	container<double> generate_rand_ort();
	container<double> generate_indent_ort();

	template<typename T>
	class container : public std::unique_ptr<T[]>
	{
	private:
		size_t _dim;
		size_t _size;
		size_t _rank;

		void _copy(const T* ptr) {
			memcpy(this->get(), ptr, _size * sizeof(T));
		}

		void _copy(const container& c) {
			this->_dim = c._dim;
			this->_size = c._size;
			this->_rank = c._rank;
			memcpy(this->get(), c.get(), _size * sizeof(T));
		}

		template<size_t N>
		void _copy(const std::array<T, N>& arr) {
			this->_dim = N;
			this->_size = N;
			this->_rank = 1;
			memcpy(this->get(), &arr, _size * sizeof(T));
		}

		void _move(container& c) {
			this->_dim  = c._dim;   c._dim = 0;
			this->_size = c._size;  c._size = 0;
			this->_rank = c._rank;  c._rank = 0;
			static_cast<std::unique_ptr<T[]>&>(*this) = std::move(c);
		}

		void _alloc() {
			if (*this) {
				throw new NoImplemetationYet();
			}
			static_cast<std::unique_ptr<T[]>&>(*this) = std::unique_ptr<T[]>(new T[_size]);
		}

		void _free() {
			this->reset();
		}

	public:
		size_t dim()  const { return _dim; };
		size_t size() const { return _size; };
		size_t rank() const { return _rank; };

		void is_consist(const container& c) const {
			if (this->_dim != c._dim || this->_rank != c._rank) {
				throw new ErrorMath::ShapeMismatch();
			};
		};

		void fill_rand() {
			for (size_t i = 0; i < _size; i++) {
				this->get()[i] = static_cast<T>(unidistr(gen));
			}
		};

		void fill_value(const T& val) {
			std::fill(this->get(), this->get() + _size, val);
		};

		container(size_t N, size_t R, std::unique_ptr<T[]>&& pointer) : std::unique_ptr<T[]>(std::move(pointer)), _dim(N), _size(R != 0 ? (size_t)pow(N, R) : N), _rank(R) {};

		container(size_t N, size_t R, FILL_TYPE type) : std::unique_ptr<T[]>(), _dim(N), _size(R != 0 ? (size_t)pow(N, R) : N), _rank(R) {
			_alloc();
			switch (type)
			{
			case tens::FILL_TYPE::ZERO:
				fill_value(T(0));
				break;
			case tens::FILL_TYPE::RANDOM:
				fill_rand();
				break;
			case tens::FILL_TYPE::RANDOMUNIT:
				fill_rand();
				*this = get_normalize(*this);
				break;
			case tens::FILL_TYPE::INDENT:
				fill_value(T(0));
				for (size_t i = 0; i < N; i++){
					this->get()[i] = T(1);
				}
				break;
			default:
				break;
			}
		};
		
		container(size_t N, size_t R) : std::unique_ptr<T[]>(), _dim(N), _size( R != 0 ? (size_t)pow(N, R) : N), _rank(R) {
			_alloc();
		};

		container(size_t N, size_t R, const T& val) : std::unique_ptr<T[]>(), _dim(N), _size(R != 0 ? (size_t)pow(N, R) : N), _rank(R) {
			_alloc();
			fill_value(val);
		};

		container(const container& c) : std::unique_ptr<T[]>(), _dim(c._dim), _size(c._size), _rank(c._rank) {
			_alloc();
			_copy(c);
		};

		container(container&& c) noexcept : std::unique_ptr<T[]>(std::move(c)), _dim(c._dim), _size(c._size), _rank(c._rank) {
		};

		operator T() const { 
			if (_size == 1) {
				return (*this)[0];
			}
			throw ErrorAccess::NoCastScalar();
		}

		inline container& operator= (const T value) {
			fill_value(value);
			return *this;
		}

		inline container& operator= (const container& rhs) {
#ifdef _DEBUG
			this->is_consist(rhs);
#endif
			this->_copy(rhs);
			return *this;
		}

		inline container& operator= (container&& rhs) noexcept {
#ifdef _DEBUG
			this->is_consist(rhs);
#endif
			this->_move(rhs);
			return *this;
		}

		[[nodiscard]] inline friend container<T> operator + (const container<T>& lhs, const container<T>& rhs) {
#ifdef _DEBUG
			lhs.is_consist(rhs);
#endif
			container<T> nhs(lhs.dim(), rhs.rank());
			for (size_t i = 0; i < lhs.size(); i++)
				nhs[i] = lhs[i] + rhs[i];
			return nhs;
		}

		[[nodiscard]] inline friend container<T> operator - (const container<T>& lhs, const container<T>& rhs) {
#ifdef _DEBUG
			lhs.is_consist(rhs);
#endif
			container nhs(lhs.dim(), rhs.rank());
			for (size_t i = 0; i < lhs.size(); i++)
				nhs[i] = lhs[i] - rhs[i];
			return nhs;
		}

		[[nodiscard]] inline friend container<T> operator * (const container<T>& lhs, const T& mul) {
			container nhs(lhs);
			for (size_t i = 0; i < lhs.size(); i++)
				nhs[i] *= mul;
			return nhs;
		}

		[[nodiscard]] inline friend container<T> operator / (const container<T>& lhs, const T& div) {
			container nhs(lhs);
			const T mul = T(1) / div;
#ifdef _DEBUG
			if (math::is_small_value(div)) {
				throw new ErrorMath::DivisionByZero();
			}
#endif
			for (size_t i = 0; i < lhs.size(); i++)
				nhs[i] *= mul;
			return nhs;
		}

		[[nodiscard]] inline friend container<T> operator * (const T& mul, const container<T>& rhs) {
			return container<T>::operator*(rhs, mul);
		}

		container& operator += (const container& rhs) {
#ifdef _DEBUG
			this->is_consist(rhs);
#endif
			container& lhs = *this;
			for (size_t i = 0; i < _size; i++)
				lhs[i] += rhs[i];
			return lhs;
		}

		container& operator -= (const container& rhs) {
#ifdef _DEBUG
			this->is_consist(rhs);
#endif
			container& lhs = *this;
			for (size_t i = 0; i < _size; i++)
				lhs[i] -= rhs[i];
			return lhs;
		}

		container& operator *= (const T& mul) {
			container& lhs = *this;
			for (size_t i = 0; i < _size; i++)
				lhs[i] *= mul;
			return lhs;
		}

		container<T>& operator /= (const T& div) {
#ifdef _DEBUG
			if (math::is_small_value(div)) {
				throw new ErrorMath::DivisionByZero();
			}
#endif
			const T mul = T(1) / div;
			container& lhs = *this;
			for (size_t i = 0; i < _size; i++)
				lhs[i] *= mul;
			return lhs;
		}

		friend bool operator == (const container<T>& lhs, const container<T>& rhs) {
#ifdef _DEBUG
			lhs.is_consist(rhs);
#endif
			const auto diff = lhs - rhs;
			if (math::is_small_value(diff.get_norm())) {
				return true;
			}
			return false;
		}

		[[nodiscard]] container<T> transpose() const {
			if (_rank == 1) return *this;
			if (_rank == 2 && _dim == 3) {
				container<T> nhs(*this);
				std::swap(nhs[3], nhs[6]);
				std::swap(nhs[4], nhs[7]);
				std::swap(nhs[5], nhs[8]);
				return nhs;
			}
			throw NoImplemetationYet();
		}
		
		T det() const {
			if (_rank != 2) {
				throw ErrorMath::ShapeMismatch();
			} else {
				if (_dim == 3) {
					return math::dim3::det_mat(this->get());
				}
				throw NoImplemetationYet();
			}
		}

		[[nodiscard]] container<T> inverse() const {
			if (_dim == 3) {
				container<T> inv_matr(3, 2);
				math::dim3::inv_mat(this->get(), inv_matr.get());
				return inv_matr;
			} else {
				// for any dim matrix
			}
			throw NoImplemetationYet();
		}

		T get_norm() const {
			T norm = T(0);
			const auto& arr = *this;
			for (size_t idx = 0; idx < _size; idx++){
				norm += arr[idx]*arr[idx];
			}
			return sqrt(norm);
		}

		[[nodiscard]] friend static container<T> get_normalize(const container<T>& m) {
			if (m.rank() == 1) {
				return m / m.get_norm();
			}
			throw NoImplemetationYet();
		}
		// lhs * rhsT
		friend  container<T> mat_scal_mat_transp(const  container<T>& lhs, const  container<T>& rhs) {
#ifdef _DEBUG
			lhs.is_consist(rhs);
#endif
			if (lhs.rank() == 2 && lhs.dim() == 3) {
				container<T> nhs(3, 2);
				math::dim3::mat_scal_mat_transp(lhs.get(), rhs.get(), nhs.get());
				return nhs;
			}
			throw NoImplemetationYet();
		}
	};

	template<typename T>
	[[nodiscard]] container<T> operator * (const container<T>& lhs, const container<T>& rhs) {
#ifdef _DEBUG
		if (lhs.dim() != rhs.dim()) {
			throw new ErrorMath::ShapeMismatch();
		}
#endif
		const size_t l_rank = lhs.rank();
		const size_t r_rank = rhs.rank();

		if (l_rank == 2 && r_rank == 2) {
			container<T> nhs(3, 2);
			math::dim3::mat_scal_mat(lhs.get(), rhs.get(), nhs.get());
			return nhs;
		} else 	if (l_rank == 1 && r_rank == 2) {
			container<T> nhs(3, 1);
			math::dim3::vect_scal_mat(lhs.get(), rhs.get(), nhs.get());
			return nhs;
		} else 	if (l_rank == 2 && r_rank == 1) {
			container<T> nhs(3, 1);
			math::dim3::mat_scal_vect(lhs.get(), rhs.get(), nhs.get());
			return nhs;
		}else if (l_rank == 1 && r_rank == 1) {
			container<T> nhs(1, 0);
			nhs[0] = math::dim3::vect_scal_vect(lhs.get(), rhs.get());
			return nhs;
		}
		throw NoImplemetationYet();
	}

	template<typename T>
	std::ostream& operator<<(std::ostream& out, const container<T>& cont) {
		out << "{ ";
		for (size_t row = 0; row < cont.size() - 1; row++) 
			out << cont[row] << ", ";
		out << cont[cont.size() - 1] << " }";
		return out;
	};
}


template <typename T>
using Basis = std::shared_ptr<const tens::container<T>>;

template<typename T>
extern const Basis<T> GLOBAL_BASIS = std::make_shared<const tens::container<T>>(3, 2, tens::FILL_TYPE::INDENT);
template<typename T>
extern const tens::container<T> IDENT_MATRIX = tens::container<T>(3, 2, tens::FILL_TYPE::INDENT);
template<typename T>
extern const tens::container<T> ZERO_MATRIX = tens::container<T>(3, 2, tens::FILL_TYPE::ZERO);

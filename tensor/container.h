#pragma once
#include <memory>
#include <array>
#include <random>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <utility>
#include "error.h"

extern std::mt19937 gen;       // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> unidistr;

namespace tens {
	template<typename T> bool is_not_small_value(T value);
	template<typename T> bool is_small_value(T value);

	enum class ARRAY_TYPE {
		ZERO,
		RANDOM,
		RANDOMUNIT
	};

	template<typename T, size_t N, size_t R>
	class container_rank;
	template<typename T, size_t N, size_t R> std::ostream& operator<< (std::ostream& o, const container_rank<T, N, R>& cont);
	
	template<typename T, size_t N>
	class container : public std::unique_ptr<T[]>
	{
	private:
		void _copy(const container& c) {
			this->_dim = c._dim;
			this->_rank = c._rank;
			this->_size = c._size;
			memcpy(this->get(), c.get(), _size * sizeof(T));
		}

		void _move(container& c) {
			this->_dim = c._dim;
			this->_rank = c._rank;
			this->_size = c._size;
			static_cast<std::unique_ptr<T[]>&>(*this) = std::move(std::move(c));
		}

	protected:
		size_t _rank;
		size_t _dim;
		size_t _size;

		void fill_rand() {
			for (size_t i = 0; i < _size; i++) {
				this->get()[i] = static_cast<T>(unidistr(gen));
			}
		};

		void fill_value(const T& val) {
			std::fill(this->get(), this->get() + _size, val);
		};
	public:
		size_t get_rank() const { return _rank; };
		size_t get_dim() const { return _dim; };
		size_t get_size() const { return _size; };

		container(size_t rank) : std::unique_ptr<T[]>(new T[(size_t)pow(N, rank)]), _dim(N), _rank(rank), _size((size_t)pow(N, rank)) {
			//_set_zero();
		};
		container(const container& c) : std::unique_ptr<T[]>(new T[c._size]), _dim(N), _rank(c._rank), _size(c._size) {
			_copy(c);
		};

		container(container&& c) noexcept : std::unique_ptr<T[]>(std::move(c)), _dim(N), _rank(c._rank), _size(c._size) {
		};

		inline container& operator= (const container& rhs) {
			this->_copy(rhs);
			return *this;
		}

		inline container& operator= (container&& rhs) noexcept {
			this->_move(rhs);
			return *this;
		}
	};


	container_rank<double, 3, 2> generate_rand_ort();
	container_rank<double, 3, 2> generate_indent_ort();

	template<typename T>
	std::unique_ptr<T[]> flat_3x3_array(const std::array<std::array<T, 3>, 3>& matrix) {
		// ------- for DIM == 3 : {00, 11, 22, 12, 02, 01, 21, 20, 10} 
		std::unique_ptr<T[]> arr = std::unique_ptr<T[]>(new T[9]);
		arr[0] = matrix[0][0];
		arr[1] = matrix[1][1];
		arr[2] = matrix[2][2];
		arr[3] = matrix[1][2];
		arr[4] = matrix[0][2];
		arr[5] = matrix[0][1];
		arr[6] = matrix[2][1];
		arr[7] = matrix[2][0];
		arr[8] = matrix[1][0];
		return arr;
	}

	template<typename T, size_t N, size_t R>
	class container_rank : public std::unique_ptr<T[]>
	{
	private:
		void _copy(const container_rank& c) {
			this->_dim = c._dim;
			this->_size = c._size;
			memcpy(this->get(), c.get(), _size * sizeof(T));
		}

		void _copy(const std::array<T, N>& arr) {
			memcpy(this->get(), &arr, _size * sizeof(T));
		}

		void _move(container_rank& c) {
			this->_dim = c._dim;
			this->_size = c._size;
			static_cast<std::unique_ptr<T[]>&>(*this) = std::move(std::move(c));
		}
	public:
		size_t _dim;
		size_t _size;

		void fill_rand() {
			for (size_t i = 0; i < _size; i++) {
				this->get()[i] = static_cast<T>(unidistr(gen));
			}
		};

		void fill_value(const T& val) {
			std::fill(this->get(), this->get() + _size, val);
		};
		size_t get_dim() const { return _dim; };
		size_t get_size() const { return _size; };

		container_rank(const std::array<std::array<T, N>, N>& matrix) : std::unique_ptr<T[]>(flat_3x3_array<T>(matrix)), _dim(N), _size((size_t)pow(N, R)) {
			assert(R == 2 || "container_rank R must be equal 2");
		};

		container_rank(const std::array<T, N>& arr) : std::unique_ptr<T[]>(new T[N]), _dim(N), _size((size_t)pow(N, R)) {
			assert(R == 1 || "container_rank R must be equal 2");
			_copy(arr);
		};

		container_rank(std::unique_ptr<T[]>&& pointer) : std::unique_ptr<T[]>(std::move(pointer)), _dim(N), _size((size_t)pow(N, R)) {
		};

		container_rank(ARRAY_TYPE type) : std::unique_ptr<T[]>(new T[(size_t)pow(N, R)]), _dim(N), _size((size_t)pow(N, R)) {
			switch (type)
			{
			case tens::ARRAY_TYPE::ZERO:
				fill_value(T(0));
				break;
			case tens::ARRAY_TYPE::RANDOM:
				fill_rand();
				break;
			case tens::ARRAY_TYPE::RANDOMUNIT:
				fill_rand();
				*this = get_normalize(*this);
				break;
			default:
				break;
			}
		};

		container_rank() : std::unique_ptr<T[]>(new T[(size_t)pow(N, R)]), _dim(N), _size((size_t)pow(N, R)) {
		};

		container_rank(const T& val) : std::unique_ptr<T[]>(new T[(size_t)pow(N, R)]), _dim(N), _size((size_t)pow(N, R)) {
			fill_value(val);
		};

		container_rank(const container_rank& c) : std::unique_ptr<T[]>(new T[c._size]), _dim(N), _size(c._size) {
			_copy(c);
		};

		container_rank(container_rank&& c) noexcept : std::unique_ptr<T[]>(std::move(c)), _dim(N), _size(c._size) {
		};

		inline container_rank& operator= (const container_rank& rhs) {
			this->_copy(rhs);
			return *this;
		}

		inline container_rank& operator= (container_rank&& rhs) noexcept {
			this->_move(rhs);
			return *this;
		}

		friend container_rank<T, N, R> operator + (const container_rank<T, N, R>& lhs, const container_rank<T, N, R>& rhs) {
			container_rank nhs;
			for (size_t i = 0; i < lhs._size; i++)
				nhs[i] = lhs[i] + rhs[i];
			return nhs;
		}

		friend container_rank<T, N, R> operator - (const container_rank<T, N, R>& lhs, const container_rank<T, N, R>& rhs) {
			container_rank nhs;
			for (size_t i = 0; i < lhs._size; i++)
				nhs[i] = lhs[i] - rhs[i];
			return nhs;
		}

		friend container_rank<T, N, R> operator * (const container_rank<T, N, R>& lhs, const T& mul) {
			container_rank nhs;
			for (size_t i = 0; i < lhs._size; i++)
				nhs[i] = mul*lhs[i];
			return nhs;
		}

		friend container_rank<T, N, R> operator / (const container_rank<T, N, R>& lhs, const T& div) {
			container_rank nhs;
			const T mul = T(1) / div;
			for (size_t i = 0; i < lhs._size; i++)
				nhs[i] = mul * lhs[i];
			return nhs;
		}

		friend container_rank<T, N, R> operator * (const T& mul, const container_rank<T, N, R>& rhs) {
			return container_rank<T, N, R>::operator*(rhs, mul);
		}

		container_rank& operator += (const container_rank& rhs) {
			container_rank& lhs = *this;
			for (size_t i = 0; i < _size; i++)
				lhs[i] += rhs[i];
			return *this;
		}

		container_rank& operator -= (const container_rank& rhs) {
			container_rank& lhs = *this;
			for (size_t i = 0; i < _size; i++)
				lhs[i] -= rhs[i];
			return *this;
		}

		container_rank operator *= (const T& mul) {
			for (size_t i = 0; i < _size; i++)
				this[i] *= mul;
			return *this;
		}

		friend bool operator == (const container_rank<T, N, R>& lhs, const container_rank<T, N, R>& rhs) {
			const auto diff = lhs - rhs;
			if (is_small_value(diff.get_norm())) {
				return true;
			}
			return false;
		}

		friend static container_rank<T, N, R> transpose(const container_rank<T, N, R>& m) {
			if (R == 1) return m;
			if (N == 3) {
				container_rank<T, 3, 2> nhs(m);
				std::swap(nhs[3], nhs[6]);
				std::swap(nhs[4], nhs[7]);
				std::swap(nhs[5], nhs[8]);
				return nhs;
			}
		}

		T get_norm() const {
			T norm = T(0);
			const auto& arr = *this;
			for (size_t idx = 0; idx < _size; idx++){
				norm += arr[idx]*arr[idx];
			}
			return sqrt(norm);
		}

		friend static container_rank<T, N, R> get_normalize(const container_rank<T, N, R>& m) {
			if (R == 1) {
				return m / m.get_norm();
			}
		}
		// lhs * rhsT
		friend  container_rank<T, N, R> mat_scal_mat_transp(const  container_rank<T, N, R>& lhs, const  container_rank<T, N, R>& rhs) {
			if (N == 3) {
				container_rank<T, N, R> nhs;
				nhs[0] = lhs[0] * rhs[0] + lhs[4] * rhs[4] + lhs[5] * rhs[5];
				nhs[1] = lhs[1] * rhs[1] + lhs[3] * rhs[3] + lhs[8] * rhs[8];
				nhs[2] = lhs[2] * rhs[2] + lhs[6] * rhs[6] + lhs[7] * rhs[7];
				nhs[3] = lhs[3] * rhs[2] + lhs[1] * rhs[6] + lhs[8] * rhs[7];
				nhs[4] = lhs[4] * rhs[2] + lhs[5] * rhs[6] + lhs[0] * rhs[7];
				nhs[5] = lhs[5] * rhs[1] + lhs[4] * rhs[3] + lhs[0] * rhs[8];
				nhs[6] = lhs[6] * rhs[1] + lhs[2] * rhs[3] + lhs[7] * rhs[8];
				nhs[7] = lhs[7] * rhs[0] + lhs[2] * rhs[4] + lhs[6] * rhs[5];
				nhs[8] = lhs[8] * rhs[0] + lhs[3] * rhs[4] + lhs[1] * rhs[5];
				return nhs;
			}
		}
	};

	template<typename T, size_t N, size_t R = 2>
	container_rank<T, N, 2> operator * (const container_rank<T, N, 2>& lhs, const container_rank<T, N, 2>& rhs) {
		container_rank<T, 3, 2> nhs;
		nhs[0] = lhs[0] * rhs[0] + lhs[4] * rhs[7] + lhs[5] * rhs[8];
		nhs[1] = lhs[1] * rhs[1] + lhs[8] * rhs[5] + lhs[3] * rhs[6];
		nhs[2] = lhs[2] * rhs[2] + lhs[6] * rhs[3] + lhs[7] * rhs[4];
		nhs[3] = lhs[3] * rhs[2] + lhs[1] * rhs[3] + lhs[8] * rhs[4];
		nhs[4] = lhs[4] * rhs[2] + lhs[5] * rhs[3] + lhs[0] * rhs[4];
		nhs[5] = lhs[5] * rhs[1] + lhs[0] * rhs[5] + lhs[4] * rhs[6];
		nhs[6] = lhs[6] * rhs[1] + lhs[7] * rhs[5] + lhs[2] * rhs[6];
		nhs[7] = lhs[7] * rhs[0] + lhs[2] * rhs[7] + lhs[6] * rhs[8];
		nhs[8] = lhs[8] * rhs[0] + lhs[3] * rhs[7] + lhs[1] * rhs[8];
		return nhs;
	}

	template<typename T, size_t N, size_t R = 1>
	container_rank<T, N, 1> operator * (const container_rank<T, N, 1>& a, const container_rank<T, N, 2>& m) {
		container_rank<T, 3, 1> nhs;
		nhs[0] = a[0] * m[0] + a[2] * m[7] + a[1] * m[8];
		nhs[1] = a[1] * m[1] + a[0] * m[5] + a[2] * m[6];
		nhs[2] = a[2] * m[2] + a[1] * m[3] + a[0] * m[4];
		return nhs;
	}

	template<typename T, size_t N, size_t R = 1>
	container_rank<T, N, 1> operator * (const container_rank<T, N, 2>& m, const container_rank<T, N, 1>& a) {
		container_rank<T, 3, 1> nhs;
		nhs[0] = a[0] * m[0] + a[2] * m[4] + a[1] * m[5];
		nhs[1] = a[1] * m[1] + a[2] * m[3] + a[0] * m[8];
		nhs[2] = a[2] * m[2] + a[1] * m[6] + a[0] * m[7];
		return nhs;
	}

	template<typename T, size_t N, size_t R>
	std::ostream& operator<<(std::ostream& out, const container_rank<T, N, R>& cont) {
		out << "{ ";
		for (size_t idx = 0; idx < cont._size-1; idx++) 
			out << cont[idx] << ", ";
		out << cont[cont._size - 1] << " }\n";
		return out;
	};


	template <typename T, size_t N>
	using Array = tens::container_rank<T, N, 1>;
	template <typename T, size_t N>
	using Matrix = tens::container_rank<T, N, 2>;
}



/* ==============================================================================================================
	containers without basis
============================================================================================================== */
namespace tens {

	template<typename T, size_t N>
	class container_array : public container_rank<T, N, 1>
	{
	public:
		container_array() {};
		container_array(const container_rank<T, N, 1>& cont) : container_rank<T, N, 1>(cont) {};
		container_array(container_rank<T, N, 1>&& cont)noexcept : container_rank<T, N, 1>(std::move(cont)) {};
	};

	template<typename T, size_t N>
	class container_matrix : public container_rank<T, N, 2>
	{
	public:
		container_matrix() {};
		container_matrix(const container_rank<T, N, 2>& cont) : container_rank<T, N, 2>(cont) {};
		container_matrix(container_rank<T, N, 2>&& cont)noexcept : container_rank<T, N, 2>(std::move(cont)){};
	};
}

template <typename T, size_t N>
using Basis = tens::container_rank<T, N, 2>;
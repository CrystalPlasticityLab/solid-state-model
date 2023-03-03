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

	template<typename T, size_t N, size_t R>
	class container_rank;

//	template<typename T, size_t N, size_t R = 2>
//	container_rank<T, N, 2> operator * (const container_rank<T, N, 2>& lhs, const container_rank<T, N, 2>& rhs);
	//template<typename T, size_t N=3, size_t R=2>
	//static container_rank<T, 3, 2> _mat_scal_mat_dim_3(const container_rank<T, 3, 2>& lhs, const container_rank<T, 3, 2>& rhs);


	template<typename T, size_t N, size_t R>
	class container_rank : public std::unique_ptr<T[]>
	{
	private:
		void _copy(const container_rank& c) {
			this->_dim = c._dim;
			this->_size = c._size;
			memcpy(this->get(), c.get(), _size * sizeof(T));
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

		container_rank(std::unique_ptr<T[]>&& pointer) : std::unique_ptr<T[]>(std::move(pointer)), _dim(N), _size((size_t)pow(N, R)) {
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


		friend static container_rank<T, N, R> transpose(const container_rank<T, N, R>& m) {
			if (R == 1) return m;
			if (N == 3) {
				container_rank<T, 3, 2> nhs;
				std::swap(nhs[3], nhs[6]);
				std::swap(nhs[4], nhs[7]);
				std::swap(nhs[5], nhs[8]);
				return nhs;
			}
		}

		friend  container_rank<T, N, R> mat_scal_mat_transp(const  container_rank<T, N, R>& lhs, const  container_rank<T, N, R>& rhs) {
			if (N == 3) {
				container_rank<T, N, R> nhs;
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
	{ // ------- for DIM == 3 : {00, 11, 22, 12, 02, 01, 21, 20, 10} 
	public:
		container_matrix() {};
		container_matrix(const container_rank<T, N, 2>& cont) : container_rank<T, N, 2>(cont) {};
		container_matrix(container_rank<T, N, 2>&& cont)noexcept : container_rank<T, N, 2>(std::move(cont)){};
	};
}

template <typename T, size_t N>
using Basis = tens::container_rank<T, N, 2>;
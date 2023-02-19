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
			for (size_t i = 0; i < _size; i++){
				this->get()[i] = static_cast<T>(unidistr(gen));
			}
		};

		void fill_value(const T& val) {
			std::fill(this->get(), this->get() + _size, val);
		};
	public:
		size_t get_rank() const { return _rank;};
		size_t get_dim() const { return _dim; };
		size_t get_size() const  { return _size; };

		container(size_t rank) : std::unique_ptr<T[]>(new T[(size_t)pow(N, rank)]), _dim(N), _rank(rank), _size((size_t)pow(N, rank)){
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

	protected:
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
		size_t get_dim() const { return _dim; };
		size_t get_size() const { return _size; };

		container_rank() : std::unique_ptr<T[]>(new T[(size_t)pow(N, R)]), _dim(N), _size((size_t)pow(N, R)) {
			fill_value(T(0));
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

		container_rank operator + (const container_rank& rhs) const {
			container_rank nhs;
			const container_rank& lhs = *this;
			for (size_t i = 0; i < _size; i++)
				nhs[i] = lhs[i] + rhs[i];
			return nhs;
		}

		container_rank operator - (const container_rank& rhs) const {
			container_rank nhs;
			const container_rank& lhs = *this->get();
			for (size_t i = 0; i < _size; i++)
				nhs[i] = lhs[i] - rhs[i];
			return nhs;
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
	};

	template<typename T> bool is_not_small_value(T value);
	template<typename T> bool is_small_value(T value);
};
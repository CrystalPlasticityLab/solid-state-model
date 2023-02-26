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

		container_rank operator + (const container_rank& rhs) const {
			container_rank nhs;
			for (size_t i = 0; i < _size; i++)
				nhs[i] = this[i] + rhs[i];
			return nhs;
		}

		container_rank operator - (const container_rank& rhs) const {
			container_rank nhs;
			for (size_t i = 0; i < _size; i++)
				nhs[i] = this[i] - rhs[i];
			return nhs;
		}

		container_rank operator * (const T& mul) const {
			container_rank nhs;
			for (size_t i = 0; i < _size; i++)
				nhs[i] = mul * this[i];
			return nhs;
		}

		container_rank& operator += (const container_rank& rhs) {
			for (size_t i = 0; i < _size; i++)
				this[i] += rhs[i];
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
	};
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
		friend T operator *(const container_array<T, N>& lhs, const container_array<T, N>& rhs) {
			T res(0);
			for (size_t col = 0; col < N; col++) {
				res += lhs[col] * rhs[col];
			}
			return res;
		}
	};

	template<typename T, size_t N>
	class container_matrix : public container_rank<T, N, 2>
	{ // ------- for DIM == 3 : {00, 11, 22, 12, 02, 01, 21, 20, 10}  

		static friend void _mat_transp_dim_3(container_matrix<T, 3>& m) {
			std::swap(m[3], m[6]);
			std::swap(m[4], m[7]);
			std::swap(m[5], m[8]);
		}
		static friend container_matrix<T, 3> _mat_transp_dim_3(const container_matrix<T, 3>& m) {
			container_matrix<T, N> nhs(m);
			std::swap(nhs[3], nhs[6]);
			std::swap(nhs[4], nhs[7]);
			std::swap(nhs[5], nhs[8]);
			return nhs;
		}
		static friend container_matrix<T, 3> _mat_scal_mat_dim_3(const container_matrix<T, 3>& lhs, const container_matrix<T, 3>& rhs) {
			container_matrix<T, N> nhs;
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
		static friend container_array<T, N> _mat_scal_vec_dim_3(const container_matrix<T, N>& m, const container_array<T, N>& a) {
			container_array<T, N> nhs;
			nhs[0] = a[0] * m[0] + a[2] * m[4] + a[1] * m[5];
			nhs[1] = a[1] * m[1] + a[2] * m[3] + a[0] * m[8];
			nhs[2] = a[2] * m[2] + a[1] * m[6] + a[0] * m[7];
			return nhs;

		}
		static friend container_array<T, N> _vec_scal_mat_dim_3(const container_array<T, N>& a, const container_matrix<T, N>& m) {
			container_array<T, N> nhs;
			nhs[0] = a[0] * m[0] + a[2] * m[7] + a[1] * m[8];
			nhs[1] = a[1] * m[1] + a[0] * m[5] + a[2] * m[6];
			nhs[2] = a[2] * m[2] + a[1] * m[3] + a[0] * m[4];
			return nhs;
		}
	public:
		container_matrix() {};
		container_matrix(const container_rank<T, N, 2>& cont) : container_rank<T, N, 2>(cont) {};
		container_matrix(container_rank<T, N, 2>&& cont)noexcept : container_rank<T, N, 2>(std::move(cont)){};
		friend container_matrix<T, N> operator * (const container_matrix<T, N>& lhs, const container_matrix<T, N>& rhs) {
			if (N == 3) {
				return _mat_scal_mat_dim_3(lhs, rhs);
			}
		}

		friend container_array<T, N> operator *(const container_matrix<T, N>& m, const container_array<T, N>& a) {
			if (N == 3) {
				return _mat_scal_vec_dim_3(m, a);
			}
		}
		friend container_array<T, N> operator *(const container_array<T, N>& a, const container_matrix<T, N>& m) {
			if (N == 3) {
				return _vec_scal_mat_dim_3(a, m);
			}
		}
		friend static container_matrix<T, N> transpose(const container_matrix<T, N>& m) {
			if (N == 3) {
				return _mat_transp_dim_3(m);
			}
		}
		void transpose() {
			if (N == 3) {
				_mat_transp_dim_3(*this);
			}
		}
	};
}


/* ==============================================================================================================
	abstract class WITH basis
============================================================================================================== */
namespace tens {
	template<typename T, size_t N, size_t R>
	class basis_object {
		std::unique_ptr<container_rank<T, N, R>> _comp;
		std::shared_ptr<container_matrix<T, N>> _basis;

		void _copy(const basis_object<T, N, R>& src) {
			if (_basis == nullptr) {
				_comp = std::make_unique<container_rank<T, N, R>>(container_rank<T, N, R>(*src._comp));
				_basis = std::make_shared<container_matrix<T, N>>(container_matrix<T, N>(*src._basis));
			} else {
				if (_basis == src._basis) {
					*_comp = *src._comp;
				} else {
					*_comp = src.get_comp_at_basis(_basis);
				}
			}
		}

		void _move(basis_object<T, N, R>&& src) {
			if (_basis == nullptr) {
				_comp = std::move(src._comp);
				_basis = std::move(src._basis);
			} else {
				if (_basis == src._basis) {
					_comp = std::move(src._comp);
				}else {
					*_comp = src.get_comp_at_basis(_basis);
					src._comp.reset();
				}
				src._basis.reset();
			}
		}

		void _reset_basis(const container_matrix<T, N>& basis) {
			_basis.reset();
			_basis = std::make_shared<container_matrix<T, N>>(container_matrix<T, N>(basis));
		}

		void _reset_basis(const std::shared_ptr<container_matrix<T, N>>& pbasis) {
			_basis.reset();
			_basis = pbasis;
		}
	protected:
		virtual container_rank<T, N, R> get_comp_at_basis(const std::shared_ptr<container_matrix<T, N>>& pbasis) const = 0;

		void change_basis(const container_matrix<T, N>& basis) {
			_comp = this->get_comp_at_basis(basis);
			_reset_basis(basis);
		}
		void move_basis(const container_matrix<T, N>& basis) {
			_reset_basis(basis);
		}
		void change_basis(const std::shared_ptr<container_matrix<T, N>>& pbasis) {
			_comp = this->get_comp_at_basis(*pbasis);
			_reset_basis(pbasis);
		}
		void move_basis(const std::shared_ptr<container_matrix<T, N>>& pbasis) {
			_reset_basis(pbasis);
		}
		virtual size_t get_rank() = 0;
		const container_rank<T, N, R>& get_comp_ref() const {
			return *this->_comp;
		}
	public:

		basis_object(const basis_object<T, N, R>& basis_obj) { // copy ctor
			_copy(basis_obj);
		}

		basis_object(basis_object<T, N, R>&& basis_obj) noexcept { // move ctor
			_move(std::move(basis_obj));
		}

		basis_object(const container_rank<T, N, R>& comp, const container_matrix<T, N>& basis) {
			_comp = std::make_unique<container_rank<T, N, R>>(comp);
			_basis = std::make_shared<container_matrix<T, N>>(basis);
		}

		basis_object(const container_rank<T, N, R>& comp, const std::shared_ptr<container_matrix<T, N>>& pbasis) {
			_comp = std::make_unique<container_rank<T, N, R>>(comp);
			_basis = pbasis;
		}

		basis_object& operator = (const basis_object<T, N, R>& rhs) { // copy assign
			_copy(rhs);
			return *this;
		}

		basis_object& operator = (basis_object<T, N, R>&& rhs) noexcept { // move assign
			_move(std::move(rhs));
			return *this;
		}

		const std::shared_ptr<container_matrix<T, N>>& get_basis_ref() const {
			return this->_basis;
		}
	};
}


/* ==============================================================================================================
	Objects WITH basis
============================================================================================================== */
namespace tens {
	template<typename T, size_t N>
	class container_vector: private basis_object<T, N, 1> {
	protected:
		virtual size_t get_rank() override {
			return 1;
		};

		virtual container_rank<T, N, 1> get_comp_at_basis(const std::shared_ptr<container_matrix<T, N>>& basis) const override {
			const auto& basis_new = basis;
			const auto& basis_cur = this->get_basis_ref();
			const container_array<T, N>& comp = static_cast<const container_array<T, N>&>(this->get_comp_ref());
			if (basis_new == basis_cur) {
				return comp;
			}
			else {
				const container_matrix<T, N>& Rl = *this->get_basis_ref();
				const container_matrix<T, N>& Rr = transpose(*basis.get());
				return static_cast<container_rank<T, N, 1>>((comp * Rl) * Rr);
			}
		}
	public:

		container_array<T, N> get_comp(const std::shared_ptr<container_matrix<T, N>>& pbasis) const {
			return static_cast<const container_array<T, N>&>(this->get_comp_at_basis(pbasis));
		}
		container_vector(const container_vector<T, N>& v) : basis_object<T, N, 1>(v) {};
		container_vector(container_vector<T, N>&& v) noexcept : basis_object<T, N, 1>(std::move(v)) {};
		container_vector(const container_array<T, N>& comp, const container_matrix<T, N>& basis) : basis_object<T, N, 1>(comp, basis) {};
		container_vector(const container_array<T, N>& comp, const std::shared_ptr<container_matrix<T, N>>& pbasis) : basis_object<T, N, 1>(comp, pbasis) {};
		template <size_t R>
		container_vector(const container_array<T, N>& comp, const basis_object<T, N, R>& object) : basis_object<T, N, 1>(comp, object.get_basis_ref()) {};

		friend T operator * (const container_vector<T, N>& lhs, const container_vector<T, N>& rhs) {
			const container_array<T, N>& lhsa = static_cast<const container_array<T, N>&>(lhs.get_comp_ref());
			const container_array<T, N>& rhsa = static_cast<const container_array<T, N>&>(rhs.get_comp_at_basis(lhs.get_basis_ref()));// container_array<T, N>(std::move(rhs.get_comp_at_basis(lhs)));
			return  lhsa * rhsa;
		}

		container_vector& operator = (const container_vector& rhs) {
			return static_cast<container_vector&>(basis_object<T, N, 1>::operator=(static_cast<const basis_object<T, N, 1>&>(rhs)));
		};

		container_vector& operator = (container_vector&& rhs) noexcept {
			return static_cast<container_vector&>(basis_object<T, N, 1>::operator=(static_cast<basis_object<T, N, 1>&&>(rhs)));
		};
	};

	template<typename T, size_t N>
	class container_tensor: private basis_object<T, N, 2> {
	protected:
		virtual size_t get_rank() override {
			return 2;
		};

		virtual container_rank<T, N, 2> get_comp_at_basis(const std::shared_ptr<container_matrix<T, N>>& basis) const override {
			const auto basis_new = basis;
			const auto basis_cur = this->get_basis_ref();
			const container_matrix<T, N>& comp = static_cast<const container_matrix<T, N>&>(this->get_comp_ref());
			if (basis_new == basis_cur) {
				return comp;
			} else {
				const container_matrix<T, N>& Rl = *this->get_basis_ref();
				const container_matrix<T, N>& Rr = transpose(*basis.get());
				container_matrix<T, N> op = Rl*Rr;
				return container_rank<T, N, 2>(transpose(op) * comp * op);// op^t * (*this) * op
			}
		}
	public:
		container_matrix<T, N> get_comp() const {
			return static_cast<container_matrix<T, N>>(this->get_comp_ref());
		}
		container_tensor(const container_tensor<T, N>& t) : basis_object<T, N, 2>(t) {};
		container_tensor(container_tensor<T, N>&& t) noexcept : basis_object<T, N, 2>(std::move(t)) {};
		container_tensor(const container_matrix<T, N>& comp, const container_matrix<T, N>& basis) : basis_object<T, N, 2>(comp, basis) {};
		container_tensor(const container_matrix<T, N>& comp, const std::shared_ptr<container_matrix<T, N>>& pbasis) : basis_object<T, N, 2>(comp, pbasis) {};
		template <size_t R>
		container_tensor(const container_matrix<T, N>& comp, const basis_object<T, N, R>& object) : basis_object<T, N, 2>(comp, object.get_basis_ref()) {};

		friend container_tensor<T, N> operator * (const container_tensor<T, N>& lhs, const container_tensor<T, N>& rhs) {
			const container_matrix<T, N>& lhsa = static_cast<const container_matrix<T, N>&>(lhs.get_comp_ref());
			if (lhs.get_basis_ref() == rhs.get_basis_ref()) {
				const container_matrix<T, N>& rhsa = static_cast<const container_matrix<T, N>&>(rhs.get_comp_ref());
				return container_tensor<T, N>(lhsa * rhsa, lhs);
			}
			else {
				const container_matrix<T, N>& rhsa = static_cast<const container_matrix<T, N>&>(rhs.get_comp_at_basis(lhs.get_basis_ref()));
				return container_tensor<T, N>(lhsa * rhsa, lhs);
			}
		}

		friend container_vector<T, N> operator * (const container_tensor<T, N>& lhs, const container_vector<T, N>& rhs) {
			container_matrix<T, N> lhsa = lhs.get_comp();
			container_array<T, N> rhsa = rhs.get_comp(lhs.get_basis_ref());
			return container_vector<T, N>(lhsa * rhsa, lhs.get_basis_ref());
		}

		friend container_vector<T, N> operator * (const container_vector<T, N>& lhs, const container_tensor<T, N>& rhs) {
			container_array<T, N> lhsa = lhs.get_comp(rhs.get_basis_ref());
			container_matrix<T, N> rhsa = rhs.get_comp();
			return container_vector<T, N>(lhsa * rhsa, rhs.get_basis_ref());
		}

		container_tensor& operator = (const container_tensor& rhs) {
			return static_cast<container_tensor&>(basis_object<T, N, 2>::operator=(static_cast<const basis_object<T, N, 2>&>(rhs)));
		}; 

		container_tensor& operator = (container_tensor&& rhs) noexcept {
			return static_cast<container_tensor&>(basis_object<T, N, 2>::operator=(static_cast<basis_object<T, N, 2>&&>(rhs)));
		};
	};

	template<typename T> bool is_not_small_value(T value);
	template<typename T> bool is_small_value(T value);
};
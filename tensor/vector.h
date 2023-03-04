#pragma once

#include <iostream>
#include "array.h"
#include "matrix.h"
#include "basis.h"

namespace tens {
	template<typename T, size_t N>
	class vector : private array<T, N>,
				   public basis<T, N>
	{
		typedef  array   <T, N>              _vector;
		typedef  matrix  <T, N>              _matrix;
		typedef  basis<T, N>  _handler;
	public:
		size_t get_rank() const override { return static_cast<const container<T, N>*>(this)->get_rank(); };
		virtual void       move_to_basis    (const _handler& m) override { basis<T, N>::move_to_basis(m); };
		virtual void       change_basis     (const _handler& m) override;
		array<T, N>  get_comp_at_basis(const _handler& m) const; // calc comp of this at basis

		vector(const _vector& comp, const _handler& basis) : _vector(comp), _handler(basis) {}
		vector(const _vector& comp, const _matrix& basis) : _vector(comp), _handler(basis) {};
		vector(const  vector& vect)  : _vector(vect), _handler(vect) {}; // copy constructor
		vector(vector&& vect) noexcept: _vector(static_cast<_vector&&>(vect)), _handler(static_cast<_handler&&>(vect))  {}; // move constructor

		inline vector  operator - ();
		inline bool    operator== (const vector<T, N>& v) const;
		inline vector& operator = (const vector<T, N>& v);
		inline vector& operator = (vector<T, N>&& v) noexcept;
		inline vector  operator + (const vector<T, N>& v) const;
		inline vector  operator - (const vector<T, N>& v) const;
		inline T       operator * (const vector<T, N>& v) const;
		inline vector  operator * (const T& mult) const;
		inline vector& operator +=(const vector<T, N>& v);
		inline vector& operator -=(const vector<T, N>& v);
		inline vector& operator *=(const T& mult);
		inline vector& operator /=(const T& mult);
		inline T norm () const{
			return array<T, N>::norm();
		}

		inline virtual void normalize() override;

		friend vector<T, N> vector_product(const vector<T, N>& _lhs, const vector<T, N>& _rhs) {
			static_assert(N == 3 && " vector_product is implemented only for 3-dimensional vectors.");
			vector<T, N> res(_lhs);
			array<T, N>& comp = static_cast<array<T, N>&>(res);
			const array<T, N>& lhs = static_cast<const array<T, N>&>(_lhs);
			const array<T, N> rhs = _rhs.get_comp_at_basis(_lhs);
			comp = array_vector_product(lhs, rhs);
			return res;
		}
	};

	template<typename T, size_t N>
	inline void vector<T, N >::normalize() {
		array<T, N>::normalize(*this);
	}

	template<typename T, size_t N>
	void vector<T, N>::change_basis(const _handler& m) {
		static_cast<array<T, N>&>(*this) = get_comp_at_basis(m);
		this->set_basis(m);
	}

	template<typename T, size_t N>
	array<T, N> vector<T, N>::get_comp_at_basis(const   _handler& m) const {
		if (*this == m) {
			return static_cast<array<T, N>> (*this);
		}
		else {
			const matrix<T, N>& Rl = *static_cast<const _handler*>(this)->get();
			const matrix<T, N>& Rr = m.get()->transpose();
			const array<T, N>& comp = static_cast<const array<T, N>&> (*this);
			return (comp * Rl) * Rr;
		}
	}

	template<typename T, size_t N>
	vector<T, N> vector<T, N>::operator - () {
		return  vector<T, N>(array<T, N>::operator-(), *this);
	}

	template<typename T, size_t N>
	bool vector<T, N>::operator== (const vector<T, N>& v) const {
		const vector<T, N> diff = (*this) - v;
		if (is_small_value(diff*diff)) {
			return true;
		}
		return false;
	};

	template<typename T, size_t N>
	vector<T, N>& vector<T, N>::operator = (const vector<T, N>& v) {
		_vector ::operator = (static_cast<const _vector &>(v));
		_handler::operator = (static_cast<const _handler&>(v));
		return *this;
	};

	template<typename T, size_t N>
	vector<T, N>& vector<T, N>::operator =  (vector&& v) noexcept {
		_vector ::operator = (static_cast<_vector &&>(v));
		_handler::operator = (static_cast<_handler&&>(v));
		return*this;
	};

	template<typename T, size_t N>
	vector<T, N> vector<T, N>::operator + (const vector<T, N>& v) const {
		return  vector<T, N>(array<T, N>::operator+(v.get_comp_at_basis(*this)), *this);
	}

	template<typename T, size_t N>
	vector<T, N> vector<T, N>::operator - (const vector<T, N>& v) const {
		return  vector<T, N>(array<T, N>::operator-(v.get_comp_at_basis(*this)), *this);
	}

	template<typename T, size_t N>
	vector<T, N>& vector<T, N>::operator +=(const vector<T, N>& v) {
		array<T, N>::operator+=(v.get_comp_at_basis(*this));
		return *this;
	}

	template<typename T, size_t N>
	vector<T, N>& vector<T, N>::operator *=(const T& mul) {
		array<T, N>::operator*=(mul);
		return *this;
	}

	template<typename T, size_t N>
	vector<T, N>& vector<T, N>::operator /=(const T& div) {
		array<T, N>::operator/=(div);
		return *this;
	}

	template<typename T, size_t N>
	vector<T, N>& vector<T, N>::operator -=(const vector<T, N>& v) {
		array<T, N>::operator-=(v.get_comp_at_basis(*this));
		return *this;
	}

	template<typename T, size_t N>
	T vector<T, N>::operator * (const vector<T, N>& v) const {
		return array<T, N>::operator*(v.get_comp_at_basis(*this));
	}

	template<typename T, size_t N>
	vector<T, N> vector<T, N>::operator * (const T& mult) const {
		return vector<T, N>(array<T, N>::operator*(mult), *this);
	}


	// output component of tensor at GLOBAL_DEFAULT_BASIS basis
	template<typename T, std::size_t N>
	std::ostream& operator<<(std::ostream& out, const vector<T, N>& v) {
		array<T, N> c = v.get_comp_at_basis(GLOBAL_DEFAULT_BASIS<T, N>);
		out << c;
		return out;
	};
};
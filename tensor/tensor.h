#pragma once
#include "matrix.h"
#include "quat.h"
#include "vector.h"
#include "basis.h"

namespace tens {
	template<typename T, std::size_t N> class Tensor;
	template<typename T, std::size_t N> class vector;
	template<typename T, std::size_t N> std::ostream& operator<< (std::ostream& o, const Tensor<T, N>& t);
	template<typename T, std::size_t N> Tensor<T,N> outer_product (const vector<T, N>& lhs, const vector<T, N>& rhs);

	template<typename T, size_t N>
	class Tensor : private matrix<T, N>,
				private shared_handler_basis<T, N>
	{
		typedef  matrix  <T, N>        _matrix;
		typedef  shared_handler_basis<T, N> _handler;
	public:
		virtual void       move_to_basis   (const _handler& m) override { shared_handler_basis<T, N>::move_to_basis(m); };
		virtual void       change_basis     (const _handler& m) override;
		matrix<T, N>       get_comp_at_basis(const _handler& m) const; // calc comp of this at basis

		Tensor(const _matrix& comp, const _handler& basis) : _matrix(comp), _handler(basis) {};
		Tensor(const  Tensor&  tens) : _matrix(tens), _handler(static_cast<const _handler&>(tens)) {}; // copy constructor
		Tensor(Tensor&& tens) noexcept : _matrix(static_cast<_matrix&&>(tens)), _handler(static_cast<_handler&&>(tens)) {}; // move constructor

		friend std::ostream& operator<< <>(std::ostream& out, const Tensor& t);
		inline Tensor& operator = (const Tensor<T,N>& t);
		inline Tensor& operator = (Tensor&& t) noexcept;
		inline bool    operator== (const Tensor& t) const;
		inline Tensor  operator + (const Tensor& t) const;
		inline Tensor  operator - (const Tensor& t) const;
		inline Tensor  operator * (const Tensor& t) const;
		inline Tensor& operator +=(const Tensor& t);
		inline Tensor& operator -=(const Tensor& t);
		inline Tensor& operator *=(const Tensor& t);

		T       convolution(const Tensor<T, N>& rhs) const;
		friend Tensor<T, N> outer_product(const vector<T, N>& lhs, const vector<T, N>& rhs);
		friend vector<T, N> operator * (const Tensor& t, const vector<T, N>& v) {
			return vector<T, N>(v.comp_at_basis(t.get_basis()) * static_cast<matrix<T, N>>(t), t.get_basis()); }

		friend vector<T, N>  operator * (const vector<T, N>& v, const Tensor& t) {
			return vector<T, N>(static_cast<matrix<T, N>>(t) * v.comp_at_basis(t.get_basis()), t.get_basis()); }
	};

	template<typename T, std::size_t N>
	bool Tensor<T, N>::operator==(const Tensor<T, N>& rhs) const {
		const Tensor<T, N> diff = rhs - *this;
		if (is_small_value(diff.convolution(diff))) {
			return true;
		}
		return false;
	}

	template<typename T, size_t N>
	void Tensor<T, N>::change_basis(const _handler& m) {
		static_cast<matrix<T, N>&>(*this) = get_comp_at_basis(m);
		this->set_basis(m);
	}

	template<typename T, size_t N>
	Tensor<T, N>& Tensor<T, N>::operator = (const Tensor<T, N>& t)  {
		_matrix ::operator = (static_cast<const _matrix &>(t));
		_handler::operator = (static_cast<const _handler&>(t));
		return *this;
	};

	template<typename T, size_t N>
	Tensor<T, N>& Tensor<T, N>::operator =  (Tensor&& t) noexcept {
		_matrix ::operator = (static_cast<_matrix &&>(t));
		_handler::operator = (static_cast<_handler&&>(t));
		return*this;
	};

	template<typename T, size_t N>
	Tensor<T, N>  Tensor<T, N>::operator + (const Tensor<T, N>& t) const {
		return Tensor<T, N>(matrix<T, N>::operator+(t.get_comp_at_basis(*this)), *this);
	}

	template<typename T, size_t N>
	Tensor<T, N>  Tensor<T, N>::operator - (const Tensor<T, N>& t) const {
		return Tensor<T, N>(matrix<T, N>::operator-(t.get_comp_at_basis(*this)), *this);
	}

	template<typename T, size_t N>
	Tensor<T, N>  Tensor<T, N>::operator * (const Tensor<T, N>& t) const  {
		return Tensor<T, N>(matrix<T, N>::operator*(t.get_comp_at_basis(*this)), *this);
	}

	template<typename T, size_t N>
	Tensor<T, N>& Tensor<T, N>::operator +=(const Tensor<T, N>& t)  {
		matrix<T, N>::operator+=(t.get_comp_at_basis(*this));
		return *this;
	}

	template<typename T, size_t N>
	Tensor<T, N>& Tensor<T, N>::operator -=(const Tensor<T, N>& t)  {
		matrix<T, N>::operator-=(t.get_comp_at_basis(*this));
		return *this;
	}

	template<typename T, size_t N>
	Tensor<T, N>& Tensor<T, N>::operator *=(const Tensor<T, N>& t){
		matrix<T, N>::operator*=(t.get_comp_at_basis(*this));
		return *this;
	}

	template<typename T, std::size_t N>
	T   Tensor<T, N>::convolution(const Tensor<T, N>& rhs) const {
		return static_cast<matrix<T, N>>(*this).convolution(rhs.get_comp_at_basis(*this));
	}
	/* return components of this in the target basis m*/
	template<typename T, size_t N>
	matrix<T, N> Tensor<T, N>::get_comp_at_basis(const   shared_handler_basis<T, N>& m) const {
		if (*this == m) {
			return static_cast<matrix<T, N>> (*this);
		}
		else {
			matrix<T, N> op   = *this->get() * m.get()->transpose();
			const matrix<T, N>& comp = static_cast<const matrix<T, N>&> (*this);
			return this->transform(TRANSPOSE::TRUE, op, TRANSPOSE::FALSE); // op^t * (*this) * op
		}
	}

	template<typename T, size_t N>
	Tensor<T, N> outer_product(const vector<T, N>& _lhs, const vector<T, N>& _rhs){
		tens::array<T, N> lhs = _lhs.get_comp_at_basis(GLOBAL_DEFAULT_BASIS<T, N>);
		tens::array<T, N> rhs = _rhs.get_comp_at_basis(GLOBAL_DEFAULT_BASIS<T, N>);
		matrix  <T, N> m = lhs.outer_product(rhs);
		return Tensor<T, N>(m, GLOBAL_DEFAULT_BASIS<T, N>);
	}

	// output component of Tensor at GLOBAL_DEFAULT_BASIS basis
	template<typename T, std::size_t N>
	std::ostream& operator<<(std::ostream& out, const Tensor<T, N>& t) {
		matrix<T, N> m = t.get_comp_at_basis(GLOBAL_DEFAULT_BASIS<T, N>);
		out << m;
		return out;
	};
};
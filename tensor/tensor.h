#pragma once
#include "matrix.h"
#include "quat.h"
#include "vector.h"
#include "basis.h"

namespace tens {
	template<typename T, std::size_t N> class tensor;
	template<typename T, std::size_t N> class vector;
	template<typename T, std::size_t N> std::ostream& operator<< (std::ostream& o, const tensor<T, N>& t);
	template<typename T, std::size_t N> tensor<T,N> outer_product (const vector<T, N>& lhs, const vector<T, N>& rhs);

	template<typename T, size_t N>
	class tensor : private matrix<T, N>,
				   public shared_handler_basis<T, N>
	{
		typedef  matrix  <T, N>        _matrix;
		typedef  shared_handler_basis<T, N> _handler;
	protected:
		virtual void       move_to_basis(const _handler& m) override { shared_handler_basis<T, N>::move_to_basis(m); };
		virtual void       change_basis(const _handler& m) override;
	public:
		matrix<T, N>       get_comp_at_basis(const _handler& m) const; // calc comp of this at basis

		tensor(const _matrix& comp, const _handler& basis) : _matrix(comp), _handler(basis) {};
		tensor(const  tensor&  tens) : _matrix(tens), _handler(static_cast<const _handler&>(tens)) {}; // copy constructor
		tensor(tensor&& tens) noexcept : _matrix(static_cast<_matrix&&>(tens)), _handler(static_cast<_handler&&>(tens)) {}; // move constructor

		friend std::ostream& operator<< <>(std::ostream& out, const tensor& t);
		inline tensor& operator = (const tensor& t);
		inline tensor& operator = (tensor&& t) noexcept;
		inline bool    operator== (const tensor& t) const;
		inline tensor  operator + (const tensor& t) const;
		inline tensor  operator - (const tensor& t) const;
		inline tensor  operator * (const tensor& t) const;
		inline tensor& operator +=(const tensor& t);
		inline tensor& operator -=(const tensor& t);
		inline tensor& operator *=(const tensor& t);

		T       convolution(const tensor<T, N>& rhs) const;
		friend tensor<T, N> outer_product(const vector<T, N>& lhs, const vector<T, N>& rhs);
		friend vector<T, N> operator * (const tensor& t, const vector<T, N>& v) {
			return vector<T, N>(v.comp_at_basis(t.get_basis()) * static_cast<matrix<T, N>>(t), t.get_basis()); }

		friend vector<T, N>  operator * (const vector<T, N>& v, const tensor& t) {
			return vector<T, N>(static_cast<matrix<T, N>>(t) * v.comp_at_basis(t.get_basis()), t.get_basis()); }
	};

	template<typename T, std::size_t N>
	bool tensor<T, N>::operator==(const tensor<T, N>& rhs) const {
		const tensor<T, N> diff = rhs - *this;
		if (is_small_value(diff.convolution(diff))) {
			return true;
		}
		return false;
	}

	template<typename T, size_t N>
	void tensor<T, N>::change_basis(const _handler& m) {
		static_cast<matrix<T, N>&>(*this) = get_comp_at_basis(m);
		this->set_basis(m);
	}

	template<typename T, size_t N>
	tensor<T, N>& tensor<T, N>::operator = (const tensor<T, N>& t)  {
		_matrix ::operator = (static_cast<const _matrix &>(t));
		_handler::operator = (static_cast<const _handler&>(t));
		return *this;
	};

	template<typename T, size_t N>
	tensor<T, N>& tensor<T, N>::operator =  (tensor&& t) noexcept {
		_matrix ::operator = (static_cast<_matrix &&>(t));
		_handler::operator = (static_cast<_handler&&>(t));
		return*this;
	};

	template<typename T, size_t N>
	tensor<T, N>  tensor<T, N>::operator + (const tensor<T, N>& t) const {
		return tensor<T, N>(matrix<T, N>::operator+(t.get_comp_at_basis(*this)), *this);
	}

	template<typename T, size_t N>
	tensor<T, N>  tensor<T, N>::operator - (const tensor<T, N>& t) const {
		return tensor<T, N>(matrix<T, N>::operator-(t.get_comp_at_basis(*this)), *this);
	}

	template<typename T, size_t N>
	tensor<T, N>  tensor<T, N>::operator * (const tensor<T, N>& t) const  {
		return tensor<T, N>(matrix<T, N>::operator*(t.get_comp_at_basis(*this)), *this);
	}

	template<typename T, size_t N>
	tensor<T, N>& tensor<T, N>::operator +=(const tensor<T, N>& t)  {
		matrix<T, N>::operator+=(t.get_comp_at_basis(*this));
		return *this;
	}

	template<typename T, size_t N>
	tensor<T, N>& tensor<T, N>::operator -=(const tensor<T, N>& t)  {
		matrix<T, N>::operator-=(t.get_comp_at_basis(*this));
		return *this;
	}

	template<typename T, size_t N>
	tensor<T, N>& tensor<T, N>::operator *=(const tensor<T, N>& t){
		matrix<T, N>::operator*=(t.get_comp_at_basis(*this));
		return *this;
	}

	template<typename T, std::size_t N>
	T   tensor<T, N>::convolution(const tensor<T, N>& rhs) const {
		return static_cast<matrix<T, N>>(*this).convolution(rhs.get_comp_at_basis(*this));
	}
	/* return components of this in the target basis m*/
	template<typename T, size_t N>
	matrix<T, N> tensor<T, N>::get_comp_at_basis(const   shared_handler_basis<T, N>& m) const {
		if (*this == m) {
			return static_cast<matrix<T, N>> (*this);
		}
		else {
			matrix<T, N> op   = *static_cast<const shared_handler_basis<T, N>*>(this)->get() * m.get()->transpose();
			const matrix<T, N>& comp = static_cast<const matrix<T, N>&> (*this);
			return this->transform(TRANSPOSE::TRUE, op, TRANSPOSE::FALSE); // op^t * (*this) * op
		}
	}

	template<typename T, size_t N>
	tensor<T, N> outer_product(const vector<T, N>& _lhs, const vector<T, N>& _rhs){
		tens::array<T, N> lhs = _lhs.get_comp_at_basis(GLOBAL_DEFAULT_BASIS<T, N>);
		tens::array<T, N> rhs = _rhs.get_comp_at_basis(GLOBAL_DEFAULT_BASIS<T, N>);
		matrix  <T, N> m = lhs.outer_product(rhs);
		return tensor<T, N>(m, GLOBAL_DEFAULT_BASIS<T, N>);
	}

	// output component of tensor at GLOBAL_DEFAULT_BASIS basis
	template<typename T, std::size_t N>
	std::ostream& operator<<(std::ostream& out, const tensor<T, N>& t) {
		matrix<T, N> m = t.get_comp_at_basis(GLOBAL_DEFAULT_BASIS<T, N>);
		out << m;
		return out;
	};
};
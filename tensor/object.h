#pragma once
#include "container.h"

namespace tens {

	template<typename T, size_t DIM, size_t RANK>
	class object;

	template<typename T, size_t DIM, size_t RANK>
	std::ostream& operator<< (std::ostream& o, const object<T, DIM, RANK>& b);

	template<typename T, size_t DIM, size_t RANK = 2>
	object<T, DIM, RANK> Tensor(const container<T, DIM, RANK>& m, const std::shared_ptr<const container<T, DIM, 2>>& _basis) {
		if (RANK != 2) {
			throw new NoImplemetationYet();
		}
		return object<T, DIM, RANK>(m, _basis);
	}

	template<typename T, size_t DIM, size_t RANK = 2>
	object<T, DIM, 2> Tensor(const object<T, DIM, RANK>& v) {
		return object<T, DIM, RANK>(v);
	}

	template<typename T, size_t DIM, size_t RANK = 2>
	object<T, DIM, 2> Tensor(const std::array<std::array<T, DIM>, DIM>& m, const std::shared_ptr<const container<T, DIM, 2>>& _basis) {
		return object<T, DIM, RANK>(Matrix<T, DIM>(m), _basis);
	}

	template<typename T, size_t DIM, size_t RANK = 1>
	object<T, DIM, RANK> Vector(const std::array<T, DIM>& _arr, const std::shared_ptr<const container<T, DIM, 2>>& _basis) {
		return object<T, DIM, 1>(Array<T, DIM>(_arr), _basis);
	}

	template<typename T, size_t DIM, size_t RANK = 1>
	object<T, DIM, RANK> Vector(const container<T, DIM, RANK>& a, const std::shared_ptr<const container<T, DIM, 2>>& _basis) {
		if (RANK != 1) {
			throw new NoImplemetationYet();
		}
		return object<T, DIM, RANK>(a, _basis);
	}

	template<typename T, size_t DIM, size_t RANK= 1>
	object<T, DIM, RANK> Vector(const object<T, DIM, RANK>& v) {
		if (RANK != 1) {
			throw new NoImplemetationYet();
		}
		return object<T, DIM, RANK>(v);
	}

	enum class DEFAULT_ORTH_BASIS {
		RANDOM,
		INDENT
	};

	template<typename T, size_t DIM, size_t RANK = 2>
	static object<T, DIM, RANK> create_basis(size_t N, DEFAULT_ORTH_BASIS type);

	template<typename T, size_t DIM, size_t RANK = 2>
	static const std::shared_ptr<container<T, DIM, RANK>> EMPTY_BASIS = std::shared_ptr<container<T, DIM, RANK>>();

	template<typename T, size_t DIM, size_t RANK>
	class object {

		std::unique_ptr<container<T, DIM, RANK>> _comp;
		Basis<T, DIM> _basis;

		void _copy(const object<T, DIM, RANK>& src) {
			if (_basis == nullptr) { // it has no object, so we create it
				_comp = std::make_unique<container<T, DIM, RANK>>(container<T, DIM, RANK>(*src._comp));
				if (src._basis) {
					_basis = std::make_shared<container<T, DIM, 2>>(container<T, DIM, 2>(*src._basis));
				}
			}
			else { // has object, recalc comp at this object if it is necessery
				if (_basis == src._basis) { // object is the same, just copy comp
					*_comp = *src._comp;
				}
				else { // other object - recalc comp
					*_comp = src.get_comp_at_basis(*this);
				}
			}
		}

		void _move(object<T, DIM, RANK>&& src) {
			if (_basis == nullptr) { // it has no object, so move all
				_comp = std::move(src._comp);
				_basis = std::move(src._basis);
			}
			else { // has object
				if (_basis == src._basis) { // the same => move only comp
					_comp = std::move(src._comp);
				}
				else { // recalc comp and reset src
					*_comp = src.get_comp_at_basis(*this);
					src._comp.reset();
				}
				src._basis.reset(); // reset src object
			}
		}

		void _reset_basis(const Basis<T, DIM>& pbasis) {
			_basis.reset();
			_basis = pbasis;
		}

		object(const container<T, DIM, RANK>& comp, const container<T, DIM, RANK>& basis) {
			_comp = std::make_unique<container<T, DIM, RANK>>(comp);
			_basis = std::make_shared<container<T, DIM, 2>>(basis);
		}
	protected:

		void move_basis(const Basis<T, DIM>& pbasis) {
			_reset_basis(pbasis);
		}
		// R.Rt
		container<T, DIM, 2> get_transform(const Basis<T, DIM>& object) const {
			return mat_scal_mat_transp(*this->_basis, *object);
		}

		container<T, DIM, RANK>& comp() {
			return *this->_comp;
		}
	public:
		bool is_empty() {
			return _comp == nullptr;
		}
		object(const object<T, DIM, RANK>& basis_obj) { // copy ctor
			_copy(basis_obj);
		}

		object(object<T, DIM, RANK>&& basis_obj) noexcept { // move ctor
			_move(std::move(basis_obj));
		}

		object(size_t N, size_t R = 0, FILL_TYPE type = FILL_TYPE::ZERO, const Basis<T, DIM>& pbasis = EMPTY_BASIS<T, DIM>) {
			_comp = std::make_unique<container<T, DIM, RANK>>(container<T, DIM, RANK>(N, R, type));
			if (R > 0) { _basis = pbasis; }
		}

		object(const container<T, DIM, RANK>& comp, const Basis<T, DIM>& pbasis = EMPTY_BASIS<T, DIM>) {
			_comp = std::make_unique<container<T, DIM, RANK>>(comp);
			if (RANK > 0) { _basis = pbasis; }
		}

		object(const container<T, DIM, RANK>& comp, Basis<T, DIM>&& pbasis = EMPTY_BASIS<T, DIM>) {
			_comp = std::make_unique<container<T, DIM, RANK>>(comp);
			if (RANK > 0) {
				std::move(pbasis); 
			} else {
				pbasis.reset();
			}
		}

		object(container<T, DIM, RANK>&& comp, const Basis<T, DIM>& pbasis = EMPTY_BASIS<T, DIM>) noexcept {
			_comp = std::make_unique<container<T, DIM, RANK>>(std::move(comp));
			if (RANK > 0) { _basis = pbasis; }
		}
		
		object(container<T, DIM, RANK>&& comp, Basis<T, DIM>&& pbasis = EMPTY_BASIS<T, DIM>) {
			_comp = std::make_unique<container<T, DIM, RANK>>(std::move(comp));
			if (RANK > 0) {
				std::move(pbasis); 
			} else {
				pbasis.reset();
			}
		}

		void change_basis(const object<T, DIM, RANK>& obj) {
			change_basis(obj.get_basis_ref());
		}

		void change_basis(const Basis<T, DIM>& pbasis) {
			*_comp = get_comp_at_basis(pbasis);
			_reset_basis(pbasis);
		}

		container<T, DIM, RANK> get_comp_at_basis(const object<T, DIM, RANK>& obj) const {
			return get_comp_at_basis(obj.get_basis_ref());
		}

		container<T, DIM, RANK> get_comp_at_basis(const Basis<T, DIM>& pbasis) const {
			const container<T, DIM, RANK>& comp = get_comp_ref();
			if (pbasis == get_basis_ref() || RANK == 0) {
				return comp;
			}
			else {
				container<T, DIM, 2> op = get_transform(pbasis);
				if (RANK == 1) {
					return operator*(comp, op);// comp * op;
				}
				else {
					return op.transpose() * comp * op;
				}
			}
		}

		const Basis<T, DIM> get_basis_ref() const {
			return this->_basis;
		}

		const container<T, DIM, RANK>& get_comp_ref() const {
			return *this->_comp;
		}

		container<T, DIM, RANK> get_basis_comp() const {
			return container<T, DIM, RANK>(*this->_basis);
		}

		container<T, DIM, RANK> get_basis() const {
			return container<T, DIM, RANK>(*this->_basis);
		}

		container<T, DIM, RANK> get_comp() const {
			return container<T, DIM, RANK>(*this->_comp);
		}

		//operator T() const {
		//	if (this->_comp->size() == 1) {
		//		return (*this->_comp)[0];
		//	}
		//	throw ErrorAccess::NoCastScalar();
		//}

		object& operator *= (const T& mul) {
			*this->_comp *= mul;
			return *this;
		}

		object& operator /= (const T& mul) {
			*this->_comp /= mul;
			return *this;
		}

		object& operator *= (const object<T, DIM, RANK>& rhs) {
			*this->_comp *= rhs.get_comp_at_basis(*this);
			return *this;
		}

		object& operator *= (const container<T, DIM, RANK>& rhs) {
			*this->_comp *= rhs;
			return *this;
		}

		object& operator += (const object<T, DIM, RANK>& rhs) {
			*this->_comp += rhs.get_comp_at_basis(*this);
			return *this;
		}

		object& operator += (const container<T, DIM, RANK>& rhs) {
			*this->_comp += rhs;
			return *this;
		}

		object& operator -= (const object<T, DIM, RANK>& rhs) {
			*this->_comp -= rhs.get_comp_at_basis(*this);
			return *this;
		}

		object& operator -= (const container<T, DIM, RANK>& rhs) {
			*this->_comp -= rhs;
			return *this;
		}

		friend bool check_ort(const container<T, DIM, RANK>& m);

		//static friend object<T, DIM, LRANK+RRANK-2> operator * <> (const object<T, DIM, LRANK>& lhs, const object<T, DIM, RRANK>& rhs);


		static friend object<T, DIM, RANK> operator + <> (const object<T, DIM, RANK>& lhs, const object<T, DIM, RANK>& rhs);
		static friend object<T, DIM, RANK> operator - <> (const object<T, DIM, RANK>& lhs, const object<T, DIM, RANK>& rhs);
		static friend object<T, DIM, RANK> operator * <> (const object<T, DIM, RANK>& lhs, const T& mul);
		static friend object<T, DIM, RANK> operator / <> (const object<T, DIM, RANK>& lhs, const T& mul);
		static friend object<T, DIM, RANK> operator * <> (const T& mul, const object<T, DIM, RANK>& rhs);

		object& operator = (const container<T, DIM, RANK>& rhs) { // copy assign
			*_comp = rhs;
			return *this;
		}

		object& operator = (const object<T, DIM, RANK>& rhs) { // copy assign
			_copy(rhs);
			return *this;
		}

		object& operator = (object<T, DIM, RANK>&& rhs) noexcept { // move assign
			_move(std::move(rhs));
			return *this;
		}

		friend bool operator == (const object<T, DIM, RANK>& lhs, const object<T, DIM, RANK>& rhs) {
			return lhs.get_comp_ref() == rhs.get_comp_at_basis(lhs.get_basis_ref());
		}

		friend static object<T, DIM, RANK> transpose(const object<T, DIM, RANK>& m) {
			return object<T, DIM, RANK>(m.get_comp_ref().transpose(), m.get_basis_ref());
		}

		friend static object<T, DIM, RANK> inverse(const object<T, DIM, RANK>& m) {
			return object<T, DIM, RANK>(m.get_comp_ref().inverse(), m.get_basis_ref());
		}

		template<typename T, size_t DIM, size_t RANK>
		friend object<T, DIM, RANK> eigen_object(const tens::container<T, DIM, RANK>& M);

		template<typename T, size_t DIM, size_t RANK>
		friend container<T, DIM, RANK> func(const tens::container<T, DIM, RANK>& M, T(&f)(T));
	};

	template<typename T, size_t DIM, size_t LRANK, size_t RRANK>
	[[nodiscard]] object<T, DIM, LRANK + RRANK - 2> operator * (const object<T, DIM, LRANK>& lhs, const object<T, DIM, RRANK>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs.get_basis_ref());
		return object<T, DIM, LRANK + RRANK - 2>(lhs.get_comp_ref() * rhsa, lhs.get_basis_ref());
	}

	template<typename T, size_t DIM, size_t RANK>
	container<T, DIM, RANK> func(const tens::container<T, DIM, RANK>& M, T(&f)(T)) {
		auto p = eigen(M);
		const size_t size = M.size();
		for (size_t i = 0; i < size; i++)
			p.first[i] = f(p.first[i]);
		auto obj = object<T, DIM, RANK>(p.first, p.second);
		obj.change_basis(GLOBAL_BASIS<T, DIM>);
		return obj.comp();
	}

	template<typename T, size_t DIM, size_t RANK = 2>
	bool check_ort(const container<T, DIM, RANK>& m) {
		const container<T, DIM, RANK> I = m * m.transpose();
		T diag = 0;
		T nondiag = 0;
		for (size_t diagIdx = 0; diagIdx < 3; diagIdx++)
			diag += I[diagIdx];
		for (size_t nonDiagIdx = 3; nonDiagIdx < 9; nonDiagIdx++)
			nondiag += I[nonDiagIdx];
		return (math::is_small_value<T>(abs(diag - (T)DIM) + abs(nondiag)) ? true : false);
	}

	template<typename T, size_t DIM, size_t RANK = 2>
	static container<T, DIM, RANK> create_orthogonal_matrix(DEFAULT_ORTH_BASIS type = DEFAULT_ORTH_BASIS::INDENT) {
		container<T, DIM, RANK> Q;
		switch (type)
		{
		case tens::DEFAULT_ORTH_BASIS::RANDOM:
			Q = generate_rand_ort();
			break;
		default:
		case tens::DEFAULT_ORTH_BASIS::INDENT:
			Q = generate_indent_ort();
			break;
		}
		if (!check_ort(Q)) {
			throw new ErrorMath::NonOrthogonal();
		}
		return Q;
	}

	template<typename T, size_t DIM, size_t RANK = 2>
	static Basis<T, DIM, RANK> create_basis(DEFAULT_ORTH_BASIS type = DEFAULT_ORTH_BASIS::INDENT) {
		return std::make_shared<const container<T, DIM, RANK>>(create_orthogonal_matrix<T, DIM>(type));
	}

	template<typename T, size_t DIM, size_t RANK = 2>
	static object<T, DIM, RANK> create_basis_object(DEFAULT_ORTH_BASIS type = DEFAULT_ORTH_BASIS::INDENT) {
		auto Q = create_orthogonal_matrix<T, DIM>(type);
		return object<T, DIM, RANK>(Q, GLOBAL_BASIS<T>);
	}

	template<typename T, size_t DIM, size_t RANK = 2>
	static Basis<T, DIM, RANK> create_basis(const container<T, DIM, RANK>& object) {
		if (!check_ort(object)) {
			throw ErrorMath::NonOrthogonal();
		}
		return std::make_shared<const container<T, DIM, RANK>>(object);
	}

	template<typename T, size_t DIM, size_t RANK>
	object<T, DIM, RANK> operator * (const object<T, DIM, RANK>& lhs, const T& mul) {
		return object<T, DIM, RANK>(*lhs._comp * mul, lhs._basis);
	}

	template<typename T, size_t DIM, size_t RANK>
	object<T, DIM, RANK> operator / (const object<T, DIM, RANK>& lhs, const T& mul) {
		return object<T, DIM, RANK>(*lhs._comp / mul, lhs._basis);
	}

	template<typename T, size_t DIM, size_t RANK>
	object<T, DIM, RANK> operator * (const T& mul, const object<T, DIM, RANK>& rhs) {
		return object<T, DIM, RANK>(*rhs._comp * mul, rhs._basis);
	}

	template<typename T, size_t DIM, size_t RANK>
	object<T, DIM, RANK> operator + (const object<T, DIM, RANK>& lhs, const object<T, DIM, RANK>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs);
		return object<T, DIM, RANK>(*lhs._comp + rhsa, lhs._basis);
	}

	template<typename T, size_t DIM, size_t RANK>
	object<T, DIM, RANK> operator - (const object<T, DIM, RANK>& lhs, const object<T, DIM, RANK>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs);
		return object<T, DIM, RANK>(*lhs._comp - rhsa, lhs._basis);
	}

	template<typename T, size_t DIM, size_t RANK>
	std::array<T, DIM> get_comp(const object<T, DIM, RANK>& vect) {
		auto arr = vect.get_comp().get();
		std::array<T, DIM> res;
		for (size_t i = 0; i < DIM; ++i) {
			res[i] = arr[i];
		}
		return res;
	}

	template<typename T, size_t DIM, size_t RANK>
	std::array<std::array<T, 3>, 3> get_comp(const object<T, DIM, RANK>& vect) {
		// {00, 11, 22, 12, 02, 01, 21, 20, 10}  
		auto arr = vect.get_comp().get();
		std::array<std::array<T, 3>, 3> res;
		res[0][0] = arr[0];
		res[1][1] = arr[1];
		res[2][2] = arr[2];
		res[1][2] = arr[3];
		res[0][2] = arr[4];
		res[0][1] = arr[5];
		res[2][1] = arr[6];
		res[2][0] = arr[7];
		res[1][0] = arr[8];
		return res;
	}

	template<typename T, size_t DIM, size_t RANK>
	object<T, DIM, RANK> eigen_object(const tens::container<T, DIM, RANK>& M) {
		auto res = eigen(M);
		return object<T, DIM, RANK>(res.first, res.second);
	}

	template<typename T, size_t DIM, size_t RANK>
	std::ostream& operator<<(std::ostream& out, const object<T, DIM, RANK>& b) {
		const auto cont = b.get_comp_ref();
		out << "{ ";
		for (size_t idx = 0; idx < cont.size() - 1; idx++)
			out << cont[idx] << ", ";
		out << cont[cont.size() - 1] << " } \n";
		return out;
	};
};
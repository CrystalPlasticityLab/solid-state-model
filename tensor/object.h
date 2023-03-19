#pragma once
#include "container.h"

namespace tens {

	template<typename T>
	class object;
	template<typename T> std::ostream& operator<< (std::ostream& o, const object<T>& b);

	template<typename T, size_t N>
	object<T> Tensor(const container<T>& m, const std::shared_ptr<const container<T>>& _basis) {
		if (m.rank() != 2) {
			throw new NoImplemetationYet();
		}
		return object<T>(m, _basis);
	}

	template<typename T, size_t N>
	object<T> Tensor(const object<T>& v) {
		if (v.rank() != 2) {
			throw new NoImplemetationYet();
		}
		return object<T>(v);
	}

	template<typename T, size_t N>
	object<T> Tensor(const std::array<std::array<T, N>, N>& m, const std::shared_ptr<const container<T>>& _basis) {
		return object<T>(Matrix<T,N>(m), _basis);
	}

	template<typename T, size_t N>
	object<T> Vector(const std::array<T, N>& _arr, const std::shared_ptr<const container<T>>& _basis) {
		return object<T>(Array<T, N>(_arr), _basis);
	}

	template<typename T, size_t N>
	object<T> Vector(const container<T>& a, const std::shared_ptr<const container<T>>& _basis) {
		if (a.rank() != 1) {
			throw new NoImplemetationYet();
		}
		return object<T>(a, _basis);
	}

	template<typename T, size_t N>
	object<T> Vector(const object<T>& v) {
		if (v.rank() != 1) {
			throw new NoImplemetationYet();
		}
		return object<T>(v);
	}

	enum class DEFAULT_ORTH_BASIS {
		RANDOM,
		INDENT
	};

	template<typename T>
	static object<T> create_basis(size_t N, DEFAULT_ORTH_BASIS type);

	template<typename T>
	static const std::shared_ptr<container<T>> EMPTY_BASIS = std::shared_ptr<container<T>>();

	template<typename T>
	class object {

		std::unique_ptr<container<T>> _comp;
		Basis<T> _basis;

		void _copy(const object<T>& src) {
			if (_basis == nullptr) { // it has no object, so we create it
				_comp = std::make_unique<container<T>>(container<T>(*src._comp));
				if (src._basis) {
					_basis = std::make_shared<container<T>>(container<T>(*src._basis));
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

		void _move(object<T>&& src) {
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

		void _reset_basis(const Basis<T>& pbasis) {
			_basis.reset();
			_basis = pbasis;
		}

	protected:

		void move_basis(const Basis<T>& pbasis) {
			_reset_basis(pbasis);
		}
		// R.Rt
		container<T> get_transform(const Basis<T>& object) const {
			return mat_scal_mat_transp(*this->_basis, *object);
		}

		container<T>& comp() {
			return *this->_comp;
		}
	public:
		bool is_empty() {
			return _comp == nullptr;
		}
		object(const object<T>& basis_obj) { // copy ctor
			_copy(basis_obj);
		}

		object(object<T>&& basis_obj) noexcept { // move ctor
			_move(std::move(basis_obj));
		}

		object(size_t N, size_t R = 0, FILL_TYPE type = FILL_TYPE::ZERO, const Basis<T>& pbasis = EMPTY_BASIS<T>) {
			_comp = std::make_unique<container<T>>(container<T>(N, R, type));
			if (R > 0) { _basis = pbasis; }
		}

		object(const container<T>& comp, const Basis<T>& pbasis = EMPTY_BASIS<T>) {
			_comp = std::make_unique<container<T>>(comp);
			if (comp.rank() > 0) { _basis = pbasis; }
		}

		object(const container<T>& comp, Basis<T>&& pbasis = EMPTY_BASIS<T>) {
			_comp = std::make_unique<container<T>>(comp);
			if (comp.rank() > 0) { 
				std::move(pbasis); 
			} else {
				pbasis.reset();
			}
		}

		object(container<T>&& comp, const Basis<T>& pbasis = EMPTY_BASIS<T>) noexcept {
			_comp = std::make_unique<container<T>>(std::move(comp));
			if (comp.rank() > 0) { _basis = pbasis; }
		}
		
		object(container<T>&& comp, Basis<T>&& pbasis = EMPTY_BASIS<T>) {
			_comp = std::make_unique<container<T>>(std::move(comp));
			if (comp.rank() > 0) { 
				std::move(pbasis); 
			} else {
				pbasis.reset();
			}
		}

		void change_basis(const object<T>& obj) {
			change_basis(obj.get_basis_ref());
		}

		void change_basis(const Basis<T>& pbasis) {
			*_comp = get_comp_at_basis(pbasis);
			_reset_basis(pbasis);
		}

		container<T> get_comp_at_basis(const object<T>& obj) const {
			return get_comp_at_basis(obj.get_basis_ref());
		}

		container<T> get_comp_at_basis(const Basis<T>& pbasis) const {
			const auto& comp = get_comp_ref();
			if (pbasis == get_basis_ref() || comp.rank() == 0) {
				return comp;
			}
			else {
				container<T> op = get_transform(pbasis);
				if (this->_comp->rank() == 1) {
					return comp * op;
				}
				else {
					return op.transpose() * comp * op;
				}
			}
		}

		const Basis<T> get_basis_ref() const {
			return this->_basis;
		}

		const container<T>& get_comp_ref() const {
			return *this->_comp;
		}

		container<T> get_basis_comp() const {
			return container<T>(*this->_basis);
		}

		container<T> get_basis() const {
			return container<T>(*this->_basis);
		}

		container<T> get_comp() const {
			return container<T>(*this->_comp);
		}

		operator T() const {
			if (this->_comp->size() == 1) {
				return (*this->_comp)[0];
			}
			throw ErrorAccess::NoCastScalar();
		}

		object& operator *= (const T& mul) {
			*this->_comp *= mul;
			return *this;
		}

		object& operator /= (const T& mul) {
			*this->_comp /= mul;
			return *this;
		}

		object& operator *= (const object& rhs) {
			*this = *this * rhs;
			return *this;
		}

		object& operator += (const object<T>& rhs) {
			auto rhsa = rhs.get_comp_at_basis(*this);
			*this->_comp += rhsa;
			return *this;
		}

		object& operator += (const container<T>& rhs) {
			*this->_comp += rhs;
			return *this;
		}

		object& operator -= (const object<T>& rhs) {
			auto rhsa = rhs.get_comp_at_basis(*this);
			*this->_comp -= rhsa;
			return *this;
		}

		friend bool check_ort(const container<T>& m);

		static friend object<T> operator * <> (const object<T>& lhs, const object<T>& rhs);
		static friend object<T> operator + <> (const object<T>& lhs, const object<T>& rhs);
		static friend object<T> operator - <> (const object<T>& lhs, const object<T>& rhs);
		static friend object<T> operator * <> (const object<T>& lhs, const T& mul);
		static friend object<T> operator / <> (const object<T>& lhs, const T& mul);
		static friend object<T> operator * <> (const T& mul, const object<T>& rhs);

		object& operator = (const container<T>& rhs) { // copy assign
			*_comp = rhs;
			return *this;
		}

		object& operator = (const object<T>& rhs) { // copy assign
			_copy(rhs);
			return *this;
		}

		object& operator = (object<T>&& rhs) noexcept { // move assign
			_move(std::move(rhs));
			return *this;
		}

		friend bool operator == (const object<T>& lhs, const object<T>& rhs) {
			return lhs.get_comp_ref() == rhs.get_comp_at_basis(lhs.get_basis_ref());
		}

		friend static object<T> transpose(const object<T>& m) {
			return object<T>(m.get_comp_ref().transpose(), m.get_basis_ref());
		}

		friend static object<T> inverse(const object<T>& m) {
			return object<T>(m.get_comp_ref().inverse(), m.get_basis_ref());
		}
	};

	template<typename T>
	bool check_ort(const container<T>& m) {
		if (m.rank() != 2) {
			throw new ErrorMath::ShapeMismatch();
		}
		const container<T> I = m * m.transpose();
		T diag = 0;
		T nondiag = 0;
		for (size_t diagIdx = 0; diagIdx < 3; diagIdx++)
			diag += I[diagIdx];
		for (size_t nonDiagIdx = 3; nonDiagIdx < 9; nonDiagIdx++)
			nondiag += I[nonDiagIdx];
		return (math::is_small_value<T>(abs(diag - (T)m.dim()) + abs(nondiag)) ? true : false);
	}

	template<typename T, size_t N>
	static container<T> create_orthogonal_matrix(DEFAULT_ORTH_BASIS type = DEFAULT_ORTH_BASIS::INDENT) {
		container<T> Q(N, 2);
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

	template<typename T, size_t N>
	static Basis<T> create_basis(DEFAULT_ORTH_BASIS type = DEFAULT_ORTH_BASIS::INDENT) {
		return std::make_shared<const container<T>>(create_orthogonal_matrix<T, N>(type));
	}

	template<typename T, size_t N>
	static object<T> create_basis_object(DEFAULT_ORTH_BASIS type = DEFAULT_ORTH_BASIS::INDENT) {
		auto Q = create_orthogonal_matrix<T, N>(type);
		return object<T>(Q, GLOBAL_BASIS<T>);
	}

	template<typename T, size_t N>
	static Basis<T> create_basis(const container<T>& object) {
		if (!check_ort(object)) {
			throw ErrorMath::NonOrthogonal();
		}
		return std::make_shared<const container<T>>(object);
	}

	template<typename T>
	object<T> operator * (const object<T>& lhs, const T& mul) {
		return object<T>(*lhs._comp * mul, lhs._basis);
	}

	template<typename T>
	object<T> operator / (const object<T>& lhs, const T& mul) {
		return object<T>(*lhs._comp / mul, lhs._basis);
	}

	template<typename T>
	object<T> operator * (const T& mul, const object<T>& rhs) {
		return object<T>(*rhs._comp * mul, rhs._basis);
	}

	template<typename T>
	object<T> operator + (const object<T>& lhs, const object<T>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs);
		return object<T>(*lhs._comp + rhsa, lhs._basis);
	}

	template<typename T>
	object<T> operator - (const object<T>& lhs, const object<T>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs);
		return object<T>(*lhs._comp - rhsa, lhs._basis);
	}

	template<typename T>
	object<T> operator * (const object<T>& lhs, const object<T>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs.get_basis_ref());
		return object<T>(lhs.get_comp_ref() * rhsa, lhs.get_basis_ref());
	}

	template<typename T, size_t N>
	std::array<T, N> get_comp(const object<T>& vect) {
		auto arr = vect.get_comp().get();
		std::array<T, N> res;
		for (size_t i = 0; i < N; ++i) {
			res[i] = arr[i];
		}
		return res;
	}

	template<typename T>
	std::array<std::array<T, 3>, 3> get_comp(const object<T>& vect) {
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

	template<typename T>
	std::ostream& operator<<(std::ostream& out, const object<T>& b) {
		const auto cont = b.get_comp_ref();
		out << "{ ";
		for (size_t idx = 0; idx < cont.size() - 1; idx++)
			out << cont[idx] << ", ";
		out << cont[cont.size() - 1] << " }";
		return out;
	};
};
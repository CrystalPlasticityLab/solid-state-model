#pragma once
#include "container.h"

namespace tens {

	template<typename T, size_t N, size_t R>
	class basis;
	template<typename T, size_t N, size_t R> std::ostream& operator<< (std::ostream& o, const basis<T, N, R>& b);

	template<typename T, size_t N>
	static Basis<T, N> create_basis(const Basis<T, N>& basis);


	enum class DEFAULT_ORTH_BASIS {
		RANDOM,
		INDENT
	};

	template<typename T, size_t N>
	static Basis<T, N> create_basis(DEFAULT_ORTH_BASIS type);

	template<typename T, size_t N>
	class basis_base { // just for abstraction, excluding R as parameter

		virtual size_t get_rank() = 0;
	};

	template<typename T, size_t N, size_t R>
	class basis : public basis_base<T,N> {

		std::unique_ptr<container<T, N, R>> _comp;
		Basis<T, N> _basis;

		void _copy(const basis<T, N, R>& src) {
			if (_basis == nullptr) { // it has no basis, so we create it
				_comp = std::make_unique<container<T, N, R>>(container<T, N, R>(*src._comp));
				_basis = std::make_shared<basis_cont<T, N>>(basis_cont<T, N>(*src._basis));
			}
			else { // has basis, recalc comp at this basis if it is necessery
				if (_basis == src._basis) { // basis is the same, just copy comp
					*_comp = *src._comp;
				}
				else { // other basis - recalc comp
					*_comp = src.get_comp_at_basis(*this);
				}
			}
		}

		void _move(basis<T, N, R>&& src) {
			if (_basis == nullptr) { // it has no basis, so move all
				_comp = std::move(src._comp);
				_basis = std::move(src._basis);
			}
			else { // has basis
				if (_basis == src._basis) { // the same => move only comp
					_comp = std::move(src._comp);
				}
				else { // recalc comp and reset src
					*_comp = src.get_comp_at_basis(*this);
					src._comp.reset();
				}
				src._basis.reset(); // reset src basis
			}
		}

		void _reset_basis(const Basis<T, N>& pbasis) {
			_basis.reset();
			_basis = pbasis;
		}

	protected:

		void move_basis(const std::shared_ptr<basis_cont<T, N>>& pbasis) {
			_reset_basis(pbasis);
		}
		// R.Rt
		basis_cont<T, N> get_transform(const Basis<T, N>& basis) const {
			return mat_scal_mat_transp(*this->_basis, *basis);
		}
	public:
		basis(const basis<T, N, R>& basis_obj) { // copy ctor
			_copy(basis_obj);
		}

		basis(basis<T, N, R>&& basis_obj) noexcept { // move ctor
			_move(std::move(basis_obj));
		}

		basis(const container<T, N, R>& comp, const basis_cont<T, N>& basis) {
			_comp = std::make_unique<container<T, N, R>>(comp);
			_basis = create_basis(basis);
		}

		basis(const std::array<T, N>& comp, const basis_cont<T, N>& basis) {
			_comp = std::make_unique<container<T, N, R>>(comp);
			_basis = create_basis(basis);
		}

		basis(container<T, N, R>&& comp, const basis_cont<T, N>& basis) noexcept {
			_comp = std::make_unique<container<T, N, R>>(std::move(comp));
			_basis = create_basis(basis);
		}

		basis(const container<T, N, R>& comp, basis_cont<T, N>&& basis) noexcept {
			_comp = std::make_unique<container<T, N, R>>(comp);
			_basis = std::move(basis);
		}

		basis(container<T, N, R>&& comp, basis_cont<T, N>&& basis) noexcept {
			_comp = std::make_unique<container<T, N, R>>(std::move(comp));
			_basis = std::move(basis);
		}

		basis(const std::array<T, N>& comp, const Basis<T, N>& pbasis) {
			_comp = std::make_unique<container<T, N, R>>(comp);
			_basis = pbasis;
		}

		basis(const container<T, N, R>& comp, const Basis<T, N>& pbasis) {
			_comp = std::make_unique<container<T, N, R>>(comp);
			_basis = pbasis;
		}

		basis(container<T, N, R>&& comp, const Basis<T, N>& pbasis) noexcept {
			_comp = std::make_unique<container<T, N, R>>(std::move(comp));
			_basis = pbasis;
		}

		void change_basis(const basis<T, N, R>& obj) {
			change_basis(obj.get_basis_ref());
		}

		void change_basis(const Basis<T, N>& pbasis) {
			_comp = std::make_unique<container<T, N, R>>(std::move(get_comp_at_basis(pbasis)));
			_reset_basis(pbasis);
		}

		size_t get_rank() override { return R;};

		container<T, N, R> get_comp_at_basis(const basis<T, N, R>& obj) const {
			return get_comp_at_basis(obj.get_basis_ref());
		}

		container<T, N, R> get_comp_at_basis(const Basis<T, N>& pbasis) const {
			const auto& comp = get_comp_ref();
			if (pbasis == get_basis_ref()) {
				return comp;
			}
			else {
				basis_cont<T, N> op = get_transform(pbasis);
				if (R == 1) {
					return comp * op;
				}
				else {
					return transpose(op) * comp * op;
				}
			}
		}

		const Basis<T, N> get_basis_ref() const {
			return this->_basis;
		}

		const container<T, N, R>& get_comp_ref() const {
			return *this->_comp;
		}

		basis_cont<T, N> get_basis_comp() const {
			return basis_cont<T, N>(*this->_basis);
		}

		container<T, N, R> get_comp() const {
			return container<T, N, R>(*this->_comp);
		}

		basis operator *= (const T& mul) {
			*this->_comp *= mul;
			return *this;
		}

		basis operator *= (const basis& rhs) {
			*this = *this * rhs;
			return *this;
		}

		basis operator += (const basis<T, N, R>& rhs) {
			auto rhsa = rhs.get_comp_at_basis(*this);
			*this->_comp += rhsa;
			return *this;
		}

		basis operator -= (const basis<T, N, R>& rhs) {
			auto rhsa = rhs.get_comp_at_basis(*this);
			*this->_comp -= rhsa;
			return *this;
		}

		static friend bool check_ort(const basis_cont<T, N>& m);

		static friend T operator * <> (const basis<T, N, 1>& lhs, const  basis<T, N, 1>& rhs);
		static friend basis<T, N, R> operator * <> (const basis<T, N, R>& lhs, const T& mul);
		static friend basis<T, N, R> operator * <> (const T& mul, const basis<T, N, R>& rhs);
		static friend basis<T, N, R> operator + <> (const basis<T, N, R>& lhs, const  basis<T, N, R>& rhs);
		static friend basis<T, N, R> operator - <> (const basis<T, N, R>& lhs, const  basis<T, N, R>& rhs);

		basis& operator = (const basis<T, N, R>& rhs) { // copy assign
			_copy(rhs);
			return *this;
		}

		basis& operator = (basis<T, N, R>&& rhs) noexcept { // move assign
			_move(std::move(rhs));
			return *this;
		}

		friend bool operator == (const basis<T, N, R>& lhs, const basis<T, N, R>& rhs) {
			return lhs.get_comp_ref() == rhs.get_comp_at_basis(lhs.get_basis_ref());
		}

		friend static basis<T, N, R> transpose(const basis<T, N, R>& m) {
			return basis<T, N, R>(transpose(m.get_comp_ref()), m.get_basis_ref());
		}
	};

	template<typename T, size_t N, size_t R = 2>
	static bool check_ort(const container<T, N, 2>& m) {
		container<T, N, 2> I = m * transpose(m);
		T diag = 0;
		T nondiag = 0;
		for (size_t diagIdx = 0; diagIdx < 3; diagIdx++)
			diag += I[diagIdx];
		for (size_t nonDiagIdx = 3; nonDiagIdx < 9; nonDiagIdx++)
			nondiag += I[nonDiagIdx];
		return (is_small_value<T>(abs(diag - (T)N) + abs(nondiag)) ? true : false);
	}

	template<typename T, size_t N>
	static Basis<T, N> create_basis(DEFAULT_ORTH_BASIS type = DEFAULT_ORTH_BASIS::INDENT) {
		basis_cont<T, N> Q;
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
		return std::make_shared<const basis_cont<T, N>>(Q);
	}

	template<typename T, size_t N>
	static Basis<T, N> create_basis(const basis_cont<T, N>& basis) {
		if (!check_ort(basis)) {
			throw ErrorMath::NonOrthogonal();
		}
		return std::make_shared<const basis_cont<T, N>>(basis);
	}

	template<typename T, size_t N, size_t R>
	basis<T, N, R> operator * (const basis<T, N, R>& lhs, const T& mul) {
		return basis<T, N, R>(*lhs._comp * mul, lhs._basis);
	}

	template<typename T, size_t N, size_t R>
	basis<T, N, R> operator * (const T& mul, const basis<T, N, R>& rhs) {
		return basis<T, N, R>(*rhs._comp * mul, rhs._basis);
	}

	template<typename T, size_t N, size_t R>
	basis<T, N, R> operator + (const basis<T, N, R>& lhs, const basis<T, N, R>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs);
		return basis<T, N, R>(*lhs._comp + rhsa, lhs._basis);
	}

	template<typename T, size_t N, size_t R>
	basis<T, N, R> operator - (const basis<T, N, R>& lhs, const basis<T, N, R>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs);
		return basis<T, N, R>(*lhs._comp - rhsa, lhs._basis);
	}

	template<typename T, size_t N, size_t Rl, size_t Rr >
	basis<T, N, Rl + Rr - 2> operator * (const basis<T, N, Rl>& lhs, const basis<T, N, Rr>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs.get_basis_ref());
		return basis<T, N, Rl + Rr - 2>(lhs.get_comp_ref() * rhsa, lhs.get_basis_ref());
	}

	template<typename T, size_t N, size_t R = 1>
	T operator * (const basis<T, N, 1>& lhs, const basis<T, N, 1>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs);
		return *lhs._comp * rhsa;
	}

	template <typename T, size_t N>
	using Vector = tens::basis<T, N, 1>;
	template <typename T, size_t N>
	using Tensor = tens::basis<T, N, 2>;


	template<typename T, size_t N>
	std::array<T, N> get_comp(const basis<T, N, 1>& vect) {
		auto arr = vect.get_comp().get();
		std::array<T, N> res;
		for (size_t i = 0; i < N; ++i) {
			res[i] = arr[i];
		}
		return res;
	}

	template<typename T>
	std::array<std::array<T, 3>, 3> get_comp(const basis<T, 3, 2>& vect) {
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

	template<typename T, size_t N, size_t R>
	std::ostream& operator<<(std::ostream& out, const basis<T, N, R>& b) {
		const auto cont = b.get_comp_ref();
		out << "{ ";
		for (size_t idx = 0; idx < cont._size - 1; idx++)
			out << cont[idx] << ", ";
		out << cont[cont._size - 1] << " }\n";
		return out;
	};
};
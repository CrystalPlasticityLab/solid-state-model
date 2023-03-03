#pragma once
#include "matrix.h"
#include "container.h"

namespace tens {

	template<typename T, size_t N>
	static std::shared_ptr<const Basis<T, N>> create_basis(const Basis<T, N>& basis);

	template<typename T, size_t N>
	static std::shared_ptr<const Basis<T, N>> create_basis();

	template<typename T, size_t N, size_t R>
	class basis_object {
		std::unique_ptr<container_rank<T, N, R>> _comp;
		std::shared_ptr<const Basis<T, N>> _basis;

		void _copy(const basis_object<T, N, R>& src) {
			if (_basis == nullptr) { // it has no basis, so we create it
				_comp = std::make_unique<container_rank<T, N, R>>(container_rank<T, N, R>(*src._comp));
				_basis = std::make_shared<Basis<T, N>>(Basis<T, N>(*src._basis));
			}
			else { // has basis, recalc comp at this basis if it is necessery
				if (_basis == src._basis) { // basis is the same, just copy comp
					*_comp = *src._comp;
				}
				else { // other basis - recalc comp
					*_comp = src.get_comp_at_basis(_basis);
				}
			}
		}

		void _move(basis_object<T, N, R>&& src) {
			if (_basis == nullptr) { // it has no basis, so move all
				_comp = std::move(src._comp);
				_basis = std::move(src._basis);
			}
			else { // has basis
				if (_basis == src._basis) { // the same => move
					_comp = std::move(src._comp);
				}
				else { // recalc comp and reset src
					*_comp = src.get_comp_at_basis(_basis);
					src._comp.reset();
				}
				src._basis.reset(); // reset src basis
			}
		}

		void _reset_basis(const std::shared_ptr<const Basis<T, N>>& pbasis) {
			_basis.reset();
			_basis = pbasis;
		}

	protected:

		void change_basis(const std::shared_ptr<Basis<T, N>>& pbasis) {
			_comp = this->get_comp_at_basis(pbasis);
			_reset_basis(pbasis);
		}
		void move_basis(const std::shared_ptr<Basis<T, N>>& pbasis) {
			_reset_basis(pbasis);
		}

		Basis<T, N> get_transform(const std::shared_ptr<const Basis<T, N>>& basis) const {
			return mat_scal_mat_transp(*this->_basis, *basis);
		}
	public:
		basis_object(const basis_object<T, N, R>& basis_obj) { // copy ctor
			_copy(basis_obj);
		}

		basis_object(basis_object<T, N, R>&& basis_obj) noexcept { // move ctor
			_move(std::move(basis_obj));
		}

		basis_object(const container_rank<T, N, R>& comp, const Basis<T, N>& basis) {
			_comp = std::make_unique<container_rank<T, N, R>>(comp);
			_basis = create_basis(basis);
		}

		basis_object(container_rank<T, N, R>&& comp, const Basis<T, N>& basis) noexcept {
			_comp = std::make_unique<container_rank<T, N, R>>(std::move(comp));
			_basis = create_basis(basis);
		}

		basis_object(const container_rank<T, N, R>& comp, Basis<T, N>&& basis) noexcept {
			_comp = std::make_unique<container_rank<T, N, R>>(comp);
			_basis = std::move(basis);
		}

		basis_object(container_rank<T, N, R>&& comp, Basis<T, N>&& basis) noexcept {
			_comp = std::make_unique<container_rank<T, N, R>>(std::move(comp));
			_basis = std::move(basis);
		}

		basis_object(const container_rank<T, N, R>& comp, const std::shared_ptr<const Basis<T, N>>& pbasis) {
			_comp = std::make_unique<container_rank<T, N, R>>(comp);
			_basis = pbasis;
		}

		basis_object(container_rank<T, N, R>&& comp, const std::shared_ptr<const Basis<T, N>>& pbasis) noexcept {
			_comp = std::make_unique<container_rank<T, N, R>>(std::move(comp));
			_basis = pbasis;
		}

		size_t get_rank() {	return R;};

		container_rank<T, N, R> get_comp_at_basis(const std::shared_ptr<const Basis<T, N>>& pbasis) const {
			const auto& comp = this->get_comp_ref();
			if (pbasis == this->get_basis_ref()) {
				return comp;
			}
			else {
				Basis<T, N> op = this->get_transform(pbasis);
				if (R == 1) {
					return comp * op;
				}
				else {
					return transpose(op) * comp * op;
				}
			}
		}
		const std::shared_ptr<const Basis<T, N>> get_basis_ref() const {
			return this->_basis;
		}

		const container_rank<T, N, R>& get_comp_ref() const {
			return *this->_comp;
		}

		Basis<T, N> get_basis_comp() const {
			return Basis<T, N>(*this->_basis);
		}

		container_rank<T, N, R> get_comp() const {
			return container_rank<T, N, R>(*this->_comp);
		}

		basis_object operator *= (const T& mul) {
			*this->_comp *= mul;
			return *this;
		}

		basis_object operator += (const basis_object<T, N, R>& rhs) {
			auto rhsa = rhs.get_comp_at_basis(this->_basis);
			*this->_comp += rhsa;
			return *this;
		}

		static bool check_ort(const Basis<T, N>& m) {
			container_rank<T, N, 2> I = m * transpose(m);
			T diag = 0;
			T nondiag = 0;
			//for (size_t row = 0; row < N; row++)
			//	for (size_t col = 0; col < N; col++)
			//		(row == col) ? diag += m[row][row] : nondiag += m[row][col];
			return (is_small_value<T>(abs(diag - (T)N) + abs(nondiag)) ? true : false);
		}

		static friend T operator * <> (const basis_object<T, N, 1>& lhs, const  basis_object<T, N, 1>& rhs);
		static friend basis_object<T, N, R> operator * <> (const basis_object<T, N, R>& lhs, const T& mul);
		static friend basis_object<T, N, R> operator * <> (const T& mul, const basis_object<T, N, R>& rhs);
		static friend basis_object<T, N, R> operator + <> (const basis_object<T, N, R>& lhs, const  basis_object<T, N, R>& rhs);
		static friend basis_object<T, N, R> operator - <> (const basis_object<T, N, R>& lhs, const  basis_object<T, N, R>& rhs);

		basis_object& operator = (const basis_object<T, N, R>& rhs) { // copy assign
			_copy(rhs);
			return *this;
		}

		basis_object& operator = (basis_object<T, N, R>&& rhs) noexcept { // move assign
			_move(std::move(rhs));
			return *this;
		}
	};

	template<typename T, size_t N>
	static std::shared_ptr<const Basis<T, N>> create_basis() {
		Basis<T, N> I;
		//check_ort(I);
		return std::make_shared<const Basis<T, N>>(I);
	}

	template<typename T, size_t N>
	static std::shared_ptr<const Basis<T, N>> create_basis(const Basis<T, N>& basis) {
		//check_ort(I);
		return std::make_shared<const Basis<T, N>>(basis);
	}

	template<typename T, size_t N, size_t R>
	basis_object<T, N, R> operator * (const basis_object<T, N, R>& lhs, const T& mul) {
		return basis_object<T, N, R>(*lhs._comp * mul, lhs._basis);
	}

	template<typename T, size_t N, size_t R>
	basis_object<T, N, R> operator * (const T& mul, const basis_object<T, N, R>& rhs) {
		return basis_object<T, N, R>(*rhs._comp * mul, rhs._basis);
	}

	template<typename T, size_t N, size_t R>
	basis_object<T, N, R> operator + (const basis_object<T, N, R>& lhs, const basis_object<T, N, R>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs._basis);
		return basis_object<T, N, R>(*lhs._comp + rhsa, lhs._basis);
	}

	template<typename T, size_t N, size_t R>
	basis_object<T, N, R> operator - (const basis_object<T, N, R>& lhs, const basis_object<T, N, R>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs._basis);
		return basis_object<T, N, R>(*lhs._comp - rhsa, lhs._basis);
	}

	template<typename T, size_t N, size_t Rl, size_t Rr >
	basis_object<T, N, Rl + Rr - 2> operator * (const basis_object<T, N, Rl>& lhs, const basis_object<T, N, Rr>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs.get_basis_ref());
		return basis_object<T, N, Rl + Rr - 2>(lhs.get_comp_ref() * rhsa, lhs.get_basis_ref());
	}

	template<typename T, size_t N, size_t R = 1>
	T operator * (const basis_object<T, N, 1>& lhs, const basis_object<T, N, 1>& rhs) {
		auto rhsa = rhs.get_comp_at_basis(lhs._basis);
		return *lhs._comp * rhsa;
	}


	template <typename T, size_t N>
	using Vector = tens::basis_object<T, N, 1>;
	template <typename T, size_t N>
	using Tensor = tens::basis_object<T, N, 2>;


	template<typename T, size_t N>
	std::array<T, N> get_comp(const basis_object<T, N, 1>& vect) {
		auto arr = vect.get_comp().get();
		std::array<T, N> res;
		for (size_t i = 0; i < N; ++i) {
			res[i] = arr[i];
		}
		return res;
	}

	template<typename T>
	std::array<std::array<T, 3>, 3> get_comp(const basis_object<T, 3, 2>& vect) {
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




	template<typename T, size_t N>
	class object : private container<T,N>
	{
		typedef  matrix<T, N> M;
		typedef  std::shared_ptr <M> _shared;
		typedef  container<T, N> _cont;

		std::shared_ptr<matrix<T, N>> basis;
	public:
		object(const _cont& c, const M& m) : basis(std::make_shared<M>(m)), _cont(c) {

		}
	};

	template<typename T, size_t N>
	class basis : public std::shared_ptr<matrix<T, N>> {
		typedef  matrix<T, N> M;
		typedef  std::shared_ptr <M> _shared;

		const basis& operator()() { return *this; };

		static const matrix<T, N> GLOBAL_DEFAULT_BASIS;
		bool owner = false;
		//size_t state = (size_t)this;
	protected:
		//size_t get_state() { return (size_t)state; };
		basis() : _shared() {};
		basis(const basis&  sh) : _shared(sh) {}; // to prevent multiply owning out of scope of vector/tensor
		basis(basis&& sh) noexcept  { _move(static_cast<basis&&>(sh)); };
		basis& operator = (const basis& sh) { _shared::operator=(sh); return *this; };
		basis& operator = (basis&& sh) noexcept { return _move(static_cast<basis&&>(sh));};
		basis& _move(basis&& rhs) {
			_shared::reset();
			static_cast<_shared&>(*this) = std::move(static_cast<_shared&&>(rhs));
			this->owner = rhs.owner;
			return *this;
		};
		void _deep_copy(const basis& rhs) {
			_shared::reset();
			static_cast<_shared&>(*this) = std::make_shared<M>(M(*rhs.get()));
			this->owner = true;
		}
		const void     set_basis(const basis& rhs) { *this = (rhs); }

		explicit basis(const matrix<T, N>& m) : _shared(std::make_shared<matrix<T, N>>(m)) {
			if (!check_ort_basis()){
				throw ErrorMath::NonOrthogonal();
			}
			owner = true;
		}
		/* WARNING : CHANGE Object: move to the target basis m, just basis changes*/
		virtual void    move_to_basis(const basis& m) { set_basis(m); };
		/* NOTE: not change object: just recalc components*/
		virtual void    change_basis(const basis& m) { set_basis(m); }; // must be overrided in chldren classes
	public:
		virtual size_t get_rank() const {return 0;};

		static basis<T,N> create_basis(const matrix<T, N>& m){
			return basis<T,N>(m);
		}

		static basis<T, N> create_basis_random() {
			return basis<T, N>(generate_rand_ort<T, N>());
		}

		bool check_ort_basis() const {
			return this->get()->check_ort();
		}

		friend std::ostream& operator<< (std::ostream& out, const basis& t) { out << *(t.get()); return out; };
	};

	template<typename T, size_t N>
	const matrix<T, N> GLOBAL_DEFAULT_BASIS = matrix<T, N>(MATRIXINITTYPE::INDENT);
};
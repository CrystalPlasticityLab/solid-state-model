#pragma once
#include "matrix.h"

namespace tens {

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
		size_t state = (size_t)this;
	protected:
		size_t get_state() { return (size_t)state; };
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
		virtual size_t get_rank() const  = 0;

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
#pragma once
#include "matrix.h"

namespace tens {
	template<typename T, size_t N>
	class shared_handler_basis : public std::shared_ptr<matrix<T, N>> {
		typedef  matrix<T, N> M;
		typedef  std::shared_ptr <M> _shared;

		const shared_handler_basis& operator()() { return *this; };

		static const shared_handler_basis<T, N> GLOBAL_DEFAULT_BASIS;
		bool owner = false;
	protected:
		shared_handler_basis() : _shared() {};
		shared_handler_basis(const shared_handler_basis&  sh) : _shared(sh) {}; // to prevent multiply owning out of scope of vector/tensor
		shared_handler_basis(shared_handler_basis&& sh) noexcept  { _move(static_cast<shared_handler_basis&&>(sh)); };
		shared_handler_basis& operator = (const shared_handler_basis& sh) { _shared::operator=(sh); return *this; };
		shared_handler_basis& operator = (shared_handler_basis&& sh) noexcept { return _move(static_cast<shared_handler_basis&&>(sh));};
		shared_handler_basis& _move(shared_handler_basis&& rhs) {
			_shared::reset();
			static_cast<_shared&>(*this) = std::move(static_cast<_shared&&>(rhs));
			this->owner = rhs.owner;
			return *this;
		};
		void _deep_copy(const shared_handler_basis& rhs) {
			_shared::reset();
			static_cast<_shared&>(*this) = std::make_shared<M>(M(*rhs.get()));
			this->owner = true;
		}
		const void     set_basis(const shared_handler_basis& rhs) { *this = (rhs); }

		explicit shared_handler_basis(const matrix<T, N>& m) : _shared(std::make_shared<matrix<T, N>>(m)) {
			if (!check_ort_basis()){
				throw ErrorMath::NonOrthogonal();
			}
			owner = true;
		}
		/* WARNING : CHANGE Object: move to the target basis m, just basis changes*/
		virtual void    move_to_basis(const shared_handler_basis& m) { set_basis(m); };
		/* NOTE: not change object: just recalc components*/
		virtual void    change_basis(const shared_handler_basis& m) { set_basis(m); }; // must be overrided in chldren classes
	public:
		static shared_handler_basis<T,N> create_basis(const matrix<T, N>& m){
			return shared_handler_basis<T,N>(m);
		}

		bool check_ort_basis() const {
			return this->get()->check_ort();
		}

		friend std::ostream& operator<< (std::ostream& out, const shared_handler_basis& t) { out << *(t.get()); return out; };
	};
	
	template<typename T, size_t N>
	shared_handler_basis<T,N> create_basis(const matrix<T, N>& m){
		return shared_handler_basis<T,N>::create_basis(m);
	}

	template<typename T, size_t N>
	const shared_handler_basis<T, N> GLOBAL_DEFAULT_BASIS = create_basis<T,N>(matrix<T, N>(MATRIXINITTYPE::INDENT));
};
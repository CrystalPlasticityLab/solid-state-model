#pragma once
#include "../tensor-matrix/tensor/object.h"
#include <unordered_map>

namespace state {
	template<typename T>
	class State;
};

namespace measure {
	using namespace state;


	class ValueUpdated : public std::exception {
	public:
		virtual const char* what() const noexcept {
			return "Measure::value was updated twice";
		};
	};

	class RateUpdated : public std::exception {
	public:
		virtual const char* what() const noexcept {
			return "Measure::rate was updated twice";
		};
	};

	class StateNotLinked : public std::exception {
	public:
		virtual const char* what() const noexcept {
			return "Measure::state not linked or alredy has been destroyed";
		};
	};

	template<typename T>
	class Measure {
		std::weak_ptr<State<T>> _state;
		std::string _name;
		typedef std::shared_ptr<tens::object<T>> object_ptr;
		bool _value_updated = false;
		bool _rate_updated = false;
		const T& _dt;
		tens::object<T> _rate;
		tens::object<T> _value;
		tens::object<T> _value_prev;
	public:
		T dt() { return _dt; };
		const tens::object<T>& rate() const { return _rate; };
		const tens::object<T>& value() const { return _value; };
		const tens::object<T>& value_prev() const { return _value_prev; };
		std::string name() { return _name; };
		// -- TODO: add move semantic 
		Measure(const tens::object<T>& obj, std::string name) :
			//_state(state),
			_name(name),
			_value(obj),
			_value_prev(obj),
			_dt(state->dt()),
			_rate(obj.get_comp_ref().dim(), obj.get_comp_ref().rank(), tens::FILL_TYPE::ZERO, obj.get_basis_ref()) {
			this->change_basis(state->basis());
			state->insert(std::unique_ptr<Measure<T>>(this));
		};

		Measure(tens::object<T>&& obj, std::string name) :
			//_state(state),
			_name(name),
			_value(obj),
			_value_prev(std::move(obj)),
			_dt(state->dt()),
			_rate(obj.get_comp_ref().dim(), obj.get_comp_ref().rank(), tens::FILL_TYPE::ZERO, obj.get_basis_ref()) {
			this->change_basis(state->basis());
			state->insert(std::unique_ptr<Measure<T>>(this));
		};

		template <typename P>
		void update_rate(P&& rate) {
			_rate_updated ? throw new RateUpdated() : _rate = std::forward<P>(rate);
			_rate_updated = true;
		};

		template <typename P>
		void update_value(P&& value) {
			_value_updated ? throw new ValueUpdated() : _value = std::forward<P>(value);
			_value_updated = true;
		};

		virtual void calc_rate() {
			_rate_updated ? throw new RateUpdated() : _rate = (_value - _value_prev) / _dt;
			_rate_updated = true;
		};

		virtual void integrate_value() {
			_value_updated ? throw new ValueUpdated() : _value = _value_prev + _rate *_dt;
			_value_updated = true;
		};

		void finalize() {
			std::swap(_value, _value_prev);
			_value_updated = false;
			_rate_updated = false;
		};

		// access by const ref to other Measures in the State
		const Measure<T>& operator[] (std::string name) {
			if (std::shared_ptr<State<T>> state = this->_state.lock()) {
				return *(*state.get())[name];
			}
			throw new StateNotLinked();
		}

		void change_basis(const Basis<T>& basis) {
			_rate.change_basis(basis);
			_value.change_basis(basis);
			_value_prev.change_basis(basis);
		}
	};

	namespace MEASURE {
		namespace STRAIN {
			const std::string DEFORM_GRADIENT = "F";
			const std::string CAUCHY_GREEN = "C";
		};
		namespace STRESS {
			const std::string CAUCHY = "S";
		}
	};

	template<typename T>
	class GradDeform : public Measure<T> {
	public:
		GradDeform() : Measure<T>(tens::object<T>(IDENT_MATRIX<T>, state->basis()), MEASURE::STRAIN::DEFORM_GRADIENT) {};

		virtual void integrate_value() override {
			auto L = this->rate();
			L *= -this->dt(); // -L*dt
			L += IDENT_MATRIX<T>; // I - L*dt
			Measure<T>::update_value(inverse(L) * this->value_prev());// (I - L * dt)^-1 * F
		}

		std::pair<tens::object<T>, tens::object<T>> polar_decomposition() {
			const auto& F = this->value();
			auto V = F * transpose(F); // TODO: take sqrt !!!
			auto R = F * inverse(V);

			return { V, R };
		}
	};

	template<typename T>
	class CaushyStress : public Measure<T> {
	public:
		CaushyStress() : Measure<T>(tens::object<T>(IDENT_MATRIX<T>, state->basis()), MEASURE::STRESS::CAUCHY) {};
	};
};

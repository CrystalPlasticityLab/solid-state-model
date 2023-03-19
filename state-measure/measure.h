#pragma once
#include "../tensor-matrix/tensor/object.h"

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
	class Measure : public tens::object<T> {
		std::string _name;
		tens::container<T>& _value;
		tens::container<T> _rate;
		tens::container<T> _value_prev;
		tens::container<T> _rate_prev;
		std::weak_ptr<State<T>> _state;
		T& _dt;
		typedef std::shared_ptr<tens::object<T>> object_ptr;
		bool _value_updated = false;
		bool _rate_updated = false;
	public:
		const tens::container<T>& rate() const { return _rate; };
		const tens::container<T>& value() const { return this->get_comp_ref(); };
		const tens::container<T>& value_prev() const { return _value_prev; };
		T dt() { return _dt; };
		std::string name() { return _name; };
		// -- TODO: add move semantic 
		Measure(std::shared_ptr<State<T>>& state, const tens::object<T>& obj, std::string name) : 
			_state(state),
			_dt(state->dt()),
			_name(name),
			_value(static_cast<tens::container<T>&>(this->comp())),
			_rate(obj.get_comp_ref().dim(), obj.get_comp_ref().rank(), tens::FILL_TYPE::ZERO),
			_value_prev(obj.get_comp_ref()),
			_rate_prev(_rate),
			tens::object<T>(obj) {
			this->change_basis(state->basis());
			insert(state, *this);
		};

		template <typename P>
		void update_rate(P&& rate) {
			if (_rate_updated) {
				throw new ValueUpdated();
			}
			else {
				std::swap(_rate, _rate_prev);
				_rate = std::forward<P>(rate);
			}
			_rate_updated = true;
		};

		template <typename P>
		void update_value(P&& value) {
			if (_value_updated) {
				throw new ValueUpdated();
			} else {
				std::swap(_value, _value_prev);
				_value = std::forward<P>(value);
			}
			_value_updated = true;
		};

		virtual void calc_rate(){
			update_rate((_value - _value_prev) / _dt);
		};

		virtual void integrate_value() {
			update_value(_value_prev + _rate * _dt);
		};

		virtual void finalize() {
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

		void set_state(std::shared_ptr<State<T>> state) {
			_state = state;
		}
	};
};

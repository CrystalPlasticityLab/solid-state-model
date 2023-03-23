#pragma once
#include "../tensor-matrix/tensor/object.h"

namespace state {
	template<typename T>
	class State;
};

namespace measure {
	using namespace state;

	namespace error {
		class StateNotLinked : public std::exception {
		public:
			virtual const char* what() const noexcept {
				return "Measure::state not linked or alredy has been destroyed";
			};
		};
	}

	template<typename Base>
	class MeasureAbstract : public Base {
	protected:
		Base& _value;
		Base _rate;
		
		template <typename P>
		void update_rate(P&& rate) {
			_rate = std::forward<P>(rate);
		};

		template <typename P>
		void update_value(P&& value) {
			_value = std::forward<P>(value);
		};

	public:
		MeasureAbstract(Base& value) : Base(value),
			_value(static_cast<Base&>(*this)),
			_rate(value) {
			_rate = 0;
		};
	};

	template<typename T>
	class StateMeasure : public tens::object<T> {
		std::string _name;
		tens::container<T>& _value;
		tens::container<T> _rate;
		tens::container<T> _value_prev;
		tens::container<T> _rate_prev;
		std::weak_ptr<State<T>> _state;
	protected:
		template <typename P>
		void update_rate(P&& rate) {
			std::swap(_rate, _rate_prev);
			_rate = std::forward<P>(rate);
		};

		template <typename P>
		void update_value(P&& value) {
			std::swap(_value, _value_prev);
			_value = std::forward<P>(value);
		};

		virtual void calc_rate(T dt) {
			update_rate((_value - _value_prev) / dt);
		};

		virtual void integrate_value(T dt) {
			update_value(_value_prev + _rate * dt);
		};

		virtual void rate_equation() {};
		virtual void finit_equation() {};

		tens::container<T>& rate() { return _rate; };
		tens::container<T>& value() { return _value; };
		tens::container<T>& rate_prev() { return _rate_prev; };
		tens::container<T>& value_prev() { return _value_prev; };
	public:
		const tens::container<T>& get_rate() const { return _rate; };
		const tens::container<T>& get_value() const { return _value; };
		std::string name() const { return _name; };
		// -- TODO: add move semantic 
		StateMeasure(std::shared_ptr<State<T>>& state, size_t dim, size_t rank, std::string name, tens::FILL_TYPE type = tens::FILL_TYPE::ZERO) :
			tens::object<T>(dim, rank, type, state->basis()),
			_state(state),
			_name(name),
			_value(this->comp()), // link ref
			_rate(dim, rank, tens::FILL_TYPE::ZERO),
			_value_prev(_value),
			_rate_prev(_rate) {
		};

		StateMeasure(StateMeasure&& measure) : 
			tens::object<T>(std::move(measure)),
			_value(this->comp()),
			_name(measure._name),
			_rate(std::move(measure._rate)),
			_value_prev(std::move(measure._value_prev)),
			_rate_prev(std::move(measure._rate_prev)),
			_state(measure._state){
		}

		// access by const ref to other Measures in the State
		const StateMeasure<T>& operator[] (std::string name) {
			if (std::shared_ptr<State<T>> state = this->_state.lock()) {
				return *(*state.get())[name];
			}
			throw new error::StateNotLinked();
		}
	};
};

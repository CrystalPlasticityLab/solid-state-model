#pragma once

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

	template<template<class> class Q, class T>
	class AbstractMeasure {
	protected:
		std::string _name;
		Q<T>& _value;
		Q<T> _rate;
		Q<T> _value_prev;
		Q<T> _rate_prev;

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

		Q<T>& rate() { return _rate; };
		Q<T>& value() { return _value; };
		Q<T>& rate_prev() { return _rate_prev; };
		Q<T>& value_prev() { return _value_prev; };
	public:
		AbstractMeasure(std::string name, Q<T>& value, Q<T>&& rate) :
			_name(name),
			_value(value),
			_rate(std::move(rate)),
			_value_prev(_value),
			_rate_prev(_rate) {};
		std::string name() const { return _name; };
		const Q<T>& get_rate() const { return _rate; };
		const Q<T>& get_value() const { return _value; };

		virtual void rate_equation() = 0;
		virtual void finit_equation() = 0;
		virtual void calc_rate(T dt) {
			update_rate((_value - _value_prev) / dt);
		};

		virtual void integrate_value(T dt) {
			update_value(_value_prev + _rate * dt);
		};
		template<template<class> class Q, class T>
		friend std::ostream& operator<<(std::ostream& out, const AbstractMeasure<Q, T>& m);
	};

	template<template<class> class Q, class T>
	std::ostream& operator<<(std::ostream& out, const AbstractMeasure<Q, T>& m) {
		out << m._name << ": value = " << m._value << ", rate = " << m._rate;
		return out;
	};

	template<typename T>
	class StateMeasure : public tens::object<T>, public AbstractMeasure<tens::container, T> {
		std::weak_ptr<State<T>> _state;
	public:
		// -- TODO: add move semantic 
		StateMeasure(std::shared_ptr<State<T>>& state, size_t dim, size_t rank, std::string name, tens::FILL_TYPE type = tens::FILL_TYPE::ZERO) :
			tens::object<T>(dim, rank, type, state->basis()),
			_state(state),
			AbstractMeasure<tens::container, T>(
				name,
				this->comp(), // link ref
				tens::container<T>(dim, rank, tens::FILL_TYPE::ZERO)) {
		};

		StateMeasure(StateMeasure&& measure) noexcept : 
			tens::object<T>(std::move(measure)),
			AbstractMeasure<tens::container, T>(
				measure._name,
				this->comp(),
				tens::container<T>(std::move(measure._rate))),
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

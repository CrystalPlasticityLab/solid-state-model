#pragma once
#include "../json/json/json.h"

namespace state {
	template<typename T>
	class MaterialPoint;
};

namespace numerical_schema {
}

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

	enum class type_schema {
		RATE_CALCULATE, // dX(n+1) := F(...), X(n+1) = X(n) + dX(n+1)*dt
		FINITE_CALCULATE // X(n+1) := G(...), dX(n+1) = (X(n+1)-X(n))/dt
	};
	extern type_schema DEFAULT_NUMERICAL_SCHEMA;

	template<template<class> class Q, class T>
	class AbstractMeasure {
		std::string _name;
		Q<T>& _value;
		Q<T> _rate;
		Q<T> _value_prev;
		Q<T> _rate_prev;
		int _lock = 0; // to prevent updates during calc step
	protected:
		// use reference rate_temp / value_temp for temporary calculation to prevent a new allocation
		// update_rate()/update_value() use rate_temp / value_temp value as a new one (for more info see declaration)
		// WARNING: rate_temp / value_temp variables updating and calling update_rate()/update_value() must be in one scope
		Q<T> rate_temp;
		Q<T> value_temp;
	public:
		AbstractMeasure(std::string name, Q<T>& value, Q<T>&& rate) :
			_name(name),
			_value(value),
			_rate(std::move(rate)),
			_value_prev(_value),
			_rate_prev(_rate),
			value_temp(_value),
			rate_temp(_rate) {};

		// const refs to private fields
		const Q<T>& rate() const { return _rate; };
		const Q<T>& value() const { return _value; };
		const Q<T>& rate_prev() const { return _rate_prev; };
		const Q<T>& value_prev() const { return _value_prev; };
		const std::string& name() const { return _name; };

		// to prevent modifying rate and value you should lock measure
		int lock() {
			return _lock = rand();
		}
		// to modify rate and value you should unlock measure
		bool unlock(int lock_value) {
			if (lock_value == _lock) {
				_lock = 0;
				return true;
			}
			return false;
		}
		void update_value(const Q<T>& value) {
			value_temp = value;
		};
		// method uses value_temp as a new one (swap pointers, more efficient)
		void update_value() {
#ifdef _DEBUG
			if (_lock) throw new std::logic_error("updates are locked");
#endif
			// _value = _value_temp, _value_temp = _value
			std::swap(_value, value_temp);
			// _value_temp = _value_prev, _value_prev = _value_temp
			std::swap(_value_prev, value_temp);
			// _value = _value_temp, _value_prev = _value, _value_temp = _value_prev
		};

		void update_rate(const Q<T>& rate) {
			rate_temp = rate;
		};
		// method uses rate_temp as a new one (swap pointers, more efficient)
		void update_rate() {
#ifdef _DEBUG
			if (_lock) throw new std::logic_error("updates are locked");
#endif
			std::swap(_rate, rate_temp);
			std::swap(_rate_prev, rate_temp);
		};

		// ============================================================================= //
		//		PERFORMANCE NOTE: to prevent unnecessery allocations use				 //
		//		[+=, -=, /= , *=] methods instead of [+, -, /, *] when possible			 //
		// ============================================================================= //

		// evolution equation in rate form
		// rate must be updated by calling update_rate
		// use reference this->rate_temp() for temporary calculation to prevent a new allocation
		// update_rate() uses this->rate_temp value as a new one (for more info see declaration)
		virtual void rate_equation(T t, T dt) = 0;

		// evolution equation in finite form
		// value must be updated by calling update_value
		// use reference this->value_temp() for temporary calculation to prevent a new allocation
		// update_value() uses this->value_temp value as a new one (for more info see declaration)
		virtual void finite_equation(T t, T dt) = 0;

		virtual T rate_intensity() const = 0;

		virtual T value_intensity() const = 0;

		// the first order schema to calculate rate, may be overriden
		// use reference this->rate_temp() for temporary calculation to prevent a new allocation
		// update_rate() uses this->rate_temp value as a new one (for more info see declaration)
		virtual void calc_rate(T dt) {
			rate_temp = _value;
			rate_temp -= _value_prev;
			rate_temp /= dt;
		};

		// the first order (Euler) schema to integrate value, may be overriden
		// use this->value_temp() for temporary calculation to prevent a new allocation
		// update_value() uses this->value_temp value as a new one (for more info see declaration)
		virtual void integrate_value(T dt) {
			value_temp = _rate;
			value_temp *= dt;
			value_temp += _value;
		};

		template<template<class> class Q, class T>
		friend std::ostream& operator<<(std::ostream& out, const AbstractMeasure<Q, T>& m);
	};

	template<template<class> class Q, class T>
	std::ostream& operator<<(std::ostream& out, const AbstractMeasure<Q, T>& m) {
		out << m._name << ": value = " << m._value << ", rate = " << m._rate;
		return out;
	};


	template <typename T>
	class AbstractSchema_ {
	protected:
		T _t = 0;
	public:
		virtual void init() = 0;
		virtual void calc(T dt) = 0;
		virtual void finalize() = 0;
		void inc_time(T dt) {
			_t += dt;
		}
		T t() { return _t; };
		virtual void step(T dt) {
			init();
			calc(dt);
			finalize();
		};
	};


	template< typename T>
	class StateMeasure : public tens::object<T>, public AbstractMeasure<tens::container, T> {
		const MaterialPoint<T>& _state;
	public:
		StateMeasure(MaterialPoint<T>& state, size_t dim, size_t rank, std::string name, tens::FILL_TYPE type = tens::FILL_TYPE::ZERO) :
			tens::object<T>(dim, rank, type, state.basis()),
			_state(state),
			AbstractMeasure<tens::container, T>(
				name,
				this->comp(), // link ref
				tens::container<T>(dim, rank, tens::FILL_TYPE::ZERO))
		{
		};

		// TODO: looks like a bit weird -> fix
		// becasue StateMeasure is AbstractMeasure
		StateMeasure(StateMeasure&& measure) noexcept : 
			tens::object<T>(std::move(measure)),
			AbstractMeasure<tens::container, T>(
				measure.name(),
				this->comp(),
				tens::container<T>(std::move(measure.rate()))),
			_state(measure._state){
		}

		// access by const ref to other Measures in the State
		const StateMeasure<T>& operator[] (const std::string& name) const {
			_state.get() ? false : new error::StateNotLinked();
			return **(*_state.get())[name];
		}

		const std::shared_ptr<const Json::Value>& param() const {
			return _state.param();
		}
	};
	

	// enum class measure::type_schema {
	// 	RATE_CALCULATE, // dX(n+1) := F(...), X(n+1) = X(n) + dX(n+1)*dt
	// 	FINITE_CALCULATE // X(n+1) := G(...), dX(n+1) = (X(n+1)-X(n))/dt
	// };

	template< typename T, template<class> class Schema = AbstractSchema_>
	class StateMeasureSchema : public StateMeasure<T>, public Schema<T> {
		StateMeasure<T>& measure;
	protected:
		const measure::type_schema _type;
	public:
		StateMeasureSchema(MaterialPoint<T>& state, size_t dim, size_t rank, std::string name, tens::FILL_TYPE type, measure::type_schema type_schema = measure::type_schema::RATE_CALCULATE) :
			StateMeasure<T>(state, dim, rank, name, type),
			Schema<T>(),
			measure(*this),
			_type(type_schema)
		{
		};

		virtual void init() {};

		virtual void calc(T dt) {
#ifdef _DEBUG
			//_key = measure.lock();
#endif
			switch (_type)
			{
			case measure::type_schema::RATE_CALCULATE:
				measure.rate_equation(AbstractSchema_<T>::_t, dt);
				measure.update_rate();
				measure.integrate_value(dt);
				measure.update_value();
				break;
			case measure::type_schema::FINITE_CALCULATE:
				measure.finite_equation(AbstractSchema_<T>::_t, dt);
				measure.update_value();
				measure.calc_rate(dt);
				measure.update_rate();
				break;
			default:
				//throw new error::UndefinedNumericalSchema();
				break;
			}
			AbstractSchema_<T>::inc_time(dt);
		};

		virtual void finalize() override {
#ifdef _DEBUG
			//measure.unlock(_key);
#endif
		};

		void set_numerical_schema(measure::type_schema type) {
			_type = type;
		}
	};
};

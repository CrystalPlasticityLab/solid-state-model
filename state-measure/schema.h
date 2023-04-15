#pragma once
#include "../tensor-matrix/tensor/object.h"
#include "../state-measure/measure.h"

namespace state {
	template<typename T>
	class State;
};

namespace numerical_schema {
	namespace error {
		class ValueUpdated : public std::exception {
		public:
			virtual const char* what() const noexcept {
				return "NumericalSchema value was updated twice";
			};
		};

		class RateUpdated : public std::exception {
		public:
			virtual const char* what() const noexcept {
				return "NumericalSchema rate was updated twice";
			};
		};

		class ValueNotUpdated : public std::exception {
		public:
			virtual const char* what() const noexcept {
				return "NumericalSchema value was not updated";
			};
		};

		class RateNotUpdated : public std::exception {
		public:
			virtual const char* what() const noexcept {
				return "NumericalSchema rate was not updated";
			};
		};

		class NotUpdated : public std::exception {
		public:
			virtual const char* what() const noexcept {
				return "NumericalSchema was not updated";
			};
		};

		class UndefinedNumericalSchema : public std::exception {
		public:
			virtual const char* what() const noexcept {
				return "NumericalSchema is undefined";
			};
		};
	}


	class AbstractSchema {
	protected:
		virtual void init() = 0;
		virtual void calc() = 0;
		virtual void finalize() = 0;
	public:
		virtual void step() {
			init();
			calc();
			finalize();
		};
	};

	enum class type_schema {
		RATE_ASSIGN,    // dX(n+1) := dx_new, X(n+1) = X(n) + dX(n+1)*dt
		RATE_CALCULATE, // dX(n+1) := F(...), X(n+1) = X(n) + dX(n+1)*dt
		FINITE_ASSIGN,  // X(n+1) := x_new, dX(n+1) = (X(n+1)-X(n))/dt
		FINITE_CALCULATE // X(n+1) := G(...), dX(n+1) = (X(n+1)-X(n))/dt
	};

	template <typename T>
	class DefaultSchema : public std::unique_ptr<measure::StateMeasure<T>>, protected AbstractSchema {
		type_schema _type;
		T& _dt;
		int _key = 0;
	public:
		DefaultSchema(type_schema type,
			std::unique_ptr<measure::StateMeasure<T>> &&measure,
			T& dt) :
			_type(type),
			_dt(dt),
			std::unique_ptr<measure::StateMeasure<T>>(std::move(measure)) {
		};

		virtual void init() {};

		virtual void calc() {
			auto& measure = *(*this);
#ifdef _DEBUG
			_key = measure.lock();
#endif
			switch (_type)
			{
			case numerical_schema::type_schema::RATE_CALCULATE:
				measure.rate_equation();
			case numerical_schema::type_schema::RATE_ASSIGN:
				measure.integrate_value(_dt);
				break;
			case numerical_schema::type_schema::FINITE_CALCULATE:
				measure.finit_equation();
			case numerical_schema::type_schema::FINITE_ASSIGN:
				measure.calc_rate(_dt);
				break;
			default:
				throw new error::UndefinedNumericalSchema();
				break;
			}
		};

		virtual void finalize() override {
			auto& measure = *(*this);
#ifdef _DEBUG
			measure.unlock(_key);
#endif
			measure.update_rate();
			measure.update_value();
		};
	};
};

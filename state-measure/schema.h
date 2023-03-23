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
	public:
		virtual void init() = 0;
		virtual void step() = 0;
		virtual void finalize() = 0;
	};

	enum class type_schema {
		RATE_ASSIGN,    // dX(n+1) := dx_new, X(n+1) = X(n) + dX(n+1)*dt
		RATE_INTEGRATE, // dX(n+1) := F(...), X(n+1) = X(n) + dX(n+1)*dt
		FINITE_ASSIGN,  // X(n+1) := x_new, dX(n+1) = (X(n+1)-X(n))/dt
		FINITE_DERIVATE // X(n+1) := G(...), dX(n+1) = (X(n+1)-X(n))/dt
	};

	template <typename T>
	class DefaultSchema : public measure::StateMeasure<T>, public AbstractSchema {
		bool _value_updated = false;
		bool _rate_updated = false;
		type_schema _type;
		T& _dt;
	public:
		DefaultSchema(type_schema type, 
						std::shared_ptr<state::State<T>>& state, 
						size_t dim, size_t rank, 
						std::string name, 
						tens::FILL_TYPE type_fill = tens::FILL_TYPE::ZERO) :
			_type(type),
			_dt(state->dt()),
			measure::StateMeasure<T>(state, dim, rank, name, type_fill) {
			link(state, std::move(*this));
		};

		DefaultSchema(type_schema type,
			measure::StateMeasure<T>&& measure,
			std::shared_ptr<state::State<T>>& state) :
			_type(type),
			_dt(state->dt()),
			measure::StateMeasure<T>(std::move(measure)) {
			link(state, std::move(*this));
		};

		virtual void init() {
			_value_updated = false;
			_rate_updated = false;
		};

		virtual void step() {
			switch (_type)
			{
			case numerical_schema::type_schema::RATE_INTEGRATE:
				rate_equation();
			case numerical_schema::type_schema::RATE_ASSIGN:
				_rate_updated ? this->integrate_value(_dt) : throw new error::RateNotUpdated();
				break;
			case numerical_schema::type_schema::FINITE_DERIVATE:
				finit_equation();
			case numerical_schema::type_schema::FINITE_ASSIGN:
				_value_updated ? this->calc_rate(_dt) : throw new error::RateNotUpdated();
				break;
			default:
				throw new error::UndefinedNumericalSchema();
				break;
			}
		};

		virtual void rate_equation() override {
			// evolution equation in rate form
			// rate must be updated by calling update_rate
		}

		virtual void finit_equation() override {
			// evolution equation in finite form
			// value must be updated by calling update_value
		};

		virtual void finalize() {
		};
	};
};

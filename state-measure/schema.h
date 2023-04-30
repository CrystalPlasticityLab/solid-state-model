#pragma once
#include "../tensor-matrix/tensor/object.h"
#include "../state-measure/measure.h"

namespace state {
	template<typename T>
	class MaterialPoint;
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


	template <typename T>
	class AbstractSchema {
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


	template <typename T>
	class DefaultSchema : public std::unique_ptr<measure::StateMeasure<T>>, protected AbstractSchema<T> {
		measure::type_schema _type;
		int _key = 0;
	protected:
	public:
		DefaultSchema(measure::type_schema type,
			std::unique_ptr<measure::StateMeasure<T>> &&measure) :
			_type(type),
			std::unique_ptr<measure::StateMeasure<T>>(std::move(measure)) {
		};

		virtual void init() {};

		virtual void calc(T dt) {
			auto& measure = *(*this);
#ifdef _DEBUG
			//_key = measure.lock();
#endif
			switch (_type)
			{
			case measure::type_schema::RATE_CALCULATE:
				measure.rate_equation(AbstractSchema<T>::_t, dt);
				measure.update_rate();
				measure.integrate_value(dt);
				measure.update_value();
				break;
			case measure::type_schema::FINITE_CALCULATE:
				measure.finite_equation(AbstractSchema<T>::_t, dt);
				measure.update_value();
				measure.calc_rate(dt);
				measure.update_rate();
				break;
			default:
				throw new error::UndefinedNumericalSchema();
				break;
			}
			AbstractSchema<T>::inc_time(dt);
		};

		virtual void finalize() override {
			auto& measure = *(*this);
#ifdef _DEBUG
			//measure.unlock(_key);
#endif
		};

		void set_numerical_schema(measure::type_schema type) {
			_type = type;
		}
	};


};

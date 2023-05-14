#pragma once
#include "../tensor-matrix/tensor/object.h"
#include "../state-measure/measure.h"

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
			try {
				init();
				calc(dt);
				finalize();
			} catch (const std::exception &e){
				std::cout << "Error during calculation. Reason: " << e.what();
				exit(1);
			}
		};
	};
};

#pragma once
#include "../state-measure/state.h"


namespace measure {
	using namespace state;

	template<typename T>
	class Stress : public StateMeasureSchema<T> {
	public:
		Stress(MaterialPoint<T>& state, measure::type_schema type_schema, const std::string& name) : StateMeasureSchema<T>(state, 3, 2, name, tens::FILL_TYPE::ZERO, type_schema) {};
	};

	namespace stress {
		const std::string CAUCHY = "S";

		template<typename T>
		class CaushyStress : public Stress<T> {
		public:
			CaushyStress(MaterialPoint<T>& state, measure::type_schema type_schema) : Stress<T>(state, type_schema, CAUCHY) {};

			// evolution equation in rate form
			virtual void rate_equation(T t, T dt) override {}

			// evolution equation in finite form
			virtual void finite_equation(T t, T dt) override {};

			virtual T rate_intensity() const override {
				const auto& dS = this->rate();
				return std::sqrt(3 * convolution_transp(dS, dS) / 2);
			}

			virtual T value_intensity() const override {
				const auto& S = this->value();
				return std::sqrt(3 * convolution_transp(S, S) / 2);
			}
			template<class T>
			friend std::ostream& operator<<(std::ostream& out, const CaushyStress<T>& m);
		};
		template<class T>
		std::ostream& operator<<(std::ostream& out, const CaushyStress<T>& m) {
			out << m.name() << ": value = " << m.value() << ", rate = " << m.rate();
			return out;
		};
	}
};

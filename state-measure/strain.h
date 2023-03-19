#pragma once
#include "../state-measure/state.h"

namespace measure {
	using namespace state;

	namespace strain {
		const std::string DEFORM_GRADIENT = "F";
		const std::string CAUCHY_GREEN = "C";

		template<typename T>
		class GradDeform : public Measure<T> {
		public:
			GradDeform(std::shared_ptr<State<T>>& state) : Measure<T>(state, tens::object<T>(IDENT_MATRIX<T>, state->basis()), DEFORM_GRADIENT) {};

			virtual void integrate_value() override {
				auto L = this->rate();
				L *= -this->dt(); // -L*dt
				L += IDENT_MATRIX<T>; // I - L*dt
				Measure<T>::update_value(L.inverse() * this->value_prev());// (I - L * dt)^-1 * F
			}

			std::pair<tens::object<T>, tens::object<T>> polar_decomposition() {
				const auto& F = this->value();
				auto V = F * transpose(F); // TODO: take sqrt !!!
				auto R = F * V.inverse();
		
				return { V, R };
			}
		};
	};
};

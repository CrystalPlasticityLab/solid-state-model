#pragma once
#include "../state-measure/state.h"

namespace measure {
	using namespace state;

	namespace strain {
		const std::string DEFORM_GRADIENT = "F";
		const std::string CAUCHY_GREEN = "C";

		template<typename T>
		class GradDeform : public StateMeasure<T> {
		public:
			GradDeform(std::shared_ptr<State<T>>& state) : StateMeasure<T>(state, 3, 2, DEFORM_GRADIENT, tens::FILL_TYPE::INDENT) {};

			virtual void integrate_value(T dt) override {
				auto L = this->rate();
				L *= -dt; // -L*dt
				L += IDENT_MATRIX<T>; // I - L*dt
				StateMeasure<T>::update_value(L.inverse() * this->value_prev());// (I - L * dt)^-1 * F
			}

			std::pair<tens::container<T>, tens::container<T>> polar_decomposition() {
				const auto& F = this->value();
				auto V = F * F.transpose(); // TODO: take sqrt !!!
				auto R = F * V.inverse();

				return { V, R };
			};
		};
	};

};

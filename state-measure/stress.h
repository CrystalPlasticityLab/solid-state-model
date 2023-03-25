#pragma once
#include "../state-measure/state.h"


namespace measure {
	using namespace state;

	namespace stress {
		const std::string CAUCHY = "S";

		template<typename T>
		class CaushyStress : public StateMeasure<T> {
		public:
			CaushyStress(std::shared_ptr<State<T>>& state) : StateMeasure<T>(state, 3, 2, CAUCHY, tens::FILL_TYPE::ZERO) {};

			virtual void rate_equation() override {
				// evolution equation in rate form
				// rate must be updated by calling update_rate
			}

			virtual void finit_equation() override {
				// evolution equation in finite form
				// value must be updated by calling update_value
			};
		};
	}
};

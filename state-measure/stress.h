#pragma once
#include "../state-measure/state.h"


namespace measure {
	using namespace state;

	namespace stress {
		const std::string CAUCHY = "S";

		template<typename T>
		class CaushyStress : public StateMeasure<T> {
		public:
			CaushyStress(State<T>& state) : StateMeasure<T>(state, 3, 2, CAUCHY, tens::FILL_TYPE::ZERO) {};

			// evolution equation in rate form
			virtual void rate_equation() override {}

			// evolution equation in finite form
			virtual void finit_equation(T t) override {};
		};
	}
};

#pragma once
#include "../state-measure/state.h"


namespace measure {
	using namespace state;

	namespace stress {
		const std::string CAUCHY = "S";

		template<typename T>
		class CaushyStress : public Measure<T> {
		public:
			CaushyStress(std::shared_ptr<State<T>>& state) : Measure<T>(state, tens::object<T>(ZERO_MATRIX<T>, state->basis()), CAUCHY) {};
		};
	}
};

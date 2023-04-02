#pragma once
#include "../state-measure/state.h"

namespace measure {
	using namespace state;

	namespace scalar {
		const std::string SCALAR_ARRAY = "SCAL_ARRAY";

		template<typename T>
		class Scalar : public StateMeasure<T> {
		public:
			Scalar(std::shared_ptr<State<T>>& state, size_t dim) : StateMeasure<T>(state, dim, 0, SCALAR_ARRAY, tens::FILL_TYPE::ZERO) {};

			virtual void rate_equation() override {
				this->rate_temp.fill_value(1.0);
				StateMeasure<T>::update_rate();
			}

			virtual void finit_equation() override {
			};
		};
	};

};

#pragma once
#include "../state-measure/state.h"

namespace measure {
	using namespace state;

	namespace scalar {
		const std::string SCALAR_ARRAY = "SCAL_ARRAY";

		template<typename T, size_t DIM, size_t RANK = 1>
		class Scalar : public StateMeasure<T, DIM, RANK> {
		public:
			Scalar(MaterialPoint<T, DIM>& state, size_t dim) : StateMeasure<T, DIM, RANK>(state, dim, 0, SCALAR_ARRAY, tens::FILL_TYPE::ZERO) {};

			virtual void rate_equation() override {
				this->rate_temp.fill_value(1.0);
			}

			virtual void finite_equation(T t) override {
			};
		};
	};

};

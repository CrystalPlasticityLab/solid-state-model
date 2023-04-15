#pragma once
#include "../state-measure/state.h"

namespace model {
	using namespace state;
	using namespace measure;


	template<typename T>
	class Elasticity : public State<T> {
	public:
		Elasticity(Basis<T>&& basis) : State<T>(std::move(basis)) {
			this->add<measure::scalar::Scalar, T, 10>(numerical_schema::type_schema::RATE_CALCULATE);
			this->add<measure::strain::GradDeform, T>(numerical_schema::type_schema::RATE_CALCULATE);
			this->add<measure::stress::CaushyStress, T>(numerical_schema::type_schema::RATE_CALCULATE);
		};
	};
}
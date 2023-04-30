#pragma once
#include "../state-measure/state.h"
#include "../tensor/object.h"
#include "./relation.h"
#include "./elasticity.h"
#include "./plasticity.h"

namespace model {
	using namespace state;
	using namespace measure;

	class ModelFactory {
	public:
		template <template<class> class ModelType, class T = double>
		static std::shared_ptr<ModelType<T>> create() {//const Json::Value &params, measure::type_schema type) {

			// auto Q = tens::create_basis<T, 3>(tens::DEFAULT_ORTH_BASIS::RANDOM);
			// return std::make_shared <model::Plasticity<strain::GradDeform, stress::CaushyStress, strain::GradDeformInelast, strain::GradDeformElast>>(std::move(Q), params, type);
			// return std::make_shared<model::Elasticity<strain::GradDeform>>(std::move(Q), params, type);
		}
	};
}

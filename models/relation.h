#pragma once
#include "../state-measure/state.h"

namespace model {
	// enum class type {
	// 	RATE,
	// 	FINIT
	// };
	using namespace state;
	using namespace measure;

	template <typename T>
	using pMeasure = std::shared_ptr<DefaultSchema<T>>;

	// binary relatoin G(x) : S(F)
	template<typename T>
	class ElasticRelation : public AbstractSchema<T> {
	protected:
		stress::CaushyStress<T>& S;
		strain::GradDeform<T>& F;
	public:
		ElasticRelation(pMeasure<T>& _S, pMeasure<T>& _F) :
			S(static_cast<stress::CaushyStress<T>&>(**_S.get())),
			F(static_cast<strain::GradDeform<T>&>(**_F.get())) {};
	};
}
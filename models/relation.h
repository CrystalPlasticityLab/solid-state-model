#pragma once
#include "../state-measure/state.h"

namespace model {
	using namespace state;
	using namespace measure;

	template <typename T>
	using MeasureType = DefaultSchema<T>;

	template<typename T>
	class Relation : public std::pair<MeasureType<T>&, MeasureType<T>&>, public AbstractSchema<T> {
	public:
		Relation(std::shared_ptr<DefaultSchema<T>>& _S, std::shared_ptr<DefaultSchema<T>>& _F) :
			std::pair<MeasureType<T>&, MeasureType<T>&>(*_S.get(), *_F.get()){};
	};

	template<typename T>
	class ElasticRelation : private Relation<T> { // dS = Ï:L
	protected:
		measure::stress::CaushyStress<T>& S;
		measure::strain::GradDeform<T>& F;
		DefaultSchema<T>& S_schema;
		DefaultSchema<T>& F_schema;
	public:
		ElasticRelation(std::shared_ptr<DefaultSchema<T>>& _S, std::shared_ptr<DefaultSchema<T>>& _F) : 
			Relation<T>(_S, _F),
			S(static_cast<measure::stress::CaushyStress<T>&>(**_S.get())),
			F(static_cast<measure::strain::GradDeform<T>&>(**_F.get())),
			S_schema(*_S.get()),
			F_schema(*_F.get()) {};
	};

}
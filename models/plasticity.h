#pragma once
#include "./elasticity.h"

namespace model {
	using namespace state;
	using namespace measure;

	/*
		Base class for elasto-plastic behavior material, inherited from Elasticity
		Plasticity<ElasticType, StrainMeasure>, where:
			- StrainMeasure - specific strain measure,
			- StressMeasure - specific stress measure, 
	*/
	
	template<
		template<class> class StrainMeasure,
		template<class> class StressMeasure,
		class T = double>
	class Plasticity : public Elasticity<StrainMeasure, StressMeasure, T> {
	protected:
		std::shared_ptr<StrainMeasure<T>> F_in;
		std::shared_ptr<StrainMeasure<T>> F_e;
	public:
		Plasticity(const json& params, measure::type_schema type) :
			Elasticity<StrainMeasure, StressMeasure, T>(params, type)
		{
			const auto curve = parse_json_value<std::vector<std::pair<T, T>>>("curve", params);
			const auto treshold = parse_json_value<T>("flow_treshold", params);
			F_in = std::make_shared<PlasticRelation<StressMeasure, StrainMeasure, T>>(type, *this, this->S, this->F, curve, this->elast_modulus[1], treshold);
			F_e = std::make_shared<StrainDecomposition<StrainMeasure, T>>(type, *this, this->F, this->F_in);
			this->reset_elastic_strain_measure(this->F_e); // change S(F) -> S(F_e)
		};

		virtual void calc(T dt) override {
			this->F->calc(dt); // loading -> L / F
			this->F_in->calc(dt);  // plastic relation -> L_in / F_in
			this->F_e->calc(dt); // strain decomposition L_e / F_e
			this->S->calc(dt); // elastic relation -> S_rate / S
		};

		virtual std::ostream& print_measures(std::ostream& out) const override {
			out << *this->F << std::endl;
			out << *this->F_in << std::endl;
			out << *this->F_e << std::endl;
			out << *this->S << std::endl;
			return out;
		}
	};
};
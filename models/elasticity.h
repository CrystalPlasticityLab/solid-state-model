#pragma once
#include "../state-measure/state.h"
#include "./relation.h"

namespace model {
	using namespace state;
	using namespace measure;

	/*
		Base class for elastic behavior material
		Elasticity<ElasticType, StrainMeasure, StressMeasure>, where:
			- StrainMeasure - specific strain measure,
			- StressMeasure - specific stress measure
	*/
	template<
		template<class> class StrainMeasure = strain::GradDeform, 
		template<class> class StressMeasure = stress::CaushyStress, class T = double>
	class Elasticity : public MaterialPoint<T, 3> {
	protected:
		std::shared_ptr<StrainMeasure<T>> F;
		std::shared_ptr<StressMeasure<T>> S;
		void parse_json_params(const Json::Value& params) {
			MaterialPoint<T, 3>::parse_json_params(params);
		}

		void reset_elastic_strain_measure(std::shared_ptr<const StrainMeasure<T>> new_F) {
			static_cast<ElasticRelation<StressMeasure, StrainMeasure, T>&>(*S).reset_elastic_strain_measure(new_F);
		}
	public:
		Elasticity(const Json::Value& params, measure::type_schema type) : 
			MaterialPoint<T, 3>(params, type),
			F(std::make_shared<StrainMeasure<T>>(*this, type)),
			S(std::make_shared<ElasticRelation<StressMeasure, StrainMeasure, T>>(type, *this, this->F))
		{
			parse_json_params(params);
		};

		virtual void calc(T dt) override {
			F->calc(dt);
			S->calc(dt);
		};

		virtual std::ostream& print_measures(std::ostream& out) const override {
			out << *F << std::endl; 
			out << *S << std::endl;
			return out;
		}
	};
}

namespace measure {
	namespace strain {
		// assignment a new rate L
		template<typename T, size_t DIM>
		void measure::strain::GradDeform<T, DIM>::rate_equation(T t, T dt) {
			auto& L = this->rate_temp;
			L.fill_value(tens::FILL_TYPE::INDENT);
		}

		// assignment a new value F 
		template<typename T, size_t DIM>
		void measure::strain::GradDeform<T, DIM>::finite_equation(T t, T dt) {
			auto& F = this->value_temp;
			F.fill_value(tens::FILL_TYPE::INDENT);
			F *= T(1) + t;
		}; 
	};
};
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
		template<typename T> class StrainMeasure,
		template<typename T> class StressMeasure, typename T>
	class Elasticity : public MaterialPoint<T, 3> {
	protected:
		std::array<T, 2> elast_modulus;
		std::shared_ptr<StrainMeasure<T>> F;
		std::shared_ptr<StressMeasure<T>> S;
		void reset_elastic_strain_measure(std::shared_ptr<const StrainMeasure<T>> new_F) {
			static_cast<ElasticRelation<StressMeasure, StrainMeasure, T>&>(*S).reset_elastic_strain_measure(new_F);
		}
	public:
		Elasticity(const json& params, measure::type_schema type) : 
			MaterialPoint<T, 3>(params, type),
			F(std::make_shared<StrainMeasure<T>>(*this, type))
		{
			elast_modulus = parse_json_value<std::array<T, 2>>("elast_modulus", params);
			S = std::make_shared<ElasticRelation<StressMeasure, StrainMeasure, T>>(type, *this, this->F, elast_modulus);
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
		template<typename T>
		void measure::strain::GradDeform<T>::rate_equation(T t, T dt) {
			auto& L = this->rate_temp;
			L.fill_value(tens::FILL_TYPE::ZERO);
			L[0] = 2;
			L[1] = -1;
			L[2] = -1;
		}

		// assignment a new value F 
		template<typename T>
		void measure::strain::GradDeform<T>::finite_equation(T t, T dt) {
			auto& F = this->value_temp;
			F.fill_value(tens::FILL_TYPE::ZERO);
			F[0] = T(1) + t;
			F[1] = F[2] = std::sqrt(T(1) / F[0]);
		}; 
	};
};
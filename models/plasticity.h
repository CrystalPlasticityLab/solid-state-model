#pragma once
#include "./elasticity.h"

namespace model {
	using namespace state;
	using namespace measure;

	template <typename T>
	using pMeasure = std::shared_ptr<DefaultSchema<T>>;

	/*
		Base class for elasto-plastic behavior material, inherited from Elasticity
		Elasticity<ElasticType, StrainMeasure>, where:
			- StrainMeasure - specific strain measure,
			- StressMeasure - specific stress measure, 
			- ElasticStrainMeasure - specific measure of elastic strain
			- PlasticStrainMeasure - specific measure of inelastic strain
	*/
	template<
		template<class> class StrainMeasure,
		template<class> class StressMeasure,
		template<class> class PlasticStrainMeasure, 
		template<class> class ElasticStrainMeasure,
		class T = double>
	class Plasticity : public Elasticity<ElasticStrainMeasure, StrainMeasure, StressMeasure, T> {
	protected:
		pMeasure<T> F_in;
	public:
		Plasticity(Basis<T>&& basis, numerical_schema::type_schema type) : Elasticity<ElasticStrainMeasure, StrainMeasure, StressMeasure, T>(std::move(basis), type) {
			F_in = this->add<PlasticStrainMeasure, T>(type);
		};

		virtual void calc(T dt) override {
			this->F->calc(dt); // L / F
			F_in->calc(dt); // L_in / F_in
			this->F_e->calc(dt); // L_e / F_e
			this->elastic_relation->calc(dt); // S rate
			this->S->calc(dt); // S
		};
	};
};

namespace measure {
	namespace strain {
		// assignment a new rate L_in
		template<typename T>
		void measure::strain::GradDeformInelast<T>::rate_equation() {
			const auto& L = (*this)[measure::strain::DEFORM_GRADIENT].rate();
			auto iS = (*this)[measure::stress::CAUCHY].value_intensity();
			auto& L_in = this->rate_temp = L;
			T threshold = 0.014;
			iS > threshold ? L_in : L_in *= (iS / threshold);
		}

		// assignment a new value F_in
		template<typename T>
		void measure::strain::GradDeformInelast<T>::finit_equation(T t) {
			const auto& F = (*this)[measure::strain::DEFORM_GRADIENT].value();
			auto iS = (*this)[measure::stress::CAUCHY].value_intensity();
			auto& F_in = this->value_temp = F;
			T threshold = 0.014;
			iS > threshold ? F_in : F_in = IDENT_MATRIX<T> + (F - IDENT_MATRIX<T>) * (iS / threshold);
		};

		// assignment a new rate L_e
		template<typename T>
		void measure::strain::GradDeformElast<T>::rate_equation() {
			const auto& F_e = (*this)[measure::strain::DEFORM_GRADIENT_ELAST].value();

			auto& L_e = this->rate_temp = (*this)[measure::strain::DEFORM_GRADIENT].rate(); // L
			auto L_in = F_e; // F_e
			L_in *= (*this)[measure::strain::DEFORM_GRADIENT_INELAST].rate(); // F_e.L_in
			L_in *= F_e.inverse(); // F_e.L_in.F_e^-1
			L_e -= L_in; //  L - F_e.L_in.F_e^-1
		}

		// assignment a new value F_e
		template<typename T>
		void measure::strain::GradDeformElast<T>::finit_equation(T t) {
			const auto& F = (*this)[measure::strain::DEFORM_GRADIENT].value();
			const auto& F_in = (*this)[measure::strain::DEFORM_GRADIENT_INELAST].value().inverse();
			this->value_temp = F * F_in;
		};
	};
};
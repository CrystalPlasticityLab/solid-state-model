#pragma once
#include "../state-measure/state.h"
#include "./relation.h"

namespace model {
	using namespace state;
	using namespace measure;
	template <typename T>
	using pMeasure = std::shared_ptr<DefaultSchema<T>>;

	/*
		Base class of model
		Model<StrainMeasure, StressMeasure>, where:
			- StrainMeasure - specific strain measure,
			- StressMeasure - specific stress measure
	*/
	template<template<class> class StrainMeasure = strain::GradDeform, template<class> class StressMeasure = stress::CaushyStress, class T = double>
	class Model : public State<T> {
	protected:
		pMeasure<T> F;
		pMeasure<T> S;
	public:
		Model(Basis<T>&& basis, numerical_schema::type_schema type) : State<T>(std::move(basis)) {
			F = this->add<StrainMeasure, T>(type);
			S = this->add<StressMeasure, T>(type);
		};
	};

	/*
		Base class for elastic behavior material
		Elasticity<ElasticType, StrainMeasure, StressMeasure>, where:
			- StrainMeasure - specific strain measure,
			- StressMeasure - specific stress measure
	*/
	template<
		template<class> class ElasticStrainMeasure,
		template<class> class StrainMeasure = strain::GradDeform, 
		template<class> class StressMeasure = stress::CaushyStress, class T = double>
	class Elasticity : public Model<StrainMeasure, StressMeasure, T> {
	protected:
		std::unique_ptr<ElasticRelation<T>> elastic_relation;
		pMeasure<T> F_e;
	public:
		Elasticity(Basis<T>&& basis, numerical_schema::type_schema type) : Model<StrainMeasure, StressMeasure, T>(std::move(basis), type) {
			F_e = this->add<ElasticStrainMeasure, T>(type);
			((type == numerical_schema::type_schema::RATE_CALCULATE) ? 
				elastic_relation = std::make_unique<HypoElastic<T>>(this->S, this->F) : 
				elastic_relation = std::make_unique<HyperElastic<T>>(this->S, this->F));
		};

		virtual void calc(T dt) override {
			this->F->calc(dt);
			elastic_relation->calc(dt);
			this->S->calc(dt);
		};
	};

	template<typename T>
	class HypoElastic : public ElasticRelation<T> { // dS = Ï:L
	public:
		HypoElastic(pMeasure<T>& _S, pMeasure<T>& _F) :
			ElasticRelation<T>(_S, _F) {
			_S->set_numerical_schema(numerical_schema::type_schema::RATE_CALCULATE);
		};

		virtual void init() override {};

		virtual void calc(T dt) override {
			const auto& L = ElasticRelation<T>::F.rate();
			ElasticRelation<T>::S.update_rate(L * 0.13);
		};
		virtual void finalize() override {};
	};

	template<typename T>
	class HyperElastic : public ElasticRelation<T> { // S = Ï:E
	public:
		HyperElastic(pMeasure<T>& _S, pMeasure<T>& _F) :
			ElasticRelation<T>(_S, _F) {
			_S->set_numerical_schema(numerical_schema::type_schema::FINITE_CALCULATE);
		};

		virtual void init() override {};

		virtual void calc(T dt) override {
			const auto& C = ElasticRelation<T>::F.lagrangian_strain_tensor();
			ElasticRelation<T>::S.update_value(C * 0.13);
		};

		virtual void finalize() override {};
	};
}

namespace measure {
	namespace strain {
		// assignment a new rate L
		template<typename T>
		void measure::strain::GradDeform<T>::rate_equation() {
			auto& L = this->rate_temp;
			L.fill_value(tens::FILL_TYPE::INDENT);
		}

		// assignment a new value F 
		template<typename T>
		void measure::strain::GradDeform<T>::finit_equation(T t) {
			auto& F = this->value_temp;
			F.fill_value(tens::FILL_TYPE::INDENT);
			F *= T(1) + t;
		};
	};
};
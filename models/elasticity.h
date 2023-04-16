#pragma once
#include "../state-measure/state.h"
#include "./relation.h"

namespace model {
	using namespace state;
	using namespace measure;

	template<template<class> class R, class T>
	class Elasticity : public State<T> {
		std::unique_ptr<R<T>> relation;
	public:
		Elasticity(Basis<T>&& basis) : State<T>(std::move(basis)) {
			this->add<measure::strain::GradDeform, T>();
			this->add<measure::stress::CaushyStress, T>();
			relation = std::make_unique<R<T>>(State<T>::state[measure::stress::CAUCHY], State<T>::state[measure::strain::DEFORM_GRADIENT]);
		};

		virtual void calc(T dt) override {
			relation->calc(dt);
		};
	};

	template<typename T>
	class HypoElasticRelation : public ElasticRelation<T> { // dS = Ï:L
	public:
		HypoElasticRelation(std::shared_ptr<DefaultSchema<T>>& _S, std::shared_ptr<DefaultSchema<T>>& _F) :
			ElasticRelation<T>(_S, _F) {
			_S->set_numerical_schema(numerical_schema::type_schema::RATE_CALCULATE);
		};

		virtual void init() override {};

		virtual void calc(T dt) override {
			ElasticRelation<T>::F_schema.calc(dt);
			const auto& L = ElasticRelation<T>::F.rate();
			ElasticRelation<T>::S.update_rate(L * 0.13);
			ElasticRelation<T>::S_schema.calc(dt);
		};
		virtual void finalize() override {};
	};

	template<typename T>
	class HyperElasticRelation : public ElasticRelation<T> { // S = Ï:E
	public:
		HyperElasticRelation(std::shared_ptr<DefaultSchema<T>>& _S, std::shared_ptr<DefaultSchema<T>>& _F) :
			ElasticRelation<T>(_S, _F) {
			_S->set_numerical_schema(numerical_schema::type_schema::FINITE_CALCULATE);
		};

		virtual void init() override {};

		virtual void calc(T dt) override {
			ElasticRelation<T>::F_schema.calc(dt);
			const auto& C = ElasticRelation<T>::F.lagrangian_strain_tensor();
			ElasticRelation<T>::S.update_value(C * 0.13);
			ElasticRelation<T>::S_schema.calc(dt);
		};

		virtual void finalize() override {};
	};
}

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
	F *= t;
};
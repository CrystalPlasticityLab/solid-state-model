#pragma once
#include "../state-measure/state.h"

namespace model {
	using namespace state;
	using namespace measure;

	// binary relatoin G(x) : S(F)
	template<template<class> class StressMeasure, template<class> class StrainMeasure, typename T>
	class ElasticRelation : public StressMeasure<T> {
	protected:
		std::shared_ptr<const StrainMeasure<T>> F_e;
	public:
		ElasticRelation(measure::type_schema type, MaterialPoint<T, 3>& state, 
			const std::shared_ptr<const StrainMeasure<T>>& _F_e) :
			StressMeasure<T>(state, type),
			F_e(_F_e) 
		{
		};

		void reset_elastic_strain_measure(std::shared_ptr<const StrainMeasure<T>>& new_F_e) {
			F_e.reset();
			F_e = new_F_e;
		}

		virtual void finite_equation(T t, T dt) {
			const auto& C = F_e->lagrangian_strain_tensor();
			this->update_value(C * 0.13);
		};
		virtual void rate_equation(T t, T dt) {
			const auto& L = F_e->rate();
			this->update_rate(L * 0.13);
		}
	};

	template<template<class> class StressMeasure, template<class> class StrainMeasure, typename T>
	class PlasticRelation : public StrainMeasure<T> {
	protected:
		std::shared_ptr<const StressMeasure<T>> S;
		std::shared_ptr<const StrainMeasure<T>> F;
	public:
		PlasticRelation(measure::type_schema type, MaterialPoint<T, 3>& state, 
			const std::shared_ptr<const StressMeasure<T>> _S, 
			const std::shared_ptr<const StrainMeasure<T>> _F) :
			StrainMeasure<T>(state, type, "F_in"),
			S(_S),
			F(_F)
		{};

		virtual void rate_equation(T t, T dt) override {
			const auto& L = F->rate();
			auto iS = S->value_intensity();
			auto& L_in = this->rate_temp = L;
			T threshold = 0.014;
			iS > threshold ? L_in : L_in *= (iS / threshold);
		};

		virtual void finite_equation(T t, T dt) override {
			const auto& F = this->F->value();
			auto iS = S->value_intensity();
			auto& F_in = this->value_temp;// = F;
			T threshold = 0.014;
			if (iS > threshold) {
				F_in = F;
			} else {
				F_in = (this->value() * 90.0 + IDENT_MATRIX<T, 3> + F * (iS / threshold))/91.0;
				//F_in = tens::func(IDENT_MATRIX<T> + F * (iS / threshold), sqrt);
			}
			//iS > threshold ? F_in : F_in = IDENT_MATRIX<T> + (F - IDENT_MATRIX<T>) * (iS / threshold);
		};
	};

	template< template<class> class StrainMeasure, typename T>
	class StrainDecomposition : public StrainMeasure<T> {
	protected:
		std::shared_ptr<const StrainMeasure<T>> F;
		std::shared_ptr<const StrainMeasure<T>> F_in;
	public:
		StrainDecomposition(measure::type_schema type, MaterialPoint<T, 3>& state, 
			const std::shared_ptr<const StrainMeasure<T>> _F, 
			const std::shared_ptr<const StrainMeasure<T>> _F_in) :
			StrainMeasure<T>(state, type, "F_e"),
			F(_F),
			F_in(_F_in)
		{};

		virtual void rate_equation(T t, T dt) override {
			const auto& F_e = this->value();
			auto& L_e = this->rate_temp = F->rate(); // L
			auto L_in = F_e; // F_e
			L_in *= this->F_in->rate(); // F_e.L_in
			L_in *= F_e.inverse(); // F_e.L_in.F_e^-1
			L_e -= L_in; //  L - F_e.L_in.F_e^-1
		};

		virtual void finite_equation(T t, T dt) override {
			this->value_temp = F->value() * F_in->value().inverse();
		};
	};
}
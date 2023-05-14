#pragma once
#include "../state-measure/state.h"
#include "./helpers/model_utils.h"

namespace model {
	using namespace state;
	using namespace measure;

	// binary relatoin G(x) : S(F)
	template<template<class> class StressMeasure, template<class> class StrainMeasure, typename T>
	class ElasticRelation : public StressMeasure<T> {
	private:
		void set_elast_modules_l_mu(T l, T mu) {
			c[0] = c[1] = c[2] = l + 2 * mu;
			c[3] = c[4] = c[5] = l;
			c[6] = c[7] = c[8] = 2 * mu;
		}
	protected:
		std::shared_ptr<const StrainMeasure<T>> F_e;
		std::array<T, 9> c;
		tens::M3x3<T> strain_mul_modules(const tens::M3x3<T> &e) {
			/*  -------------- s = c.e ------------------
			    c0 c5 c4  0  0  0   |e0|   |s0|
				c5 c1 c3  0  0  0   |e1|   |s1|
				c4 c3 c2  0  0  0   |e2|   |s2|
				 0  0  0  c6 0  0 * |e3| = |s3|
				 0  0  0  0  c7 0   |e4|   |s4|
				 0  0  0  0  0  c8  |e5|   |s5|
			*/
			return tens::M3x3<T>({ c[0] * e[0] + c[5] * e[1] + c[4] * e[2],
								   c[5] * e[0] + c[1] * e[1] + c[3] * e[2],
								   c[4] * e[0] + c[3] * e[1] + c[2] * e[2],
								   c[6] * (e[3] + e[6]),  c[7] * (e[4] + e[7]),  c[8] * (e[5] + e[8]),
						           c[6] * (e[3] + e[6]),  c[7] * (e[4] + e[7]),  c[8] * (e[5] + e[8])});
		}
	public:
		ElasticRelation(measure::type_schema type, MaterialPoint<T, 3>& state, 
			const std::shared_ptr<const StrainMeasure<T>>& _F_e, const std::array<T, 2>& _c) :
			StressMeasure<T>(state, type), 
			F_e(_F_e),
			c{T(0)}
		{
			set_elast_modules_l_mu(_c[0], _c[1]);
		};

		void reset_elastic_strain_measure(std::shared_ptr<const StrainMeasure<T>>& new_F_e) {
			F_e.reset();
			F_e = new_F_e;
		}

		virtual void finite_equation(T t, T dt) {
			const auto& H = F_e->right_hencky();
			this->update_value(strain_mul_modules(H));
		};
		virtual void rate_equation(T t, T dt) {
			const auto& L = F_e->rate();
			this->update_rate(strain_mul_modules(L));
		}
	};

	template<template<class> class StressMeasure, template<class> class StrainMeasure, typename T>
	class PlasticRelation : public StrainMeasure<T> {
	protected:
		std::shared_ptr<const StressMeasure<T>> S;
		std::shared_ptr<const StrainMeasure<T>> F;
		const Curve<T> curve;
		const T mu;
		const T flow_treshold;
	public:
		PlasticRelation(measure::type_schema type, MaterialPoint<T, 3>& state, 
			const std::shared_ptr<const StressMeasure<T>> _S, 
			const std::shared_ptr<const StrainMeasure<T>> _F,
			const std::vector<std::pair<T, T>>& _curve, T _mu, T _flow_treshold) :
			StrainMeasure<T>(state, type, "F_in"),
			S(_S), F(_F), curve(_curve), mu(_mu), flow_treshold(_flow_treshold)
		{
		};

		virtual void rate_equation(T t, T dt) override {
			const auto& L = this->F->rate();
			const auto iS = S->value_intensity();
			auto& L_in = this->rate_temp;
			if (iS < flow_treshold) {
				L_in.fill(0);
			} else {
				L_in = L;
				const auto iE = this->F->value_intensity();
				const T plastic_part = curve.value(iE);
				L_in *= plastic_part;
			}
		};

		virtual void finite_equation(T t, T dt) override {
			const auto& F = this->F->value();
			const auto iS = S->value_intensity();
			auto& F_in = this->value_temp = I3x3<T>;
			if (iS > flow_treshold) {
				const auto iE = this->F->value_intensity();
				const T plastic_part = curve.value(iE);
				F_in += (F - I3x3<T>) * plastic_part;
			}
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
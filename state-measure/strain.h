#pragma once
#include "../state-measure/state.h"

namespace measure {
	using namespace state;

	namespace strain {
		const std::string DEFORM_GRADIENT = "F";

		// Deformation gradient tensor
		// see https://en.wikipedia.org/wiki/Finite_strain_theory for more information
		template<typename T>
		class GradDeform : public StateMeasure<T> {
			tens::object<T> E;  //  (Ft*F-I)/2
			tens::object<T> dE; //  dE/dt = Ft*(L+Lt)*F/2
		public:
			GradDeform(std::shared_ptr<State<T>>& state) : 
				StateMeasure<T>(state, 3, 2, DEFORM_GRADIENT, tens::FILL_TYPE::INDENT), 
				 E(3, 2, tens::FILL_TYPE::ZERO, state->basis()),
				dE(3, 2, tens::FILL_TYPE::ZERO, state->basis()) {};

			// calc a new value F
			virtual void integrate_value(T dt) override {
				auto& df = this->value_temp = this->rate();
				(df *= -dt) += IDENT_MATRIX<T>; // I - L*dt
				inverse(df); // (I - L * dt)^-1 * F
				df *= this->value();
			};

			// calc a new rate L
			virtual void calc_rate(T dt) override {
				auto& L = this->rate_temp = this->value_prev(); // fn_1
				L *= this->value().inverse(); // fn_1 * fn^-1
				(L -= IDENT_MATRIX<T>) /= (-dt); // (I - fn_1 * fn^-1)/ dt
			};

			// assignment a new rate L
			virtual void rate_equation() override {
				auto& L = this->rate_temp;
				L.fill_value(tens::FILL_TYPE::INDENT);
			};

			// assignment a new value F 
			virtual void finit_equation() override {
				auto& F = this->value_temp;
				F.fill_value(tens::FILL_TYPE::INDENT);
			};

			// ---------------------------------- helper functions ----------------------------------------
			std::pair<tens::container<T>, tens::container<T>> polar_decomposition() {
				const auto& F = this->value();
				auto V = left_stretch_tensor(); 
				auto R = F * V.inverse();
				return { V, R };
			};

			// The right Cauchy–Green deformation tensor, Ft.F
			tens::container<T> right_cauchy_green() {
				const auto& F = this->value();
				return F * F.transpose();
			}
			// V the right stretch tensor
			tens::container<T> left_stretch_tensor() {
				return func(right_cauchy_green(), sqrt); // sqrt(F.Ft)
			}

			// The left Cauchy–Green deformation tensor, F.Ft
			tens::container<T> left_cauchy_green() {
				const auto& F = this->value();
				return F.transpose() * F;
			}
			// U the left stretch tensor
			tens::container<T> right_stretch_tensor() {
				return func(left_cauchy_green(), sqrt); // sqrt(Ft.F)
			}

			//  (Ft*F-I)/2
			const tens::object<T>& lagrangian_strain_tensor() {
				E = this->value_temp;
				const auto& F = this->value();
				E = F.transpose(); // F
				E *= F; // Ft*F
				E -= IDENT_MATRIX<T>;
				return E *= T(0.5);
			}
		
			//  dE/dt = Ft*(L+Lt)*F/2
			const tens::object<T>& lagrangian_strain_rate_tensor() {
				const auto& F = this->value();
				dE = F.transpose();
				dE *= (this->rate().symmetrize() *= F);
				return dE *= T(0.5);
			}
		};
	};
};

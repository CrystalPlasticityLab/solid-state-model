#pragma once
#include "../state-measure/state.h"

namespace measure {
	using namespace state;

	namespace strain {
		const std::string DEFORM_GRADIENT = "F";

		// Deformation gradient tensor
		// see https://en.wikipedia.org/wiki/Finite_strain_theory for more information
		template<typename T = double>
		class GradDeform : public StateMeasure<T> {
			tens::container<T> E;  //  (Ft*F-I)/2
			tens::container<T> dE; //  dE/dt = Ft*(L+Lt)*F/2
		public:
			GradDeform(State<T>& state, const std::string& name = DEFORM_GRADIENT) :
				StateMeasure<T>(state, 3, 2, name, tens::FILL_TYPE::INDENT), 
				 E(3, 2, tens::FILL_TYPE::ZERO),
				dE(3, 2, tens::FILL_TYPE::ZERO) {};

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
			virtual void rate_equation() override;

			// assignment a new value F 
			virtual void finit_equation(T t) override;

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
			const tens::container<T>& lagrangian_strain_tensor() {
				E = this->value_temp;
				const auto& F = this->value();
				E = F.transpose(); // F
				E *= F; // Ft*F
				E -= IDENT_MATRIX<T>;
				return E *= T(0.5);
			}
		
			//  dE/dt = Ft*(L+Lt)*F/2
			const tens::container<T>& lagrangian_strain_rate_tensor() {
				const auto& F = this->value();
				dE = F.transpose();
				dE *= (this->rate().symmetrize() *= F);
				return dE *= T(0.5);
			}
		};

		const std::string DEFORM_GRADIENT_ELAST = "F_e";
		template<typename T = double>
		class GradDeformElast : public GradDeform<T> {
		public:
			GradDeformElast(State<T>& state) : GradDeform<T>(state, DEFORM_GRADIENT_ELAST) {};
			// assignment a new rate L
			virtual void rate_equation() override;

			// assignment a new value F 
			virtual void finit_equation(T t) override;
		};

		const std::string DEFORM_GRADIENT_INELAST = "F_in";
		template<typename T = double>
		class GradDeformInelast : public GradDeform<T> {
		public:
			GradDeformInelast(State<T>& state) : GradDeform<T>(state, DEFORM_GRADIENT_INELAST) {};
			// assignment a new rate L
			virtual void rate_equation() override;

			// assignment a new value F 
			virtual void finit_equation(T t) override;
		};
	};
};

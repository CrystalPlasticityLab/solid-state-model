#pragma once
#include "../state-measure/state.h"
#include "../state-measure/stress.h"

namespace measure {
	using namespace state;

	namespace strain {
		const std::string DEFORM_GRADIENT = "F";

		// Deformation gradient tensor
		// see https://en.wikipedia.org/wiki/Finite_strain_theory for more information
		template<typename T = double>
		class GradDeform : public StateMeasureSchema<T> {
			const std::shared_ptr<tens::container<T>> E;  //  (Ft*F-I)/2
			const std::shared_ptr<tens::container<T>> dE; //  dE/dt = Ft*(L+Lt)*F/2
		public:
			GradDeform(MaterialPoint<T>& state, measure::type_schema type_schema, const std::string& name = DEFORM_GRADIENT) :
				StateMeasureSchema<T>(state, 3, 2, name, tens::FILL_TYPE::INDENT, type_schema),
				 E(std::make_shared<tens::container<T>>(3, 2, tens::FILL_TYPE::ZERO)),
				dE(std::make_shared<tens::container<T>>(3, 2, tens::FILL_TYPE::ZERO)) {};

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
			virtual void rate_equation(T t, T dt) override;

			// assignment a new value F 
			virtual void finite_equation(T t, T dt) override;

			// ---------------------------------- helper const methods ----------------------------------------
			virtual T rate_intensity() const override {
				const auto& L = this->rate();
				return std::sqrt(2*convolution_transp(L, L)/3);
			}

			virtual T value_intensity() const override {
				const auto E = lagrangian_strain_tensor();
				return std::sqrt(2 * convolution_transp(E, E) / 3);
			}

			std::pair<tens::container<T>, tens::container<T>> polar_decomposition() const {
				const auto& F = this->value();
				auto V = left_stretch_tensor(); 
				auto R = F * V.inverse();
				return { V, R };
			};

			// The right Cauchy–Green deformation tensor, Ft.F
			tens::container<T> right_cauchy_green() const {
				const auto& F = this->value();
				return F * F.transpose();
			}
			// V the right stretch tensor
			tens::container<T> left_stretch_tensor() const {
				return func(right_cauchy_green(), std::sqrt); // sqrt(F.Ft)
			}

			// The left Cauchy–Green deformation tensor, F.Ft
			tens::container<T> left_cauchy_green() {
				const auto& F = this->value();
				return F.transpose() * F;
			}
			// U the left stretch tensor
			tens::container<T> right_stretch_tensor() const {
				return func(left_cauchy_green(), std::sqrt); // sqrt(Ft.F)
			}

			//  (Ft*F-I)/2
			const tens::container<T>& lagrangian_strain_tensor() const {
				const auto& F = this->value();
				*E = F.transpose() * F; // F.Ft
				*E -= IDENT_MATRIX<T>;
				return *E *= T(0.5);
			}
		
			//  dE/dt = Ft*(L+Lt)*F/2
			const tens::container<T>& lagrangian_strain_rate_tensor() const {
				const auto& F = this->value();
				*dE = F.transpose();
				*dE *= (this->rate().symmetrize() *= F);
				return dE *= T(0.5);
			}
			template<class T>
			friend std::ostream& operator<<(std::ostream& out, const GradDeform<T>& m);
		};

		template<class T>
		std::ostream& operator<<(std::ostream& out, const GradDeform<T>& m) {
			out << m.name() << ": value = " << m.value() << ", rate = " << m.rate();
			return out; 
		};
		/*
		const std::string DEFORM_GRADIENT_INELAST = "F_in";
		const std::string DEFORM_GRADIENT_ELAST = "F_e";
		template<typename T = double>
		class GradDeformElast : public GradDeform<T> {
			//const GradDeform<T>& F;
			//const GradDeform<T>& F_in;
		public:
			GradDeformElast(State<T>& state) : 
				GradDeform<T>(state, DEFORM_GRADIENT_ELAST)//,
				//F(static_cast<const GradDeform<T>&>((*this)[measure::strain::DEFORM_GRADIENT])),
				//F_in(static_cast<const GradDeform<T>&>((*this)[measure::strain::DEFORM_GRADIENT_INELAST]))
			  {};
			// assignment a new rate L
			virtual void rate_equation(T t, T dt) override;

			// assignment a new value F 
			virtual void finite_equation(T t, T dt) override;
		};

		template<typename T = double>
		class GradDeformInelast : public GradDeform<T> {
			//const GradDeform<T>& F;
			//const measure::Stress<T>& S;
		public:
			GradDeformInelast(State<T>& state) :
				GradDeform<T>(state, DEFORM_GRADIENT_INELAST)//,
				//F(static_cast<const GradDeform<T>&>((*this)[measure::strain::DEFORM_GRADIENT])),
				//S(static_cast<const measure::Stress<T>&>((*this)[measure::stress::CAUCHY]))
			{};
			// assignment a new rate L 
			virtual void rate_equation(T t, T dt) override;

			// assignment a new value F 
			virtual void finite_equation(T t, T dt) override;
		};
		*/
	};
};

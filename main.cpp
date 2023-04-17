#pragma once
#include "tensor/test/expect.h"
#include "tensor/test/test.h"
#include "tensor/object.h"
#include "tensor/quat.h"
#include "./state-measure/schema.h"
#include "./state-measure/state.h"
#include "./state-measure/measure.h"
#include "./state-measure/scalar.h"
#include "./state-measure/strain.h"
#include "./state-measure/stress.h"
#include "./models/elasticity.h"
#include "./models/plasticity.h"

const size_t DIM = 3;
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> unidistr = std::uniform_real_distribution<double>((double)0, (double)1);

numerical_schema::type_schema numerical_schema::DEFAULT_NUMERICAL_SCHEMA = numerical_schema::type_schema::RATE_CALCULATE;
int main() 
{
#ifdef _DEBUG
	std::cout << " ------------------------------ Running in DEBUG mode --------------------------- \n";
#endif

	const auto gl = GLOBAL_BASIS<double>;
	//try 
	{
		using namespace measure;
		std::unordered_map<int, void*> map;
		map.insert({ 1, new int[10] });
		map.insert({ 2, new float[10] });
		map.insert({ 3, new double[10] });
		map.insert({ 4, new std::string[10] });
		map.insert({ 5, new std::array<double, 10>() });
		const auto type = typeid(std::array<double, 10>).name();
		auto xmap = static_cast<std::array<double, 10>*>(map[5]);
		using namespace tens;
		const auto qwew = container<double>(30, 0);
		//func(qwew, sqrt);
		//run_test();
		//const auto yy = container<double>(30, 2, std::move(std::unique_ptr<double>(arr)));
		const auto xx = container<double>(30, 0);
		auto scalar = container<double>(1, 1, 0.4534535);
		auto scalar_array = container<double>(10, 2, 0.4534535);
		auto scalar_array1 = object<double>(10, 0, FILL_TYPE::RANDOM);
		auto scalar_array2 = object<double>(10, 0, FILL_TYPE::RANDOM);
		scalar_array2 += scalar_array1;
		const auto scal(scalar_array2);
		double value = scalar;

		auto object = create_basis<double, 3>();
		const auto b1 = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
		const auto b2 = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
		const auto m1 = Matrix<double, 3>(FILL_TYPE::RANDOM);
		const auto m2 = Matrix<double, 3>(FILL_TYPE::RANDOM);
		auto v6 = Vector<double, 3>(std::array<double, 3>{1, 0, 0}, object);
		auto a3 = std::array<double, 3>{1, 0, 0};
		auto v1 = Vector<double, 3>(a3, b1);
		auto v2 = Vector<double, 3>(a3, b1);
		auto v3 = v1;
		auto v5 = std::move(v1);
		auto v4 = Vector<double, 3>(a3, b2);
		auto t1 = Tensor<double, 3>(m1, b2);
		auto t2 = Tensor<double, 3>(m2, b2);
		auto t3 = Tensor<double, 3>(m1, b1);

		// t2 = t1 / 1e-18;
		auto res = m1 * m1;
		auto vres = v4 + v3;
		v1 = v4 - v3;
		v1 = v4 + v3;
		auto vres1 = v1 * v3;
		auto tvres2 = v4 * t1;
		t1 = t2 * t2;
		t1 = t2 * 2.0;
		
		{
			auto Q = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
			auto elasticity = model::Elasticity<strain::GradDeform>(std::move(Q), numerical_schema::type_schema::RATE_CALCULATE);
			std::cout << elasticity;
			for (size_t i = 0; i < 1000; i++) {
				elasticity.step(1e-4);
			}
			std::cout << "Result of hypoelasic model \n";
			std::cout << elasticity;
		}
		{
			auto Q = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
			auto elasticity = model::Elasticity<strain::GradDeform>(std::move(Q), numerical_schema::type_schema::FINITE_CALCULATE);
			std::cout << elasticity;
			for (size_t i = 0; i < 10000; i++) {
				elasticity.step(1e-5);
			}
			std::cout << "Result of hyperelasic model \n";
			std::cout << elasticity;
		}
		{
			auto Q = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
			auto plasticity = model::Plasticity<strain::GradDeform, stress::CaushyStress, strain::GradDeformInelast, strain::GradDeformElast>(std::move(Q), numerical_schema::type_schema::RATE_CALCULATE);
			std::cout << plasticity;
			for (size_t i = 0; i < 10000; i++) {
				plasticity.step(1e-5);
			}
			std::cout << "Result of hypoplastic model \n";
			std::cout << plasticity;
		}
	}
	return 0;
}
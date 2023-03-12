#pragma once
#include "tensor/test/expect.h"
#include "tensor/test/test.h"
#include "tensor/object.h"
#include "tensor/quat.h"
#include "../state-measure/state_measure.h"

const size_t DIM = 3;
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> unidistr = std::uniform_real_distribution<double>((double)0, (double)1);

int main() 
{
#ifdef _DEBUG
	std::cout << " ------------------------------ Running in DEBUG mode --------------------------- \n";
#endif

	const auto gl = GLOBAL_BASIS<double>;
	//try 
	{
		run_test();
		using namespace tens;
		const auto xx = container<double>(30, 0);
		const auto scalar = container<double>(1, 1, 0.4534535);
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
		double vresd = vres1;
		auto tvres2 = v4 * t1;
		t1 = t2 * t2;
		t1 = t2 * 2.0;

		auto mes1 = state::Measure(t1, "F");
		auto mes2 = state::Measure(std::move(t3), "S");
		auto mes3 = mes1;

		//
		auto B = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
		auto state = state::create(std::move(B));// State<double>(b1);
		auto& pmes = mes1;
		state::insert(state, pmes);
		state::insert(state, mes2);
		//std::cout << (*state)["F"];
		auto state1 = state;
		std::cout << *state;
		return 0;
	}
}
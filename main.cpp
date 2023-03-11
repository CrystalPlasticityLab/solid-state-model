#pragma once
#include "tensor/test/expect.h"
#include "tensor/test/test.h"
#include "tensor/object.h"
#include "tensor/quat.h"
#include "factory/factory.h"
#include "../state-measure/state_measure.h"

const size_t DIM = 3;
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> unidistr = std::uniform_real_distribution<double>((double)0, (double)1);

int main() 
{
	try {
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
		auto res = m1 * m1;
		auto vres = v4 + v3;
		v1 = v4 - v3;
		v1 = v4 + v3;
		double vres1 = v1 * v3;
		//auto tvres2 = v4 * t1;
		t1 = t2 * t2;
		t1 = t2 * 2.0;

		auto state1 = StateMeasure(t1);
		auto state2 = StateMeasure(std::move(t1));
		return 0;
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
/*	try
	{
		//run_test();
		StateMeasure sm(t2);

		auto pbas = t2.get_basis_ref();
		//bas = m2;
		//auto vres2 = v4 * v3;
		auto vc = get_comp(v3);
		auto tc = get_comp(t1);
		v3 += v4;
		t2 = std::move(t1); 
		auto s = factory::state(b1);
		auto& ob1 = s.push<1>("vector");
		auto& ob2 = s.push<2>("Indent", FILL_TYPE::INDENT);
		auto& ob3 = s.push<2>("Tens_rand", FILL_TYPE::RANDOM);
		//s.remove("Tens_rand");
		//auto& ob4 = s.push<2>("Tens_rand", FILL_TYPE::RANDOM);
		//auto ob31 = s["Tens_rand"];
		auto& ob31 = s.get<2>("Tens_rand");
		//ob31 + ob31;
		//std::shared_ptr<basis_base<double, 3>> ob = *ob31;
		std::cout << ob3;
		std::cout << ob31;

		std::vector<boost::any> vect;
		vect.push_back(ob1);
		vect.push_back(ob3);
		vect.push_back(b1);
		std::cout << vect[0].type().name();
		return 0;
/*		std::shared_ptr<object<double, 3>> tp;
		{
			auto a = new array<double, 3>(FILL_TYPE::RANDOM);
			auto s = new factory::state<double, 3>(b);
			s->push("matrix", &m);
			s->push("vector", a);
			t = s->get<tensor<double, 3>>("matrix");
			tp = (*s)["matrix"];// s->get<tensor<double, 3>>("matrix");
			delete s;
			m = m;
		}
		//for (size_t i = 0; i < 100; i++)
		{
			//auto m = generate_rand<double, 3>();
			auto a = new array<double, 3>(FILL_TYPE::RANDOM);
			auto s = new factory::state<double, 3>(b);
			s->push("matrix", &m);
			s->push("vector", a);
			auto s_copy = *s;
			std::swap(s_copy, *s);
			auto s_move = std::move(s_copy);
			s_copy = std::move(s_move);
			auto t = s->get<tensor<double,3>>("matrix");
			//auto mq = s["matrix"];

			std::cout << t;
			//mp.move_to_basis(b1);
			std::cout << t;
			delete s;
			delete a;
		}
		
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
	*/
	return 0;
}
#pragma once
#include "tensor/test/expect.h"
#include "tensor/test/test.h"
#include "tensor/basis.h"
#include "tensor/quat.h"
#include "factory/factory.h"

const size_t DIM = 3;
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> unidistr = std::uniform_real_distribution<double>((double)0, (double)1);

int main() 
{
	using namespace tens;
	try
	{
		run_test();
		auto basis = create_basis<double, 3>();
		auto v6 = Vector<double, 3>(std::array<double, 3>{1, 0, 0}, basis);
		auto a3 = tens::container_array<double, 3>(0);
		auto m1 = Basis<double, 3>(0);
		auto m2 = Basis<double, 3>(0);
		auto pm1 = std::make_shared<Basis<double, 3>>(m1);
		//auto y = a3 * m1;
		auto v1 = Vector<double, 3>(a3, m1);
		auto v2 = Vector<double, 3>(a3, pm1);
		auto v3 = v1;
		auto v5 = std::move(v1);
		auto v4 = Vector<double, 3>(a3, m1);
		auto t1 = Tensor<double, 3>(m1, m2);
		auto t2 = Tensor<double, 3>(m1, pm1);
		auto res = m1 * m1;
		auto m1t = transpose(m1);
		auto vres = v4 + v3;
		v1 = v4 - v3;
		v1 = v4 + v3;
		auto tvres2 = v4 * t1;
		t1 = t2 * t2;
		t1 = t2 * 2.0;
		auto pbas = t2.get_basis_ref();
		//bas = m2;
		//auto vres2 = v4 * v3;
		auto vc = get_comp(v3);
		auto tc = get_comp(t1);
		v3 += v4;
		t2 = std::move(t1); 

/*		std::shared_ptr<basis<double, 3>> tp;
		{
			auto a = new array<double, 3>(ARRAY_TYPE::RANDOM);
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
			auto a = new array<double, 3>(ARRAY_TYPE::RANDOM);
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
		*/
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
	
	return 0;
}
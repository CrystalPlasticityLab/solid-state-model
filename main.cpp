#pragma once
#include "tensor/test/expect.h"
#include "tensor/test/test.h"
#include "tensor/basis.h"
#include "tensor/vector.h"
#include "tensor/tensor.h"
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
		auto basis = create_basis<double, 3>();
		auto c = tens::container<double, 3>(2);
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
		//double x = a1 * a2;
		//auto c3 = c1 += c2 * c1;
		auto cb = std::move(c);
		//std::cout << cb[8];
		c = std::move(cb);
		//const auto b = basis<double, 3>::create_basis_random();
		//const auto b1 = basis<double, 3>::create_basis_random();
		auto b = generate_rand_ort<double,3>();

		auto m = matrix<double, 3>::generate_rand();
		auto o1 = object<double, 3>(m, b);
		auto o2 = object<double, 3>(m, b);

		tensor<double, 3> t(m,b);
/*		std::shared_ptr<basis<double, 3>> tp;
		{
			auto a = new array<double, 3>(ARRAYTTYPE::RANDOM);
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
			auto a = new array<double, 3>(ARRAYTTYPE::RANDOM);
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
		run_test();
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
	
	return 0;
}
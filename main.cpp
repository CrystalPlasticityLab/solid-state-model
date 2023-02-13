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
	try
	{
		const auto b = create_basis(generate_rand_ort<double, 3>());
		const auto b1 = create_basis(generate_rand_ort<double, 3>());
		//for (size_t i = 0; i < 100; i++)
		{
			auto s = factory::state<double, 3>(b);
			auto m = generate_rand<double, 3>();
			auto a = array<double, 3>(ARRAYTTYPE::RANDOM);
			s.push("matrix", m);
			s.push("vector", a);
			auto mp = s.get<tensor<double,3>>("matrix");

			std::cout << mp;
			mp.move_to_basis(b1);
			std::cout << mp;
		}
		
		run_test();
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
	
	return 0;
}
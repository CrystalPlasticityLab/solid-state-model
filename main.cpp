#include "tensor/test/expect.h"
#include "tensor/test/test.h"
#include "tensor/basis.h"
#include "tensor/vector.h"
#include "tensor/tensor.h"
#include "tensor/quat.h"
#include "state/state.h"

const size_t DIM = 3;
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> unidistr = std::uniform_real_distribution<double>((double)0, (double)1);

int main() 
{
	try
	{
		const auto b = create_basis(generate_rand_ort<double, 3>());
		auto s = state::state<double, 3>(b);
		for (size_t i = 0; i < 100; i++)
		{
			auto m = generate_rand<double, 3>();
			auto a = array<double, 3>(ARRAYTTYPE::RANDOM);
			auto t = s.add_object(m);
			auto v = s.add_object(a);

			const auto x = 0.1;
		}
		
		run_test();
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
	
	return 0;
}
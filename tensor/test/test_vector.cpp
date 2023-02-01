#include "test.h"

void test_vector(){
    using namespace tens;

    std::cout << " =================== Start testing tens::vector ===================" << std::endl;
    int all_tests = 0;
    int pass_tests = 0;
    auto a_zero = array<double, 3>(0.0);
    auto a1 = random_array<double, 3>();
    auto a2 = random_array<double, 3>();
    auto a3 = random_array<double, 3>();
    auto a4 = random_array<double, 3>();

	auto basis1 = shared_handler_basis<double, 3>(generate_rand_ort<double, 3>());
	auto basis2 = shared_handler_basis<double, 3>(generate_rand_ort<double, 3>());
    auto v_zero1 = Vector<double, 3>(a_zero, basis1);
    auto v_zero2 = Vector<double, 3>(a_zero, basis2);
    auto v11 = Vector<double, 3>(a1, basis1);
    auto v12 = Vector<double, 3>(a2, basis1);
    auto v21 = Vector<double, 3>(a3, basis2);
    auto v22 = Vector<double, 3>(a4, basis2);
    {
        auto vr = tens::array<double, 3>(std::array<double, 3>{0,0,0});
        pass_tests += expect((v_zero1==v_zero2), "equal vectors");
        all_tests++;
    }
    
    std::cout << " Test passed : " << std::to_string(pass_tests) << "/" << std::to_string(all_tests) << std::endl;
    std::cout << " ==================== End Testing tens::vector ====================" << std::endl;
}
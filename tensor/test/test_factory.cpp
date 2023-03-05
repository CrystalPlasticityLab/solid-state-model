#include "test.h"

void test_factory() {
    using namespace tens;
    using namespace factory;

    std::cout << " =================== Start testing Factory ===================" << std::endl;
    int all_tests = 0;
    int pass_tests = 0;


    std::cout << " Test passed : " << std::to_string(pass_tests) << "/" << std::to_string(all_tests) << std::endl;
    std::cout << " ==================== End Testing Tensor ====================" << std::endl;
}
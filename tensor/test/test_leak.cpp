#include "test.h"

void test_leak(){
    using namespace tens;
    for (size_t i = 0; i < 1e6; i++)
    {
        {
            auto b1 = matrix<double, 3>::generate_rand();
            auto b2 = matrix<double, 3>::generate_rand();
            auto f1 = factory::state<double, 3>(b1);
            auto f2 = factory::state<double, 3>(b1);
            {
                auto m1 = matrix<double, 3>::generate_rand();
                auto m2 = matrix<double, 3>::generate_rand();

                auto t11 = tensor<double, 3>(m1, b1);
                auto t12 = tensor<double, 3>(m2, b2);
                auto t21 = t11;
                auto t31 = std::move(t21);
                t11 = std::move(t12);
            }
            {
                auto a1 = array<double, 3>(ARRAYTTYPE::RANDOM);
                auto a2 = array<double, 3>(ARRAYTTYPE::RANDOM);
                auto v11 = vector<double, 3>(a1, b1);
                auto v12 = vector<double, 3>(a2, b1);
                auto v21 = v11;
                auto v31 = std::move(v21);
                v11 = std::move(v12);
            }
        }
    }
}
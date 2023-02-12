#include "test.h"

void test_leak(){
    using namespace tens;
    for (size_t i = 0; i < 1e6; i++)
    {
        {
            auto basis1 = create_basis(generate_rand_ort<double, 3>());
            auto basis2 = create_basis(generate_rand_ort<double, 3>());
            {
                auto m1 = generate_rand<double, 3>();
                auto m2 = generate_rand<double, 3>();

                auto t11 = Tensor<double, 3>(m1, basis1);
                auto t12 = Tensor<double, 3>(m2, basis1);
                auto t21 = t11;
                auto t31 = std::move(t21);
                t11 = std::move(t12);
            }
            {
                auto a1 = array<double, 3>(ARRAYTTYPE::RANDOM);
                auto a2 = array<double, 3>(ARRAYTTYPE::RANDOM);
                auto v11 = vector<double, 3>(a1, basis1);
                auto v12 = vector<double, 3>(a2, basis1);
                auto v21 = v11;
                auto v31 = std::move(v21);
                v11 = std::move(v12);
            }
        }
    }
}
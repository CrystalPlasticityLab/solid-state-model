#include "test.h"

void test_leak(){
/*    using namespace tens;
    for (size_t i = 0; i < 1e6; i++)
    {
        {
            auto b1 = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
            auto b2 = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
            //auto f1 = factory::state<double, 3>(b1);
            //auto f2 = factory::state<double, 3>(b1);
            {
                auto m1 = Matrix<double, 3>(); m1.fill_rand();
                auto m2 = Matrix<double, 3>(); m2.fill_rand();

                auto t11 = Tensor<double, 3>(m1, b1);
                auto t12 = Tensor<double, 3>(m2, b2);
                auto t21 = t11;
                auto t31 = std::move(t21);
                t11 = std::move(t12);
            }
            {
                auto a1 = Array<double, 3>(FILL_TYPE::RANDOM);
                auto a2 = Array<double, 3>(FILL_TYPE::RANDOM);
                auto v11 = Vector<double, 3>(a1, b1);
                auto v12 = Vector<double, 3>(a2, b1);
                auto v21 = v11;
                auto v31 = std::move(v21);
                v11 = std::move(v12);
            }
        }
    }
    */
}
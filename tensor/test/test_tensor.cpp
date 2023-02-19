#include "test.h"

void test_tensor(){
    using namespace tens;

    std::cout << " =================== Start testing tensor ===================" << std::endl;
    int all_tests = 0;
    int pass_tests = 0;
    auto m_zero = matrix<double, 3>(MATRIXINITTYPE::ZERO);
    const auto m1 = matrix<double, 3>::generate_rand();
    const auto m2 = matrix<double, 3>::generate_rand();
    const auto m3 = matrix<double, 3>::generate_rand();
    const auto m4 = matrix<double, 3>::generate_rand();

	const auto basis1 = generate_rand_ort<double, 3>();
	const auto basis2 = generate_rand_ort<double, 3>();
    const auto t_zero1 = tensor<double, 3>(m_zero, basis1);
    const auto t_zero2 = tensor<double, 3>(m_zero, basis2);
    const auto t11 = tensor<double, 3>(m1, basis1);
    const auto t12 = tensor<double, 3>(m2, basis1);
    const auto t21 = tensor<double, 3>(m1, basis2);
    const auto t22 = tensor<double, 3>(m2, basis2);

    {
        pass_tests += expect((t_zero1==t_zero2), "equal zero tensors");
        all_tests++;
    }
    {
        pass_tests += expect((t11.get_comp_at_basis(t11) == m1), "equal tensor components");
        all_tests++;
    }
    {
        pass_tests += expect((t21+t_zero1 == t21) && (t22+t_zero2 == t22), "T+T_zero = T");
        all_tests++;
    }
    {
        auto tr = t11 + t12;
        auto mr = m1 + m2;
        pass_tests += expect((tr.get_comp_at_basis(t11) == mr), "check sum at different basis");
        all_tests++;
    }
    {
        auto tr = t11;
        tr += t21;
        pass_tests += expect((t11+t21 == tr), "check += at different basis");
        all_tests++;
    }
   {
       auto tr = t12;
       tr -= t22;
       pass_tests += expect((t12-t22 == tr), "check -= at different basis");
       all_tests++;
   }
   t11* t22;
   {
       auto tr = t11 - t12;
       auto mr = m1 - m2;
       pass_tests += expect((tr.get_comp_at_basis(t11) == mr), "check sub at different basis");
       all_tests++;
   }
   {
       auto tr1 = t11;
       //tr1.change_basis(basis2);
       pass_tests += expect((t11 == tr1), "tensor has not changed after changing basis");
       all_tests++;
   }
   {
       auto tr1 = t11;
       //tr1.change_basis(basis2);
       pass_tests += expect((m1 == tr1.get_comp_at_basis(t11)), "component has not changed at the same basis");
       all_tests++;
   }
   {
       auto tr1 = t11;
       //tr1.move_to_basis(basis2);
       pass_tests += expect(!(t11 == tr1), "tensor has not changed after changing basis");
       all_tests++;
       //tr1.move_to_basis(basis1);
       pass_tests += expect((t11 == tr1), "tensor has changed after changing basis to the origin basis");
       all_tests++;
   }
    
    std::cout << " Test passed : " << std::to_string(pass_tests) << "/" << std::to_string(all_tests) << std::endl;
    std::cout << " ==================== End Testing tensor ====================" << std::endl;
}
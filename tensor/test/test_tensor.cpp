#include "test.h"

double forward(double x) {
    return x * x * x;
};
double backward(double x) {
    if (x < 0)
        return -pow(-x, 0.333333333333333333);
    else
        return pow(x, 0.333333333333333333);
};

void test_tensor(){
    using namespace tens;

    std::cout << " =================== Start testing Tensor ===================" << std::endl;
    const auto gl = GLOBAL_BASIS<double, 3>;
    int all_tests = 0;
    int pass_tests = 0;
    const auto m_I = Matrix<double, 3>(FILL_TYPE::INDENT);
    const auto m_zero = Matrix<double, 3>(FILL_TYPE::ZERO);
    const auto m1 = Matrix<double, 3>(FILL_TYPE::RANDOM);
    const auto m2 = Matrix<double, 3>(FILL_TYPE::RANDOM);
    const auto m3 = Matrix<double, 3>(FILL_TYPE::RANDOM);
    const auto m4 = Matrix<double, 3>(FILL_TYPE::RANDOM);
    const auto a1 = Array<double, 3>(FILL_TYPE::RANDOM);

    const auto basis1 = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
    const auto basis2 = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
    const auto t_zero1 = Tensor<double, 3>(m_zero, basis1);
    const auto t_zero2 = Tensor<double, 3>(m_zero, basis2);
    const auto t11 = Tensor<double, 3>(m1, basis1);
    const auto t12 = Tensor<double, 3>(m2, basis1);
    const auto t21 = Tensor<double, 3>(m1, basis2);
    const auto t22 = Tensor<double, 3>(m2, basis2);
    const auto v11 = Vector<double, 3>(a1, basis1);
    const auto v21 = Vector<double, 3>(a1, basis2);

    {
        pass_tests += expect((m1*m1.inverse() == m_I) && (m1.inverse() * m1 == m_I), "check inverse func m*m^-1 = m^-1*m = I");
        all_tests++;
    }
    {
        pass_tests += expect((t_zero1==t_zero2), "equal zero Tensors");
        all_tests++;
    }
    {
        pass_tests += expect((t11.get_comp_at_basis(t11) == m1), "equal Tensor components");
        all_tests++;
    }
    {
        pass_tests += expect((t21+t_zero1 == t21) && (t22+t_zero2 == t22), "T+T_zero = T");
        all_tests++;
    }
    {
        auto tr = t11;
        pass_tests += expect((tr == t11), "check t->t1, t1 == t");
        all_tests++;
    }
    {
        auto tr = t11 + t11;
        pass_tests += expect((tr == 2.0 * t11), "check t+t = 2*t");
        all_tests++;
    }
    {
        auto t1l = t11;
        t1l.change_basis(t22);
        pass_tests += expect((t1l == t11), "check t = t after changing basis");
        all_tests++;
    }
    {
        auto t1l = t11;
        t1l.change_basis(t22);
        t1l.change_basis(t11);
        pass_tests += expect((t1l.get_comp_at_basis(t11) == m1), "check t.comp = t.comp after changing basis");
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
    {
        auto tr = t11 - t12;
        auto mr = m1 - m2;
        pass_tests += expect((tr.get_comp_at_basis(t11) == mr), "check sub at different basis");
        all_tests++;
    }
    {
        auto tr1 = t11;
        tr1.change_basis(basis2);
        pass_tests += expect((t11 == tr1), "Tensor has not changed after changing basis");
        all_tests++;
    }
    {
        auto tr1 = t11;
        tr1.change_basis(basis2);
        pass_tests += expect((m1 == tr1.get_comp_at_basis(t11)), "component has not changed at the same basis");
        all_tests++;
    }
    {
        auto tv1 = t11 * v11;
        auto tv2 = t21 * v21;
        pass_tests += expect((tv1.get_comp_at_basis(basis1) == tv2.get_comp_at_basis(basis2)), "check scal t.v");
        all_tests++;
    }
    {
        auto vt1 = v11 * t11;
        auto vt2 = v21 * t21;
        pass_tests += expect((vt1.get_comp_at_basis(basis1) == vt2.get_comp_at_basis(basis2)), "check scal v*t");
        all_tests++;
    }
    {
        auto vt1 = v11 * t21;
        auto vt2 = transpose(t21) * v11;
        pass_tests += expect((vt1 == vt2), "check scal v*t = tT*v");
        all_tests++;
    }
    {
        auto vt1 = t11 * v21;
        auto vt2 = v21 * transpose(t11);
        pass_tests += expect((vt1 == vt2), "check scal t*v = v*tT");
        all_tests++;
    }
    {
        auto vt1 = v11 * t11;
        auto vt2 = t11 * v11;
        pass_tests += expect(!(vt1 == vt2), "check scal v*t != t*v");
        all_tests++;
    }
    {
        const auto sm = container<double, 3, 2>(FILL_TYPE::RANDOMSYMM);
        const auto ref_obj = object<double, 3, 2>(sm, gl);
        const auto eig_obj = eigen_object(sm);
        pass_tests += expect((ref_obj == eig_obj), "eigen test");
        all_tests++;
    }
    {
        auto I2 = IDENT_MATRIX<double, 3> * 2.0;

        auto I2sqrt = func(I2, sqrt) - (IDENT_MATRIX<double, 3> * sqrt(2.0));
        auto err = I2sqrt.trace();
        pass_tests += expect(math::is_small_value(err), "func of matrix test #1");
        all_tests++;
    }
    {
        // TODO: check accuracy of eigen
        // auto m1sqrt = func(m1, forward);
        // auto m1log = func(m1sqrt, backward);
        // auto err = (m1- m1log);
        // std::cout << m1;
        // std::cout << m1sqrt;
        // std::cout << m1log;
        // pass_tests += expect(err == m_zero, "func of matrix test #2");
        // all_tests++;
    }
    
   //{
   //    auto tr1 = t11;
   //    //tr1.move_to_basis(basis2);
   //    pass_tests += expect(!(t11 == tr1), "Tensor has not changed after changing basis");
   //    all_tests++;
   //    //tr1.move_to_basis(basis1);
   //    pass_tests += expect((t11 == tr1), "Tensor has changed after changing basis to the origin basis");
   //    all_tests++;
   //}
    
    std::cout << " Test passed : " << std::to_string(pass_tests) << "/" << std::to_string(all_tests) << std::endl;
    std::cout << " ==================== End Testing Tensor ====================" << std::endl;
}
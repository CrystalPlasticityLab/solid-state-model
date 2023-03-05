#include "test.h"

void test_vector(){
    using namespace tens;

    std::cout << " =================== Start testing Vector ===================" << std::endl;
    int all_tests = 0;
    int pass_tests = 0;
    auto a_zero = Array<double, 3>(0.0);
    const auto a1 = Array<double, 3>(FILL_TYPE::RANDOM);
    const auto a2 = Array<double, 3>(FILL_TYPE::RANDOM);
    const auto a3 = Array<double, 3>(FILL_TYPE::RANDOM);
    const auto a4 = Array<double, 3>(FILL_TYPE::RANDOM);

	const auto basis1 = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
	const auto basis2 = create_basis<double, 3>(DEFAULT_ORTH_BASIS::RANDOM);
    const auto v_zero1 = Vector<double, 3>(a_zero, basis1);
    const auto v_zero2 = Vector<double, 3>(a_zero, basis2);
    const auto v11 = Vector<double, 3>(a1, basis1);
    const auto v12 = Vector<double, 3>(a2, basis1);
    const auto v21 = Vector<double, 3>(a1, basis2);
    const auto v22 = Vector<double, 3>(a2, basis2);

    {
        pass_tests += expect((v_zero1==v_zero2), "equal zero vectors");
        all_tests++;
    }
    {
        pass_tests += expect((v11.get_comp_at_basis(v11) == a1), "equal Vector components");
        all_tests++;
    }
    {
        pass_tests += expect((v21+v_zero1 == v21) && (v22+v_zero2 == v22), "v+v_zero = v");
        all_tests++;
    }
    {
        auto vr = v11;
        pass_tests += expect((vr == v11), "check v->v1, v1 == v");
        all_tests++;
    }
    {
        auto vr = v11 + v11;
        pass_tests += expect((vr == 2.0*v11), "check v+v = 2*v");
        all_tests++;
    }
    {
        auto v1l = v11;
        v1l.change_basis(v22);
        pass_tests += expect((v1l == v11), "check v = v after changing basis");
        all_tests++;
    }
    {
        auto v1l = v11;
        v1l.change_basis(v22);
        v1l.change_basis(v11);
        pass_tests += expect((v1l.get_comp_at_basis(v11) == a1), "check v.comp = v.comp after changing basis");
        all_tests++;
    }
    {
        auto v1l = v11;
        v1l.change_basis(v22);
        auto vr = v11 + v1l;
        pass_tests += expect((vr == 2.0 * v11), "check v+v = 2*v");
        all_tests++;
    }
    {
        auto vr = v11 + v12;
        auto ar = a1 + a2;
        pass_tests += expect((vr.get_comp_at_basis(v11) == ar), "check sum at different basis");
        all_tests++;
    }
    {
        auto vr = v11;
        vr += v21;
        auto vl = v11 + v21;
        pass_tests += expect((vl == vr), "check += at different basis");
        all_tests++;
    }
    {
        auto vr = v12;
        vr -= v22;
        auto vl = v12 - v22;
        pass_tests += expect((vl == vr), "check -= at different basis");
        all_tests++;
    }
    {
        auto vr = v11 - v12;
        auto ar = a1 - a2;
        pass_tests += expect((vr.get_comp_at_basis(v11) == ar), "check sub at different basis");
        all_tests++;
    }
    {
        Vector<double,3> vr1 = v11;
        vr1.change_basis(v21);
        pass_tests += expect((v11 == vr1), "Vector has not changed after changing basis");
        all_tests++;
    }
    {
        Vector<double,3> vr1 = v11;
        vr1.change_basis(v21);
        pass_tests += expect((a1 == vr1.get_comp_at_basis(v11)), "component has not changed at the same basis");
        all_tests++;
    }
  /* {
        Vector<double,3> vr1 = v11;
        vr1.move_to_basis(v21);
        pass_tests += expect(!(v11 == vr1), "Vector has changed after changing basis");
        all_tests++;
        vr1.move_to_basis(v11);
        pass_tests += expect((v11 == vr1), "Vector has not changed after changing basis to the origin basis");
        all_tests++;
    }
    {
        pass_tests += expect((v11*v12 == v21*v22), "scalar product is invariant");
        all_tests++;
    }
    {
        const auto cross1 =  vector_product(v11,v21);
        const auto cross2 = -vector_product(v21,v11);
        pass_tests += expect((cross1==cross2), "v2 x v1 = -(v1 x v2)");
        all_tests++;
    } 
    {
        const auto cross = vector_product(v11,v21);
        pass_tests += expect(( !is_not_small_value(cross*v11) && !is_not_small_value(cross*v21)), "v1 x v2 ia ort to v1 and v2");
        all_tests++;
    }
    */
    std::cout << " Test passed : " << std::to_string(pass_tests) << "/" << std::to_string(all_tests) << std::endl;
    std::cout << " ==================== End Testing Vector ====================" << std::endl;
}
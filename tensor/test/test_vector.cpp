#include "test.h"

void test_vector(){
    using namespace tens;

    std::cout << " =================== Start testing vector ===================" << std::endl;
    int all_tests = 0;
    int pass_tests = 0;
    auto a_zero = array<double, 3>(0.0);
    const auto a1 = array<double, 3>(ARRAYTTYPE::RANDOM);
    const auto a2 = array<double, 3>(ARRAYTTYPE::RANDOM);
    const auto a3 = array<double, 3>(ARRAYTTYPE::RANDOM);
    const auto a4 = array<double, 3>(ARRAYTTYPE::RANDOM);

	const auto basis1 = generate_rand_ort<double, 3>();
	const auto basis2 = generate_rand_ort<double, 3>();
    const auto v_zero1 = vector<double, 3>(a_zero, basis1);
    const auto v_zero2 = vector<double, 3>(a_zero, basis2);
    const auto v11 = vector<double, 3>(a1, basis1);
    const auto v12 = vector<double, 3>(a2, basis1);
    const auto v21 = vector<double, 3>(a1, basis2);
    const auto v22 = vector<double, 3>(a2, basis2);

    {
        auto vr = array<double, 3>(std::array<double, 3>{0,0,0});
        pass_tests += expect((v_zero1==v_zero2), "equal zero vectors");
        all_tests++;
    }
    {
        pass_tests += expect((v11.get_comp_at_basis(v11) == a1), "equal vector components");
        all_tests++;
    }
    {
        pass_tests += expect((v21+v_zero1 == v21) && (v22+v_zero2 == v22), "v+v_zero = v");
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
        pass_tests += expect((v11+v21 == vr), "check += at different basis");
        all_tests++;
    }
    {
        auto vr = v12;
        vr -= v22;
        pass_tests += expect((v12-v22 == vr), "check -= at different basis");
        all_tests++;
    }
    {
        auto vr = v11 - v12;
        auto ar = a1 - a2;
        pass_tests += expect((vr.get_comp_at_basis(v11) == ar), "check sub at different basis");
        all_tests++;
    }
    {
        vector<double,3> vr1 = v11;
        vr1.change_basis(v21);
        pass_tests += expect((v11 == vr1), "vector has not changed after changing basis");
        all_tests++;
    }
    {
        vector<double,3> vr1 = v11;
        vr1.change_basis(v21);
        pass_tests += expect((a1 == vr1.get_comp_at_basis(v11)), "component has not changed at the same basis");
        all_tests++;
    }
    {
        vector<double,3> vr1 = v11;
        vr1.move_to_basis(v21);
        pass_tests += expect(!(v11 == vr1), "vector has changed after changing basis");
        all_tests++;
        vr1.move_to_basis(v11);
        pass_tests += expect((v11 == vr1), "vector has not changed after changing basis to the origin basis");
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
    
    std::cout << " Test passed : " << std::to_string(pass_tests) << "/" << std::to_string(all_tests) << std::endl;
    std::cout << " ==================== End Testing vector ====================" << std::endl;
}
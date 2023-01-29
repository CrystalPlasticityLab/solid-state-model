#include "test_vector.h"

void test_vector(){
    std::cout << " =================== Start testing vect_base ===================" << std::endl;
    int all_tests = 0;
    int pass_tests = 0;
    auto v_zero = vect_base<double, 3>(0.0);
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0,0,0});
        pass_tests += expect((v_zero==vr), "equal vectors");
        all_tests++;
    }
    {
        auto vl = vect_base<double, 3>(std::array<double, 3>{0.,0.1,0});
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,0,0});
        pass_tests += expect(!(vl==vr), "not equal vectors");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect(((vr+v_zero)==vr), "(v+zero_vector) = zero_vector");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect(((vr-v_zero)==vr), "(v-zero_vector) = zero_vector");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect(((vr-vr)==v_zero), "(v-v) = zero_vector (operator -)");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect(((vr-=vr)==v_zero), "(v-v) = 2*v (operator -=)");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect(((vr+vr)==vr*2), "(v+v) = 2*v (operator +)");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect(((vr+=vr)==vr), "(v+v) = 2*v (operator +=)");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect(((vr*1.0)==vr), "1*v = v");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect(((vr*0.0)==v_zero), "0*v = zero_vector");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        double x = 1.23456;
        pass_tests += expect(((vr*x).norm()==vr.norm()*x), "|x*v| = x*|v| (operator *)");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        double x = 1.23456;
        pass_tests += expect(((vr/x).norm()==vr.norm()/x), "|v/x| = |v|/x (operator /)");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        const double vr_norm = vr.norm();
        const double x = 1.23456;
        vr *= x;
        pass_tests += expect((vr.norm()==vr.norm())&&(vr.norm()==vr_norm*x), "|x*v| = x*|v| (operator *=)");
        all_tests++;
    }
    {
        auto arr = std::array<double, 3>{0.01,-1.234,3.853};
        auto vr = vect_base<double, 3>(arr);
        const double vr_norm = vr.norm();
        const double x = 1.23456;
        vr /= x;
        pass_tests += expect((vr.norm()==vr.norm())&&(vr.norm()==vr_norm/x), "|v/x| = |v|/x (operator /=)");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect((vr*v_zero)==0.0, "v.zero_vector = 0");
        all_tests++;
    }
    {
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        pass_tests += expect(sqrt(vr*vr)==vr.norm(), "v.v = |v|^2");
        all_tests++;
    }

    // TODO: check copy and move ctors and operators =

    {
        auto vl = vect_base<double, 3>(std::array<double, 3>{123.123,7.234,-45.2356});
        auto vr = vect_base<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        auto vl_add_vr = vect_base<double, 3>(std::array<double, 3>{123.133, 6., -41.3826});
        auto vl_sub_vr = vect_base<double, 3>(std::array<double, 3>{123.113, 8.468, -49.0886});
        auto vr_sub_vl = vect_base<double, 3>(std::array<double, 3>{-123.113, -8.468, 49.0886});
        auto vl_cros_vr = vect_base<double, 3>(std::array<double, 3>{-27.948128399999995, -474.845275, -152.006122});
        double vl_scal_vr = -181.9882928;
        {
            pass_tests += expect(((vl+vr)==vl_add_vr), "v1+v2 check");
            all_tests++;
        }
        {
            pass_tests += expect(((vr+vl)==vl_add_vr), "v2+v1 check");
            all_tests++;
        }
        {
            pass_tests += expect(((vl-vr)==vl_sub_vr), "v1-v2 check");
            all_tests++;
        }
        {
            pass_tests += expect(((vr-vl)==vr_sub_vl), "v2-v1 check");
            all_tests++;
        }
        {
            pass_tests += expect((vr-vl)==-(vl-vr), "v2-v1 = v1-v2");
            all_tests++;
        }
        {
            pass_tests += expect((vl*vr==vl_scal_vr), "v1.v2 check");
            all_tests++;
        }
        {
            pass_tests += expect((vr*vl==vl_scal_vr), "v2.v1 check");
            all_tests++;
        }
        {
           pass_tests += expect((vector_product(vl,vr)==vl_cros_vr), "v2 x v1 check");
           all_tests++;
        }
        {
           pass_tests += expect((vector_product(vl,vr)==-vector_product(vr,vl)), "v2 x v1 = -(v1 x v2)");
           all_tests++;
        }
    }
    
    std::cout << " Test passed : " << std::to_string(pass_tests) << "/" << std::to_string(all_tests) << std::endl;
    std::cout << " ==================== End Testing vect_base ====================" << std::endl;
}
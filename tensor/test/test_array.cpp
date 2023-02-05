#include "test.h"

void test_array(){
    using namespace tens;

    std::cout << " =================== Start testing array ===================" << std::endl;
    int all_tests = 0;
    int pass_tests = 0;
    auto v_zero = tens::array<double, 3>(0.0);
    auto v1 = tens::array<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
    {
        auto vr = tens::array<double, 3>(std::array<double, 3>{0,0,0});
        pass_tests += expect((v_zero==vr), "equal vectors");
        all_tests++;
    }
    {   
        try{
            v_zero.normalize();
        } catch (std::exception &e){
            pass_tests += expect(std::string(e.what()) == ErrorMessage::DIV_BY_ZERO_MSG, "throw div zero at normalize");
        }
        all_tests++;
    }
    {
        auto vl = tens::array<double, 3>(std::array<double, 3>{0.,0.1,0});
        auto vr = tens::array<double, 3>(std::array<double, 3>{0.01,0,0});
        pass_tests += expect(!(vl==vr), "not equal vectors");
        all_tests++;
    }
    {
        pass_tests += expect(((v1+v_zero)==v1), "(v+zero_vector) = zero_vector");
        all_tests++;
    }
    {
        pass_tests += expect(((v1-v_zero)==v1), "(v-zero_vector) = zero_vector");
        all_tests++;
    }
    {
        pass_tests += expect(((v1-v1)==v_zero), "(v-v) = zero_vector (operator -)");
        all_tests++;
    }
    {
        pass_tests += expect(((v1-=v1)==v_zero), "(v-v) = 2*v (operator -=)");
        all_tests++;
    }
    {
        pass_tests += expect(((v1+v1)==v1*2), "(v+v) = 2*v (operator +)");
        all_tests++;
    }
    {
        pass_tests += expect(((v1+=v1)==v1), "(v+v) = 2*v (operator +=)");
        all_tests++;
    }
    {
        pass_tests += expect(((v1*1.0)==v1), "1*v = v");
        all_tests++;
    }
    {
        pass_tests += expect(((v1*0.0)==v_zero), "0*v = zero_vector");
        all_tests++;
    }
    {
        double x = 1.23456;
        pass_tests += expect(((v1*x).norm()==v1.norm()*x), "|x*v| = x*|v| (operator *)");
        all_tests++;
    }
    {
        double x = 1.23456;
        pass_tests += expect(((v1/x).norm()==v1.norm()/x), "|v/x| = |v|/x (operator /)");
        all_tests++;
    }
    {
        const double vr_norm = v1.norm();
        const double x = 1.23456;
        v1 *= x;
        pass_tests += expect((v1.norm()==v1.norm())&&(v1.norm()==vr_norm*x), "|x*v| = x*|v| (operator *=)");
        all_tests++;
    }
    {
        auto vr = tens::array<double, 3>(v1);
        const double vr_norm = vr.norm();
        const double x = 1.23456;
        vr /= x;
        pass_tests += expect((vr.norm()==vr.norm())&&(vr.norm()==vr_norm/x), "|v/x| = |v|/x (operator /=)");
        all_tests++;
    }
    {
        pass_tests += expect((v1*v_zero)==0.0, "v.zero_vector = 0");
        all_tests++;
    }
    {
        pass_tests += expect(sqrt(v1*v1)==v1.norm(), "v.v = |v|^2");
        all_tests++;
    }

    // TODO: check copy and move ctors and operators =

    {
        auto vl = tens::array<double, 3>(std::array<double, 3>{123.123,7.234,-45.2356});
        auto vr = tens::array<double, 3>(std::array<double, 3>{0.01,-1.234,3.853});
        auto vl_add_vr = tens::array<double, 3>(std::array<double, 3>{123.133, 6., -41.3826});
        auto vl_sub_vr = tens::array<double, 3>(std::array<double, 3>{123.113, 8.468, -49.0886});
        auto vr_sub_vl = tens::array<double, 3>(std::array<double, 3>{-123.113, -8.468, 49.0886});
        auto vl_cros_vr = tens::array<double, 3>(std::array<double, 3>{-27.948128399999995, -474.845275, -152.006122});
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
           pass_tests += expect((array_vector_product(vl,vr)==vl_cros_vr), "v2 x v1 check");
           all_tests++;
        }
        {
           pass_tests += expect((array_vector_product(vl,vr)==-array_vector_product(vr,vl)), "v2 x v1 = -(v1 x v2)");
           all_tests++;
        }
    }
    
    std::cout << " Test passed : " << std::to_string(pass_tests) << "/" << std::to_string(all_tests) << std::endl;
    std::cout << " ==================== End Testing array ====================" << std::endl;
}
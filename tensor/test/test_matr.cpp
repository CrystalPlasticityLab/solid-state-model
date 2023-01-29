#include "test_matr.h"

void test_matr(){
    std::cout << " =================== Start testing matr_base ===================" << std::endl;
    int all_tests = 0;
    int pass_tests = 0;
    auto m_zero = matrix_base<double, 3>(MATRIXINITTYPE::ZERO);
    auto m_indent = matrix_base<double, 3>(MATRIXINITTYPE::INDENT);
	auto m_rand = matrix_base<double, 3>(matrix_generator::generate_rand<double, 3>());
    
    {
        pass_tests += expect((m_indent==m_indent), "I = I");
        all_tests++;
    }
    {
        pass_tests += expect(!(m_indent==m_rand), "I != M_rand");
        all_tests++;
    }
    {
        pass_tests += expect((m_indent+m_zero==m_indent), "I+0 = 0");
        all_tests++;
    }
    {
        pass_tests += expect((m_rand-m_rand==m_zero), "M-M = 0");
        all_tests++;
    }
    {
        pass_tests += expect((m_rand+m_rand==m_rand*2), "M+M = 2*M");
        all_tests++;
    }
    {
        pass_tests += expect((m_rand==m_rand*1), "M = 1*M");
        all_tests++;
    }
    {
        m_rand*=0;
        pass_tests += expect((m_rand==m_zero), "0 = 0*M (operator *=)");
        all_tests++;
    }
    {
        m_rand*=1;
        pass_tests += expect((m_rand==m_rand), "M = 1*M (operator *=)");
        all_tests++;
    }
    {
        m_rand/=1;
        pass_tests += expect((m_rand==m_rand), "M = M/1 (operator /=)");
        all_tests++;
    }

    auto m1 = std::array<std::array<double, 3>, 3>{
        std::array<double, 3>{0.212852, -0.366117, 0.905898},
        std::array<double, 3>{0.97707, 0.0746825, -0.199392},
        std::array<double, 3>{0.00534623, 0.927567, 0.373618} };

    auto m2 = std::array<std::array<double, 3>, 3>{
        std::array<double, 3>{-0.498728, 0.834428, 0.234522},
        std::array<double, 3>{0.854559, 0.518604, -0.0279043},
        std::array<double, 3>{-0.144908, 0.186496, -0.97171}  };
    {
        auto ml = matrix_base<double, 3>(m1);
        auto mr = matrix_base<double, 3>(m2);
        auto ml_add_mr = std::array<std::array<double, 3>, 3>{
                         std::array<double, 3>{-0.285876, 0.468311, 1.14042}, 
                         std::array<double, 3>{1.831629,  0.5932865, -0.2272963}, 
                         std::array<double, 3>{-0.13956177,  1.114063, -0.598092}};
        {
            pass_tests += expect((ml+mr==matrix_base<double, 3>(ml_add_mr)), "M1+M2 check");
            all_tests++;
        }
        auto ml_sub_mr = std::array<std::array<double, 3>, 3>{
                         std::array<double, 3>{0.71158, -1.200545, 0.671376}, 
                         std::array<double, 3>{0.122511, -0.4439215, -0.1714877}, 
                         std::array<double, 3>{0.15025423, 0.741071, 1.345328}};
        {
            pass_tests += expect((ml-mr==matrix_base<double, 3>(ml_sub_mr)), "M1-M2 check");
            all_tests++;
        }
        {
            pass_tests += expect((mr-ml==-matrix_base<double, 3>(ml_sub_mr)), "M1-M2 = -(M2-M1)");
            all_tests++;
        }
        {
            pass_tests += expect(((ml-mr)+(ml+mr)==ml*2), "M1-M2 + (M2-M1) = 2*M1");
            all_tests++;
        }
        auto ml_mul_mr = std::array<std::array<double, 3>, 3>{
                         std::array<double, 3>{-0.550295697043, 0.156686281396, -0.8201354302328999}, 
                         std::array<double, 3>{-0.3945780685065, 0.8168393987579999, 0.42081164797525006}, 
                         std::array<double, 3>{0.7358541762135601, 0.5551792630024399, -0.38767764606604}};
        auto mr_mul_ml = std::array<std::array<double, 3>, 3>{
                         std::array<double, 3>{0.71039312225606, 0.46244483626, -0.530553324924}, 
                         std::array<double, 3>{0.6884578197422109, -0.3000210420111, 0.6603122514565999},
                         std::array<double, 3>{0.1461807039507, -0.834344859814, -0.5315060245959999}};
        {
            pass_tests += expect((ml*mr==matrix_base<double, 3>(ml_mul_mr)), "M1.M2 check");
            all_tests++;
        }
        {
            pass_tests += expect((mr*ml==matrix_base<double, 3>(mr_mul_ml)), "M2.M1 check");
            all_tests++;
        }
        {
            pass_tests += expect((mr.transpose().transpose()==mr), "(Mt)t = M");
            all_tests++;
        }
        {
            pass_tests += expect((mr.transpose()*ml.transpose()==matrix_base<double, 3>(ml_mul_mr).transpose()), "M2t.M1t = (M1.M2)t");
            all_tests++;
        }
        {
            pass_tests += expect(((mr.transpose()*ml.transpose()).transpose()==matrix_base<double, 3>(ml_mul_mr)), "(M2t.M1t)t = M1.M2");
            all_tests++;
        }
        {
            pass_tests += expect(!ml.check_ort(), "M rand is not orthohonal");
            all_tests++;
        }
	    auto m_rand_ort = matrix_base<double, 3>(matrix_generator::generate_rand_ort<double, 3>());{
            pass_tests += expect(m_rand_ort.check_ort(), "M ort is orthohonal");
            all_tests++;
        }
    }
    // TODO: m1:m2
    // TODO: check copy and move ctors and operators =
    // TODO: transform with ort M

    std::cout << " Test passed : " << std::to_string(pass_tests) << "/" << std::to_string(all_tests) << std::endl;
    std::cout << " ==================== End Testing matr_base ====================" << std::endl;
}
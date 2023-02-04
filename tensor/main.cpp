#include "test/expect.h"
#include "test/test.h"
#include "vector.h"
#include "tensor.h"
#include "quat.h"

const size_t DIM = 3;
std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
extern std::uniform_real_distribution<double> unidistr = std::uniform_real_distribution<double>((double)0, (double)1);


void quternion_unit_test()
{
	//vector<double, 3> v2( 1.0, GLOBAL_DEFAULT_BASIS<double, 3>);
	//vector<double, 3> v2( 2, GLOBAL_DEFAULT_BASIS<double, 3>);
	//vector<double, 3> v3(GLOBAL_DEFAULT_BASIS<double, 3>);
	//v1 *= 1.4;
	//v1 * 0.5;
	//v1 = -v2;
	//std::cout << v1;
	//auto x = v1* v2;
	//v3 = (v1 - v2*0.5) * 0.3;
	//vector<double, 3> v5 ( v1 + v2);
	//v5.normalize();
	////v5 = v1 + v2 * 0.5;
	//v3 = v1;
	//v3 = v1.vector_product(v2);
	//quat<double> q0;
	//quat<double> q1(0.1, v1);
	//quat<double> q2(-0.3, v2);
	//quat<double> q3(0.3, v2);
	//quat<double> q4(q3);
	//q2 = *q1;
	//q0 = q1;
	//q0 = q2 / 0.1;
	//q0 /= 0.1;
	//q3 = !q0;
	//q2 = q0 + q1;
	//q3 = q0 * q1;
	//q4 = q1 * 1.1 - q2;
	//q3 += q0 * 1.5;
	////vector<double, 4> v4;
	////q1 = v4;
	//v1 = v2 + v2;
	//q3 = q1;
	//q2 = q1 * 0.5;
	//q3 = q1 + q2;
}
void matrix_unit_test()
{
	//matrix<double, 3> m;
	//std::array<double, 3>& a = m[0];
	//double el = m[0][1];
	//m[0][2] = 1;
	//return;
}
void tensor_unit_test()
{
	auto M1 = tens::matrix<double, DIM>(tens::MATRIXINITTYPE::INDENT);
	auto M2 = tens::matrix<double, DIM>(tens::MATRIXINITTYPE::INDENT);
	M1 = std::move(M2);
	std::vector<std::vector<double>> V1(10);
	std::vector<std::vector<double>> V2(20);
	V1 = std::move(V2);

	//q1.vector_product(q2);
	auto mm = tens::matrix<double, DIM>(tens::MATRIXINITTYPE::INDENT);
	tens::shared_handler_basis<double, DIM> rr(mm);
	//typedef std::shared_ptr<> shared_handler;
	std::shared_ptr<double> xx = std::make_shared<double>();
	tens::matrix<double, DIM>* R = new tens::matrix<double, DIM>(tens::generate_rand_ort<double, DIM>());
	tens::shared_handler_basis<double, DIM> MM0(*R);
	for (size_t i = 0; i < 10000000; i++)
	{
		tens::matrix<double, DIM>* R = new tens::matrix<double, DIM>(tens::generate_rand_ort<double, DIM>());
		tens::matrix<double, DIM>* Q = new tens::matrix<double, DIM>(tens::generate_rand_ort<double, DIM>());
		tens::matrix<double, DIM>* P = new tens::matrix<double, DIM>(tens::generate_rand_ort<double, DIM>());
		auto sm1 = new tens::shared_handler_basis<double, DIM>(*R); 
		auto sm2 = new tens::shared_handler_basis<double, DIM>(*Q); 
		auto sm3 = new tens::shared_handler_basis<double, DIM>(*P);
		tens::Tensor<double, DIM>* t0 = new tens::Tensor<double, DIM>(*Q, *sm1); bool ort = sm1->check_ort_basis();
		tens::Tensor<double, DIM>* t1 = new tens::Tensor<double, DIM>(*P, *sm2);  ort = sm2->check_ort_basis();
		delete R;
		delete Q;
		tens::Tensor<double, DIM> t4(std::move(*t1)); delete t1;
		t1 = new tens::Tensor<double, DIM>(*P, *sm2);
		delete P;

		tens::Tensor<double, DIM> t3(*t0);
		*t0 =  t4;
		*t1 = t4 + *t0;
		t4 -= *t1;
		t4 *= *t0;
		t4 += *t0;
		t4.change_basis(*sm2);
		t1->change_basis(*sm3);
		delete sm1;
		delete sm2;
		delete sm3;
		delete t0;
		delete t1;
	}
}

void tens_vect_unit_test()
{
	tens::matrix<double, DIM>* R = new tens::matrix<double, DIM>(tens::generate_rand_ort<double, DIM>());
	auto sm1 = new  tens::shared_handler_basis<double, DIM>(*R);

	tens::array<double, DIM> v0(0.0);
	v0[0] = 1.0;
	v0[1] = 1.5;
	v0[2] = -0.5;
	tens::vector<double, 3> V0(v0, *sm1);
	tens::vector<double, 3> V1(v0, *sm1);
	tens::vector<double, 3> V3 = V0 + V1;
	tens::Tensor<double, DIM>* t0 = new tens::Tensor<double, DIM>(*R, *sm1);
	tens::Tensor<double, DIM>* t1 = new tens::Tensor<double, DIM>(*R, *sm1);
	tens::Tensor<double, DIM>  t2 = tens::outer_product(V0, V3); std::cout << t2;
	tens::Tensor<double, DIM>  t3 = tens::outer_product(V3, V0); std::cout << t3;
}

void vector_unit_test()
{
	tens::matrix<double, DIM>* R1 = new tens::matrix<double, DIM>(tens::generate_rand_ort<double, DIM>());
	tens::matrix<double, DIM>* R2 = new tens::matrix<double, DIM>(tens::generate_rand_ort<double, DIM>());
	auto sm1 = new  tens::shared_handler_basis<double, DIM>(tens::generate_rand_ort<double, DIM>());
	auto sm2 = new  tens::shared_handler_basis<double, DIM>(*R2);
	delete R1;
	delete R2;
	tens::array<double, DIM> v0(0.0);
	v0[0] = 1.0;
	v0[1] = 1.5;
	v0[2] = -0.5;
	tens::vector<double, 3> V0(v0, *sm1);
	tens::vector<double, 3> V1(v0, *sm2);
	{
		tens::vector<double, 3> V3 = V0 + V1;
		tens::vector<double, 3> V4 = V0 - V1;
		auto V5 = V0 * V1;
	}
	tens::vector<double, 3> V3 = std::move(V0);
	//V3 = std::move(V1);
	delete sm1;
	std::cout << V3 << " " << V1;
	//vector<double, DIM>  V4 = vector_product(V3, V1);
	std::cout << vector_product(V3, V1);
	//tens::array<double, 40> VV(0.0);
	//VV.normalize();
	return;
}
int main() 
{
	tens::is_small_value(1e-10);
	run_test();
	//expect(true, "this is good expext message");
	//expect(false, "this is bad expext message");
	return 0;
	//tens_vect_unit_test();
	//vector_unit_test();
	//tens::Tensor_unit_test();
	//matrix_unit_test();

	return 0;
}
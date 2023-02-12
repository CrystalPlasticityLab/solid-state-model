#pragma once
#include "../tensor/tensor.h"
#include "../tensor/vector.h"
#include "../tensor/quat.h"
#include <vector>

using namespace tens;

namespace state {
    template<typename T, size_t N>
    class state : private shared_handler_basis<T, N>
    {
		typedef std::reference_wrapper<shared_handler_basis<T,N>>  ref_wrap;
		typedef shared_handler_basis<T, N>  handler;
    private:
        std::vector<ref_wrap> object;
    public:
        Tensor<T, N> add_object(matrix<T,N>& m){
            auto t = Tensor<T, N>(m, static_cast<const handler&>(*this));
            object.push_back(t);
            return t;
        }
        vector<T, N> add_object(array<T,N>& a){
            auto v = vector<T, N>(a, static_cast<const handler&>(*this));
            object.push_back(v);
            return v;
        }
        state(const handler& basis) : handler(basis){

        };
        ~state(){};
    };
}
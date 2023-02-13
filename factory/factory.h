#pragma once
#include "../tensor/tensor.h"
#include "../tensor/vector.h"
#include "../tensor/quat.h"
#include "../tensor/error.h"
#include <vector>
#include <unordered_map>
#include <algorithm>

using namespace tens; 

namespace factory {
    std::hash<std::string> hasher;
    struct StringHasher {
        size_t operator()(const std::string& t) const {
            return hasher(t);
        }
    };

    template<typename T, size_t N>
    class factory : private shared_handler_basis<T, N>
    {
		typedef shared_handler_basis<T, N>  handler;
    private:
        factory() {};
    protected:
        tensor<T, N>* new_object(matrix<T, N>& m) {
            return new tensor<T, N>(m, static_cast<const handler&>(*this));
        }
        vector<T, N>* new_object(array<T, N>& a) {
            return new vector<T, N>(a, static_cast<const handler&>(*this));
        }
        factory(const handler& basis) : handler(basis) { };
        ~factory() {};
    public:
    };


    template<typename T, size_t N>
    class state : protected factory<T, N>, 
                  private std::unordered_map<const std::string, shared_handler_basis<T, N>*, StringHasher>
    {
        typedef shared_handler_basis<T, N>  handler;
    private:
        state() {};
        void push(const std::string& name, handler* item) {
            if (this->find(name) != this->end()) {
                throw ErrorAccess::Exists();
            };
            this->insert({ name, item });
        };
        state(const state& s) { }; // copy ctor
        state(state&& m)noexcept { };  // move ctor
        state& operator=(const state& s) {};  // copy assign
        state& operator= (state&& rhs) noexcept { }; // move assign
    public:
        state(const handler& basis) : factory<T,N>(basis) {};
        ~state() {
            // TODO : does not work -> fix
            for (auto it = this->begin(); it != this->end(); ++it) {
                //delete it->second;
                //std::cout << " " << it->first << ":" << *it->second;
            }
        };

        void push(const std::string& name, array<T, N>& a) {
            push(name, factory<T, N>::new_object(a));
        };

        void push(const std::string& name, matrix<T, N>& m) {
            push(name, factory<T, N>::new_object(m));
        };

        template<typename type>
        type& get(const std::string& name) {
            // TODO check type and item type
            // throw ErrorAccess::WrongCast(); 
            return static_cast<type&>(*(*this)[name]);
        }
    };
}


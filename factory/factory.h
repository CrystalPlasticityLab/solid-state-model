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
    class state : private shared_handler_basis<T, N>,
                  public std::unordered_map<const std::string, std::shared_ptr<shared_handler_basis<T, N>>, StringHasher>
    {
        typedef shared_handler_basis<T, N>  handler;
        typedef std::shared_ptr<shared_handler_basis<T, N>>  handler_ptr;
        typedef std::unordered_map<const std::string, std::shared_ptr<shared_handler_basis<T, N>>, StringHasher> map;
    private:
        enum class TYPEOBJECT {
            VECTOR,
            TENSOR,
            UNDEFINED
        };
        state() {}; 
        TYPEOBJECT get_type_name(const std::string& name) {
            if (name.find("tens::vector") != std::string::npos) { return TYPEOBJECT::VECTOR; }
            if (name.find("tens::tensor") != std::string::npos) { return TYPEOBJECT::TENSOR; };
            return TYPEOBJECT::UNDEFINED;
        }

        TYPEOBJECT get_type_name(const handler_ptr& obj) {
            return get_type_name(typeid(*obj).name());
        }
        void push(const std::string& name, handler_ptr&& item) {
            if (this->find(name) != this->end()) {
                throw ErrorAccess::Exists();
            };
            this->insert({ name, std::move(item) });
        };

        std::shared_ptr<tensor<T, N>> new_object(const matrix<T, N>& m) {
            return std::make_shared<tensor<T, N>>(m, static_cast<const handler&>(*this));
        }

        std::shared_ptr<vector<T, N>> new_object(const array<T, N>& a) {
            return std::make_shared<vector<T, N>>(a, static_cast<const handler&>(*this));
        }

        void _deep_copy(const map& s) {
            auto& t = *this;
            const auto& basis = static_cast<const handler&>(*this);
            for (auto const& [key, val] : s) {
                TYPEOBJECT type = get_type_name(val);
                if (type == TYPEOBJECT::VECTOR){
                    t[key] = new_object(static_cast<const vector<T, N>&>(*val).get_comp_at_basis(basis));
                }
                if (type == TYPEOBJECT::TENSOR) {
                    t[key] = new_object(static_cast<const tensor<T, N>&>(*val).get_comp_at_basis(basis));
                }
            }
        }

    public:
        state(const handler& basis) {
            handler::_deep_copy(basis);
        };
        state(const handler& basis, const state& s) {
            handler::_deep_copy(basis);
            state::_deep_copy(s);
        };

        state(const state& s) : handler(s), map(s) { // copy ctor
            handler::_deep_copy(s);
            state::_deep_copy(s);
        }; 
        state& operator= (const state& s) {  // copy assign
            handler::_deep_copy(s);
            state::_deep_copy(s);
            return *this;
        };

        state(state&& s)noexcept : handler(std::move(s)), map(std::move(s)) { };  // move ctor

        state& operator= (state&& rhs) noexcept {  // move assign
            handler::operator = (static_cast<handler&&>(rhs));
            map::operator = (static_cast<map&&>(rhs));
            return *this;
        };
        ~state() {};

        void push(const std::string& name, array<T, N>& a) {
            push(name, new_object(a));
        };

        void push(const std::string& name, matrix<T, N>& m) {
            push(name, new_object(m));
        };

        template<typename type>
        type& get(const std::string& name) {
            TYPEOBJECT t = get_type_name((*this)[name]);
            if (typeid(*(*this)[name]).name() != typeid(type).name()) {
                const std::string msg = std::string(typeid((*this)[name]).name()) + " to " + std::string(typeid(type).name());
                throw ErrorAccess::WrongCast(msg);
            }
            return static_cast<type&>(*(*this)[name]);
        }
    };
}


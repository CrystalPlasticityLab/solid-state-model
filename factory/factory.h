#pragma once
#include <vector>
#include <unordered_map>
#include <algorithm>
#include "../tensor/tensor.h"
#include "../tensor/vector.h"
#include "../tensor/quat.h"
#include "../tensor/error.h"

using namespace tens; 

namespace factory {
    struct StringHasher {
        std::hash<std::string> hasher;
        size_t operator()(const std::string& t) const {
            return hasher(t);
        }
    };

    template<typename T, size_t N>
    class state : private basis<T, N>,
                  public std::unordered_map<const std::string, std::shared_ptr<basis<T, N>>, StringHasher>
    {
        typedef basis<T, N>  handler;
        typedef std::shared_ptr<basis<T, N>>  handler_ptr;
        typedef std::unordered_map<const std::string, std::shared_ptr<basis<T, N>>, StringHasher> map;
    private:
        state() {}; 
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
                const size_t rank_object = val->get_rank();
                switch (rank_object)
                {
                case 1:
                    t[key] = new_object(static_cast<const vector<T, N>&>(*val).get_comp_at_basis(basis));
                    break;
                case 2:
                    t[key] = new_object(static_cast<const tensor<T, N>&>(*val).get_comp_at_basis(basis));
                    break;
                default:
                    break;
                }
            }
        }

    public:
        size_t get_rank() const override { return 0; };
        state(const matrix<T,N>& m) : handler(m) {
            //handler::_deep_copy(basis);
        };

        state(const matrix<T, N>& m, const state& s) : handler(m) {
            //handler::_deep_copy(basis);
            state::_deep_copy(s);
        };

        //state(const handler& basis) {
        //    handler::_deep_copy(basis);
        //};
        //
        //state(const handler& basis, const state& s) {
        //    handler::_deep_copy(basis);
        //    state::_deep_copy(s);
        //};

        state(const state& s) : handler(static_cast<const basis<T, N>&>(s)), map(s) { // copy ctor
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

        void push(const std::string& name, const container<T, N>* c) {
            switch (c->get_rank())
            {
            case 1:
                push(name, new_object(static_cast<const array<T, N>&>(*c)));
                break;
            case 2:
                push(name, new_object(static_cast<const matrix<T, N>&>(*c)));
                break;
            default:
                break;
            }
        };

        template<typename type>
        type& get(const std::string& name) {
            if (typeid(*(*this)[name]).name() != typeid(type).name()) {
                const std::string msg = std::string(typeid((*this)[name]).name()) + " to " + std::string(typeid(type).name());
                throw ErrorAccess::WrongCast(msg);
            }
            return static_cast<type&>(*(*this)[name]);
        }
    };
}


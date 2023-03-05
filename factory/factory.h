#pragma once
#include <Vector>
#include <unordered_map>
#include <algorithm>
#include "../tensor/basis.h"
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

    template <typename T, size_t N>
    using object = std::shared_ptr<basis_base<T, N>>;


    template<typename T, size_t N>
    class state : public std::unordered_map<const std::string, object<T,N>, StringHasher>
    {
        typedef std::unordered_map<const std::string, Basis<T, N>, StringHasher> map;

    private:
        Basis<T, N> _basis;
    public:
        state(const Basis<T, N>& basis) : _basis(basis) {};

        void remove(const std::string& name) {
            this->erase(name);
        }

        template<size_t R>
        basis<T, N, R>& push(const std::string& name, FILL_TYPE fill_type = FILL_TYPE::ZERO) {
            if (this->find(name) != this->end()) {
                throw ErrorAccess::Exists();
            };
            container<T, N, R> c(fill_type);
            auto obj = std::make_shared<basis<T, N, R>>(basis<T, N, R>(std::move(c), _basis));
            this->insert({ name, obj });
            return *obj.get();// get<R>(name);
        }

        template<size_t R>
        basis<T, N, R>& get(const std::string& s) {
            const auto obj_ptr = (*this)[s];
            if (obj_ptr) {
                auto& obj = static_cast<basis<T, N, R>&>(*obj_ptr);
                if (R != obj.get_rank()) {
                    throw ErrorAccess::WrongTemplateType();
                }
                return obj;
            } else {
                throw ErrorAccess::NotExists();
            }
            throw UnHandled();
        };
    };
}


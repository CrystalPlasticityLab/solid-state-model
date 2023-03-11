#pragma once
#include <Vector>
#include <unordered_map>
#include <algorithm>
#include "../tensor/object.h"
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

    template <typename T>
    using object = Basis<T>;// std::shared_ptr<basi<T>>;


    template<typename T>
    class state : public std::unordered_map<const std::string, object<T>, StringHasher>
    {
        typedef std::unordered_map<const std::string, Basis<T>, StringHasher> map;

    private:
        Basis<T> _basis;
    public:
        state(const Basis<T>& object) : _basis(object) {};

        void remove(const std::string& name) {
            this->erase(name);
        }

        object<T>& push(const std::string& name, size_t N, size_t R, FILL_TYPE fill_type = FILL_TYPE::ZERO) {
            if (this->find(name) != this->end()) {
                throw ErrorAccess::Exists();
            };
            container<T> c(fill_type);
            auto obj = std::make_shared<object<T>>(object<T>(std::move(c), _basis));
            this->insert({ name, obj });
            return *obj.get();// get<R>(name);
        }

        object<T>& get(const std::string& s) {
            const auto obj_ptr = (*this)[s];
            if (obj_ptr) {
                auto& obj = static_cast<object<T>&>(*obj_ptr);
                //if (R != obj.rank()) {
                //    throw ErrorAccess::WrongTemplateType();
                //}
                return obj;
            } else {
                throw ErrorAccess::NotExists();
            }
            throw UnHandled();
        };
    };
}


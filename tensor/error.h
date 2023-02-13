#pragma once
#include <exception>
namespace ErrorMessage {
    const std::string DIV_BY_ZERO_MSG = std::string("Division by zero");
    const std::string NON_ORT_MATRIX = std::string("Matrix is not orthogonal");
    const std::string HAS_NO_ACCESS = std::string("Access is not permitted: called object is not an owner");
    const std::string UNORDERED_MAP_ITEM_EXIST = std::string("Item with this key alredy exists");
    const std::string WRONG_TYPE_CAST = std::string("Unable to cast");
}


namespace ErrorMath {
    class DivisionByZero: public std::exception {
        public:
        virtual const char* what() const noexcept{
            return ErrorMessage::DIV_BY_ZERO_MSG.c_str();
        };
    };

    class NonOrthogonal: public std::exception {
        public:
        virtual const char* what() const noexcept{
            return ErrorMessage::NON_ORT_MATRIX.c_str();
        };
    };
};

namespace ErrorAccess{
    class NoAccess: public std::exception {
        public:
        virtual const char* what() const noexcept{
            return ErrorMessage::HAS_NO_ACCESS.c_str();
        };
    };

    class Exists : public std::exception {
    public:
        virtual const char* what() const noexcept {
            return ErrorMessage::UNORDERED_MAP_ITEM_EXIST.c_str();
        };
    };

    class WrongCast : public std::exception {
    public:
        virtual const char* what() const noexcept {
            return ErrorMessage::WRONG_TYPE_CAST.c_str();
        };
    };
}
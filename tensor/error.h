#pragma once
#include <exception>
namespace ErrorMessage {
    const std::string UNHANDLED_BEHAVIOR = std::string("Non handled behavior, possible it is just not implemented yet ");
    const std::string DIV_BY_ZERO_MSG = std::string("Division by zero");
    const std::string NON_ORT_MATRIX = std::string("Matrix is not orthogonal");
    const std::string HAS_NO_ACCESS = std::string("Access is not permitted: called object is not an owner");
    const std::string UNORDERED_MAP_ITEM_EXIST = std::string("Item with this key alredy exists");
    const std::string UNORDERED_MAP_ITEM_NOT_EXIST = std::string("Item with this key is not exists");
    const std::string WRONG_TEMPLATE_CAST = std::string("Wrong template value");
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

class UnHandled : public std::exception {
public:
    virtual const char* what() const noexcept {
        return ErrorMessage::UNHANDLED_BEHAVIOR.c_str();
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

    class NotExists : public std::exception {
    public:
        virtual const char* what() const noexcept {
            return ErrorMessage::UNORDERED_MAP_ITEM_NOT_EXIST.c_str();
        };
    };

    class WrongTemplateType : public std::exception {
    public:
        virtual const char* what() const noexcept {
            return ErrorMessage::WRONG_TEMPLATE_CAST.c_str();
        };
    };
}
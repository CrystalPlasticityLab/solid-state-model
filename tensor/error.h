#pragma once
#include <exception>
namespace ErrorMessage {
    const std::string UNHANDLED_BEHAVIOR = std::string("This feature is not implemented yet");
    const std::string DIV_BY_ZERO_MSG = std::string("Division by zero");
    const std::string NON_ORT_MATRIX = std::string("Matrix is not orthogonal");
    const std::string NO_CAST_TO_SCALAR = std::string("Object has no cast to scalar");
    const std::string UNORDERED_MAP_ITEM_EXIST = std::string("Item with this key alredy exists");
    const std::string UNORDERED_MAP_ITEM_NOT_EXIST = std::string("Item with this key is not exists");
    const std::string WRONG_TEMPLATE_CAST = std::string("Wrong template value");
    const std::string SHAPE_MISMATCH = std::string("Containers have dofferent shapes");
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

    class ShapeMismatch : public std::exception {
    public:
        virtual const char* what() const noexcept {
            return ErrorMessage::SHAPE_MISMATCH.c_str();
        };
    };
};

class NoImplemetationYet : public std::exception {
public:
    virtual const char* what() const noexcept {
        return ErrorMessage::UNHANDLED_BEHAVIOR.c_str();
    };
};

namespace ErrorAccess{
    class NoCastScalar: public std::exception {
        public:
        virtual const char* what() const noexcept{
            return ErrorMessage::NO_CAST_TO_SCALAR.c_str();
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
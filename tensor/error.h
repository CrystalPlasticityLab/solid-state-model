
#include <exception>
namespace ErrorMessage {
    const std::string DIV_BY_ZERO_MSG = std::string("Division by zero");
    const std::string NON_ORT_MATRIX = std::string("Matrix is not orthogonal");
    const std::string HAS_NO_ACCESS = std::string("access is not permitted: called object is not an owner");
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
}
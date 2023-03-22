#include <boost/operators.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <Eigen/Dense>

// Quotient
template <class NT_>
class Toto
    : boost::additive2 < Toto<NT_>, NT_ >
{
public:
    NT_ num, den;

    Toto() {}
    Toto(const NT_& a) {}

    Toto(const NT_& a, const NT_& b) {}

    template <class T>
    Toto operator+=(const T&) { return *this; }
};

// Sqrt_extension:
template <class NT_>
class Toto2
    : boost::additive2 < Toto2<NT_>, NT_ >
{
public:
    NT_ num, den;
    Toto2() {}
    Toto2(const NT_& nt)
    {}

    template <class T>
    Toto2 operator+=(const T&) { return *this; }
};




void toto(Toto<boost::multiprecision::cpp_int>& x)
{
    x.den;
}

// No idea if needed
namespace Eigen {
    template<class> struct NumTraits;
    template<class NT> struct NumTraits<Toto<NT> >
    {
        typedef Toto<NT> Real;
        typedef Toto<NT> NonInteger;
        typedef Toto<NT> Nested;
        typedef Toto<NT> Literal;

        static inline Real epsilon() { return NumTraits<NT>::epsilon(); }
        static inline Real dummy_precision() { return NumTraits<NT>::dummy_precision(); }

        enum {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = NumTraits<NT>::RequireInitialization,
            ReadCost = 2 * NumTraits<NT>::ReadCost,
            AddCost = 150,
            MulCost = 100
        };
    };

    template<class NT> struct NumTraits<Toto2<NT> >
    {
        typedef Toto2<NT> Real;
        typedef Toto2<NT> NonInteger;
        typedef Toto2<NT> Nested;
        typedef Toto2<NT> Literal;

        static inline Real epsilon() { return NumTraits<NT>::epsilon(); }
        static inline Real dummy_precision() { return NumTraits<NT>::dummy_precision(); }

        enum {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = NumTraits<NT>::RequireInitialization,
            ReadCost = 2 * NumTraits<NT>::ReadCost,
            AddCost = 150,
            MulCost = 100
        };
    };
}


int main() {
    typedef Toto2<int> NT;
    Eigen::Matrix<NT, 3, 3> m(3, 3);

    m + m;

    return 0;
}

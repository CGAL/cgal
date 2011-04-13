//#line 213 "testtypes.fw"

#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#include <CGAL/basic.h>
#include <iostream>
//#include <iostream.h>
#include <CGAL/double.h>

#ifndef PEDANTIC
#include <iostream>
#endif

//#line 158 "testtypes.fw"

CGAL_BEGIN_NAMESPACE

class TestfieldC {
public:
    TestfieldC() {}
    TestfieldC(int) ;
    TestfieldC(unsigned char, signed char, double d) {_d = d;}
    TestfieldC operator-() const {return TestfieldC(-_d);}
    TestfieldC operator+(TestfieldC tf) const
                    { return TestfieldC(_d + tf._d);}
    TestfieldC operator-(TestfieldC tf) const
                    { return TestfieldC(_d - tf._d);}
    TestfieldC operator*(TestfieldC tf) const
                    { return TestfieldC(_d * tf._d);}
    TestfieldC operator/(TestfieldC tf) const
                    { return TestfieldC(_d / tf._d);}
    bool operator==(TestfieldC tf) const {return _d == tf._d;}
    bool operator!=(TestfieldC tf) const {return _d != tf._d;}
    bool operator<(TestfieldC tf) const {return _d < tf._d;}
    bool operator>(TestfieldC tf) const {return _d > tf._d;}
    bool operator<=(TestfieldC tf) const {return _d <= tf._d;}
    bool operator>=(TestfieldC tf) const {return _d >= tf._d;}
friend inline double to_double(TestfieldC tf);
protected:
    TestfieldC(float) ;    // not implemented.
    TestfieldC(double d) { _d = d;}
    TestfieldC(int, int, double) ;    // not implemented.
    double _d;
};

inline TestfieldC::TestfieldC(int d)
{
    if (d < 0 || d > 127)
        std::cerr << "Non-standard use of number constructor: RT("<<d<<").\n";
    _d = double(d);
}

inline double to_double(TestfieldC tf)
{
    return tf._d;
}

inline bool is_finite(TestfieldC tf)
{
    return is_finite(to_double(tf));
}

inline bool is_valid(TestfieldC tf)
{
    return is_valid(to_double(tf));
}
//#line 222 "testtypes.fw"

//#line 15 "testtypes.fw"


class TestrepH {
public:
    TestrepH() {}
    TestrepH(int d) ;   // only for testrep(0) and testrep(1).
    TestrepH(unsigned char, signed char, double d) {_d = d;}
    TestrepH operator-() const {return TestrepH(-_d);}
    TestrepH operator+(TestrepH tr) const
                    {return TestrepH(_d + tr._d);}
    TestrepH operator-(TestrepH tr) const
                    {return TestrepH(_d - tr._d);}
    TestrepH operator*(TestrepH tr) const
                    {return TestrepH(_d * tr._d);}
//    double to_double() const {return _d;}
    bool operator==(TestrepH tr) const {return _d == tr._d;}
    bool operator!=(TestrepH tr) const {return _d != tr._d;}
    bool operator<(TestrepH tr) const {return _d < tr._d;}
    bool operator>(TestrepH tr) const {return _d > tr._d;}
    bool operator<=(TestrepH tr) const {return _d <= tr._d;}
    bool operator>=(TestrepH tr) const {return _d >= tr._d;}
friend inline double to_double(TestrepH tf);
#ifndef PEDANTIC
// operators added to make life for Stefan easier...
    TestrepH operator/(TestrepH tr) const
                    {return TestrepH(_d / tr._d);}
    void operator*=(TestrepH tr)
                    {_d *= tr._d;}
#endif
protected:
    TestrepH(float ) ;  // not implemented
    TestrepH(double d) {_d = d;}
    TestrepH(int, int, double) ;    // not implemented
    double _d;
};

inline TestrepH::TestrepH(int d)
{
    if (d < 0 || d > 127)
        std::cerr << "Non-standard use of number constructor: RT("<<d<<").\n";
    _d = double(d);
}

inline bool is_finite(TestrepH tr)
{
    return is_finite(to_double(tr));
}

inline bool is_valid(TestrepH tr)
{
    return is_valid(to_double(tr));
}

inline double to_double(TestrepH tr)
{
    return tr._d;
}

//#line 223 "testtypes.fw"

//#line 103 "testtypes.fw"

class TestfieldH {
public:
    TestfieldH() {}
    TestfieldH(int) ;
    TestfieldH(unsigned char, signed char, double d) {_d = d;}
    TestfieldH(TestrepH tr) {_d = to_double(tr);}
    TestfieldH operator-() const {return TestfieldH(-_d);}
    TestfieldH operator+(TestfieldH tf) const
                    { return TestfieldH(_d + tf._d);}
    TestfieldH operator-(TestfieldH tf) const
                    { return TestfieldH(_d - tf._d);}
    TestfieldH operator*(TestfieldH tf) const
                    { return TestfieldH(_d * tf._d);}
    TestfieldH operator/(TestfieldH tf) const
                    { return TestfieldH(_d / tf._d);}
    bool operator==(TestfieldH tf) const {return _d == tf._d;}
    bool operator!=(TestfieldH tf) const {return _d != tf._d;}
    bool operator<(TestfieldH tf) const {return _d < tf._d;}
    bool operator>(TestfieldH tf) const {return _d > tf._d;}
    bool operator<=(TestfieldH tf) const {return _d <= tf._d;}
    bool operator>=(TestfieldH tf) const {return _d >= tf._d;}
friend inline double to_double(TestfieldH tf);
protected:
    TestfieldH(float) ;    // not implemented
    TestfieldH(double d) { _d = d;}
    TestfieldH(int, int, double) ;    // not implemented
    double _d;
};

inline TestfieldH::TestfieldH(int d)
{
    if (d < 0 || d > 127)
        std::cerr << "Non-standard use of number constructor: RT("<<d<<").\n";
    _d = double(d);
}

inline double to_double(TestfieldH tf)
{
    return tf._d;
}

inline bool is_finite(TestfieldH tf)
{
    return is_finite(to_double(tf));
}

inline bool is_valid(TestfieldH tf)
{
    return is_valid(to_double(tf));
}
//#line 224 "testtypes.fw"

#ifndef PEDANTIC
std::istream & operator>>(std::istream &is, TestrepH &tr)
{
    double d;
    unsigned char ucdum = 0;
    signed char scdum = 0;
    is >> d;
    tr = TestrepH(ucdum, scdum, d);
    return is;
}

std::ostream & operator<<(std::ostream &os, TestrepH tr)
{
    os << to_double(tr);
    return os;
}
#endif

CGAL_END_NAMESPACE

#endif












#ifndef CGAL_BOX_INTERSECTION_D_BOX_LIMITS_H
#define CGAL_BOX_INTERSECTION_D_BOX_LIMITS_H

#include <CGAL/basic.h>
#include <cfloat>
#include <climits>

CGAL_BEGIN_NAMESPACE

template<class T>
struct box_limits {};

template<>
struct box_limits<int> {
    static int inf() { return INT_MIN; }
    static int sup() { return INT_MAX; }
};

template<>
struct box_limits<unsigned int> {
    static unsigned int inf() { return 0; }
    static unsigned int sup() { return UINT_MAX; }
};

template<>
struct box_limits<float> {
    struct X {
        X() : u( 0x7f800000 ) {}
        union {
            float f;
            unsigned int u;
        };
    };
    static X x;
    static float inf() { return -sup(); }
    static float sup() { return x.f; }
    //static double sup() { return FLT_MAX; }
};

template<>
box_limits<float>::X box_limits<float>::x;

template<>
struct box_limits<double> {
    struct X {
        X() : u( 0x7FF0000000000000ull ) {}
        union {
            double d;
            unsigned long long u;
        };
    };
    static X x;
    static double inf() { return -sup(); }
    static double sup() { return x.d; }
    //static double sup() { return DBL_MAX; }
};

template<>
box_limits<double>::X box_limits<double>::x;


CGAL_END_NAMESPACE

#endif

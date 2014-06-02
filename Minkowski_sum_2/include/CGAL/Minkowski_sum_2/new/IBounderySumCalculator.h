#ifndef IBOUNDERY_SUM_CALCULATOR_HEADER
#define IBOUNDERY_SUM_CALCULATOR_HEADER

#include <CGAL/basic.h>
#include <CGAL/Polygon_2.h>

template <class Kernel_, class Container_> class IBounderySumCalculator {

protected:
    typedef Kernel_ Kernel;
    typedef CGAL::Polygon_2<Kernel, Container_> Polygon_2;
public:
    virtual void calc_sum(Polygon_2 &a, Polygon_2 &b, Polygon_2 &res_poly) = 0;

};

#endif

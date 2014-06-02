#ifndef BOUNDERYCONCALC_HEADER
#define BOUNDERYCONCALC_HEADER
#include "IBounderySumCalculator.h"
#include "MyMink.h"
#include <CGAL/Env_default_diagram_1.h>
#include <CGAL/envelope_2.h>

template <class Kernel_, class Container_> class BounderyConvCalc: public IBounderySumCalculator< Kernel_,  Container_> {

public:

    typedef Arr_segment_traits_2<Kernel>                    Traits_2;//Segment_traits_2;
    typedef typename Traits_2::X_monotone_curve_2   Segment_2;
    typedef std::list<Segment_2>                    Segments_list;

    virtual void calc_sum(Polygon_2 &a, Polygon_2 &b, Polygon_2 &res_poly) {
        Segments_list reduced_conv;
        _mink->buildReducedConvolution(a, b, reduced_conv);

    }

private:

    Minkowski_sum_by_convolution_lien_2<Kernel, Container_> *_mink;

};
#endif
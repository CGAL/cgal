// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : point_generators_3.h
// chapter       : $CGAL_Chapter: Geometric Object Generators $
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// source        : generators.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// 3D Point Generators
// ============================================================================

#ifndef CGAL_POINT_GENERATORS_3_H
#define CGAL_POINT_GENERATORS_3_H 1
#ifndef CGAL_GENERATORS_H
#include <CGAL/generators.h>
#endif
#ifndef CGAL_POINT_GENERATORS_2_H
#include <CGAL/point_generators_2.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class P, class Creator = Creator_uniform_3<double,P> >
#else
template < class P, class Creator >
#endif
class Random_points_in_sphere_3 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef Random_points_in_sphere_3<P,Creator> This;
    Random_points_in_sphere_3( double r = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed in the open sphere with radius r, i.e. |`*g'| < r .
        // Three random numbers are needed from `rnd' for each point.
    : Random_generator_base<P>( r, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_in_sphere_3<P,Creator>::
generate_point() {
   do {
       Creator creator;
       d_item = creator( d_range * ( 2 * _rnd.get_double() - 1.0),
                         d_range * ( 2 * _rnd.get_double() - 1.0),
                         d_range * ( 2 * _rnd.get_double() - 1.0));
   } 
   while (CGAL::to_double(d_item.x() * d_item.x() + d_item.y() * d_item.y()) >=
          d_range * d_range);
}


#ifndef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class P, class Creator = Creator_uniform_3<double,P> >
#else
template < class P, class Creator >
#endif
class Random_points_on_sphere_3 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef Random_points_on_sphere_3<P,Creator> This;
    Random_points_on_sphere_3( double r = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the circle with radius r, i.e. |`*g'| == r . A
        // single random number is needed from `rnd' for each point.
    : Random_generator_base<P>( r, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_on_sphere_3<P,Creator>::
generate_point() {
    double alpha = _rnd.get_double() * 2.0 * M_PI;
    double z     = 2 * _rnd.get_double() - 1.0;
    double r     = std::sqrt( 1 - z * z);
    Creator creator;
    d_item = creator( d_range * r * CGAL_CLIB_STD::cos(alpha),
                      d_range * r * CGAL_CLIB_STD::sin(alpha),
                      d_range * z);
}


#ifndef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class P, class Creator = Creator_uniform_3<double,P> >
#else
template < class P, class Creator >
#endif
class Random_points_in_cube_3 : public Random_generator_base<P>{
    void generate_point();
public:
    typedef Random_points_in_cube_3<P,Creator> This;
    Random_points_in_cube_3( double a = 1, Random& rnd = default_random)
    : Random_generator_base<P>( a, rnd) { generate_point(); }
    This& operator++()    {
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P, class Creator >
void
Random_points_in_cube_3<P,Creator>::
generate_point() {
    Creator creator;
    d_item = creator( d_range * ( 2 * _rnd.get_double() - 1.0),
                      d_range * ( 2 * _rnd.get_double() - 1.0),
                      d_range * ( 2 * _rnd.get_double() - 1.0));
}

CGAL_END_NAMESPACE
#endif // CGAL_POINT_GENERATORS_3_H //
// EOF //

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
// file          : point_generators_2.h
// chapter       : $CGAL_Chapter: Geometric Object Generators $
// package       : $CGAL_Package: Generator 2.12 (28 Jul 1999) $
// source        : generators.fw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Lutz Kettner  <kettner@inf.ethz.ch>
//
// coordinator   : INRIA, Sophia Antipolis
//
// 2D Point Generators
// ============================================================================

#ifndef CGAL_POINT_GENERATORS_2_H
#define CGAL_POINT_GENERATORS_2_H 1
#ifndef CGAL_GENERATORS_H
#include <CGAL/generators.h>
#endif
#include <CGAL/ch_value_type.h>

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class P, class Creator = Creator_uniform_2<double,P> >
#else
template < class P, class Creator >
#endif
class Random_points_in_disc_2 : public Random_generator_base<P>{
    void generate_point();
public:
    typedef Random_points_in_disc_2<P,Creator> This;
    Random_points_in_disc_2( double r = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed in the open disc with radius r, i.e. |`*g'| < r .
        // Two random numbers are needed from `rnd' for each point.
    : Random_generator_base<P>(r, rnd) { generate_point(); }
    This& operator++() {
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
Random_points_in_disc_2<P,Creator>::
generate_point() {
    double alpha = _rnd.get_double() * 2.0 * M_PI;
    double r     = d_range * CGAL_CLIB_STD::sqrt( _rnd.get_double());
    Creator creator;
    d_item = creator( r * CGAL_CLIB_STD::cos(alpha), 
                      r * CGAL_CLIB_STD::sin(alpha));
}


#ifndef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class P, class Creator = Creator_uniform_2<double,P> >
#else
template < class P, class Creator >
#endif
class Random_points_on_circle_2 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef Random_points_on_circle_2<P,Creator> This;
    Random_points_on_circle_2( double r = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the circle with radius r, i.e. |`*g'| == r . A
        // single random number is needed from `rnd' for each point.
    : Random_generator_base<P>(r, rnd) { generate_point(); }
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
Random_points_on_circle_2<P,Creator>::
generate_point() {
    double a = _rnd.get_double() * 2.0 * M_PI;
    Creator creator;
    d_item = creator( d_range * CGAL_CLIB_STD::cos(a), 
                      d_range * CGAL_CLIB_STD::sin(a));
}


#ifndef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class P, class Creator = Creator_uniform_2<double,P> >
#else
template < class P, class Creator >
#endif
class Random_points_in_square_2 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef Random_points_in_square_2<P,Creator> This;
    Random_points_in_square_2( double a = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed in the half-open square with side length a,
        // centered around the origin, i.e. \forall p = `*g': -\frac{a}{2}
        // <= p.x() < \frac{a}{2} and -\frac{a}{2} <= p.y() < \frac{a}{2}
        // . Two random numbers are needed from `rnd' for each point.
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
Random_points_in_square_2<P,Creator>::
generate_point() {
    Creator creator;
    d_item = creator( d_range * (2 * _rnd.get_double() - 1.0),
                      d_range * (2 * _rnd.get_double() - 1.0));
}


#ifndef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class P, class Creator = Creator_uniform_2<double,P> >
#else
template < class P, class Creator >
#endif
class Random_points_on_square_2 : public Random_generator_base<P> {
    void generate_point();
public:
    typedef Random_points_on_square_2<P,Creator> This;
    Random_points_on_square_2( double a = 1, Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the boundary of the square with side length a,
        // centered around the origin, i.e. \forall p = `*g': one
        // coordinate is either \frac{a}{2} or -\frac{a}{2} and for the
        // other coordinate c holds -\frac{a}{2} <= c < \frac{a}{2} . A
        // single random number is needed from `rnd' for each point.
    : Random_generator_base<P>( a, rnd)  { generate_point(); }
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
Random_points_on_square_2<P,Creator>::
generate_point() {
    double d = _rnd.get_double() * 4.0;
    int    k = int(d);
    d = d_range * (2 * (d - k) - 1.0);
    CGAL_assertion( - d_range <= d && d < d_range);
    Creator creator;
    switch (k) {
    case 0:
        d_item = creator(        d, -d_range);
        break;
    case 1:
        d_item = creator(        d,  d_range);
        break;
    case 2:
        d_item = creator( -d_range,        d);
        break;
    case 3:
        d_item = creator(  d_range,        d);
        break;
    }
}


#ifndef CGAL_CFG_NO_DEFAULT_PREVIOUS_TEMPLATE_ARGUMENTS
template < class P, class Creator = Creator_uniform_2<double,P> >
#else
template < class P, class Creator >
#endif
class Random_points_on_segment_2 : public Random_generator_base<P> {
    P _p;
    P _q;
    void generate_point();
public:
    typedef Random_points_on_segment_2<P,Creator> This;
    Random_points_on_segment_2( const P& p = P( -1, 0),
                                const P& q = P(  1, 0),
                                Random& rnd = default_random)
        // g is an input iterator creating points of type `P' uniformly
        // distributed on the segment from p to q except q, i.e. `*g' ==
        // \lambda p + (1-\lambda)\, q where 0 <= \lambda < 1 . A single
        // random number is needed from `rnd' for each point.
    : Random_generator_base<P>( max( max( to_double(p.x()), to_double(q.x())),
                                          max( to_double(p.y()),
                                               to_double(q.y()))),
                                 rnd) , _p(p), _q(q) {
        generate_point();
    }
    const P&  source() const { return _p; }
    const P&  target() const { return _q; }
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
Random_points_on_segment_2<P,Creator>::
generate_point() {
    double la = _rnd.get_double();
    double mu = 1.0 - la;
    Creator creator;
    d_item = creator( mu * to_double(_p.x()) + la * to_double(_q.x()),
                      mu * to_double(_p.y()) + la * to_double(_q.y()));
}

template < class P >
class Points_on_segment_2 : public Generator_base<P> {
    P _p;
    P _q;
    std::size_t  d_i;
    std::size_t  d_mx;
    void generate_point();
public:
    typedef Points_on_segment_2<P> This;
    Points_on_segment_2() {}
    Points_on_segment_2( const P& p, const P& q,
                         std::size_t mx, std::size_t i = 0)
    : Generator_base<P>( max( max( to_double(p.x()), to_double(q.x())),
                                   max( to_double(p.y()), to_double(q.y())))),
                          _p(p), _q(q), d_i(i), d_mx(mx) {
        generate_point();
    }
    const P&  source() const { return _p; }
    const P&  target() const { return _q; }
    // Sufficient equality test.
    bool operator==( const This& base) const { return ( d_i == base.d_i); }
    bool operator!=( const This& base) const { return ! operator==(base); }
    This& operator++()    {
        d_i++;
        generate_point();
        return *this;
    }
    This  operator++(int) {
        This tmp = *this;
        ++(*this);
        return tmp;
    }
};

template < class P >
void
Points_on_segment_2<P>::
generate_point() { d_item = _p + (_q-_p) * d_i / (d_mx-1); }

template <class OutputIterator, class Creator>
OutputIterator
points_on_square_grid_2( double a, std::size_t n, OutputIterator o,
                         Creator creator)
{
    if  (n == 0)
        return o;
    int m = int(CGAL_CLIB_STD::ceil(std::sqrt(static_cast<double>(n))));
    double base = -a;  // Left and bottom boundary.
    double step = (2*a)/(m - 1);
    int j = 0;
    double px = base;
    double py = base;
    *o++ = creator( px, py);
    for (std::size_t i = 1; i < n; i++) {
        j++;
        if ( j == m) {
            px = base;
            py = py + step;
            j = 0;
        } else {
            px = px + step;
        }
        *o++ = creator( px, py);
    }
    return o;
}

template <class OutputIterator, class P>
OutputIterator
_points_on_square_grid_2( double a, std::size_t n, OutputIterator o,
                          const P*)
{
    return points_on_square_grid_2(a, n, o, Creator_uniform_2<double,P>());
}

// Works only if value_type is defined for the OutputIterator.
template <class OutputIterator>
OutputIterator
points_on_square_grid_2( double a, std::size_t n, OutputIterator o)
{
    return _points_on_square_grid_2( a, n, o, ch_value_type(o));
}


template <class P, class OutputIterator>
OutputIterator
points_on_segment_2( const P& p, const P& q, std::size_t n,
                     OutputIterator o)
    // creates n points regular spaced on the segment from p to q, i.e.
    // \forall i: 0 <= i < n: o[i] := \frac{n-i-1}{n-1} p + \frac{i}{n-1
    // } q.
{
    for (std::size_t i = 0; i < n; i++) {
        *o++ = p + (q-p) * i / (n-1);
    }
    return o;
}
template <class ForwardIterator, class Creator>
void perturb_points_2( ForwardIterator first,
                       ForwardIterator last,
                       double xeps,
                       double yeps,
                       Random& rnd,
                       Creator creator)
    // perturbs the points in the range [`first',`last') by replacing
    // each point with a random point from the rectangle `xeps' \times
    // `yeps' centered around the original point. Two random numbers are
    // needed from `rnd' for each point. Precondition:
    // The expression `to_double((*first).x())' and `to_double((
    // *begin).y())' must be legal.
{
    xeps *= 2.0;
    yeps *= 2.0;
    for ( ; first != last; ++first) {
        double x = to_double( (*first).x());
        double y = to_double( (*first).y());
        x += xeps * (rnd.get_double() - 0.5);
        y += yeps * (rnd.get_double() - 0.5);
        *first = creator( x, y);
    }
}

template <class ForwardIterator, class P>
void _perturb_points_2( ForwardIterator first,
                        ForwardIterator last,
                        double xeps,
                        double yeps,
                        Random& rnd,
                        const P*)
{
    perturb_points_2( first, last, xeps, yeps, rnd,
                      Creator_uniform_2<double,P>());
}

template <class ForwardIterator>
void perturb_points_2( ForwardIterator first,
                       ForwardIterator last,
                       double xeps,
                       double yeps,
                       Random& rnd)
{
    _perturb_points_2( first, last, xeps, yeps, rnd, ch_value_type(first));
}

template <class ForwardIterator>
inline
void perturb_points_2( ForwardIterator first,
                       ForwardIterator last,
                       double xeps,
                       Random& rnd)
{
    perturb_points_2( first, last, xeps, xeps, rnd);
}

template <class ForwardIterator>
void perturb_points_2( ForwardIterator first,
                       ForwardIterator last,
                       double xeps,
                       double yeps)
{
    perturb_points_2( first, last, xeps, yeps, default_random);
}

template <class ForwardIterator>
void perturb_points_2( ForwardIterator first,
                       ForwardIterator last,
                       double xeps)
{
    perturb_points_2( first, last, xeps, xeps, default_random);
}
template <class RandomAccessIterator, class OutputIterator, class Creator>
OutputIterator random_collinear_points_2(
                       RandomAccessIterator first,
                       RandomAccessIterator last,
                       std::size_t n,
                       OutputIterator first2,
                       Random& rnd,
                       Creator creator)
{
    typedef typename Creator::result_type Point;
    int m = last - first;
    for ( std::size_t i = 0; i < n; i++) {
        const Point& p = first[ rnd.get_int( 0, m-1)];
        const Point& q = first[ rnd.get_int( 0, m-1)];
        double la = rnd.get_double();
        double mu = 1.0 - la;
        *first2++ = creator(mu * to_double(p.x()) +
                            la * to_double(q.x()),
                            mu * to_double(p.y()) +
                            la * to_double(q.y()));
    }
    return first2;
}

template <class RandomAccessIterator, class OutputIterator, class P>
OutputIterator _random_collinear_points_2(
                       RandomAccessIterator first,
                       RandomAccessIterator last,
                       std::size_t n,
                       OutputIterator first2,
                       Random& rnd,
                       P*)
{
    return random_collinear_points_2( first, last, n, first2, rnd,
                                      Creator_uniform_2<double,P>());
}

template <class RandomAccessIterator, class OutputIterator>
OutputIterator random_collinear_points_2(
                       RandomAccessIterator first,
                       RandomAccessIterator last,
                       std::size_t n,
                       OutputIterator first2,
                       Random& rnd)
    // choose two random points from the range [`first',`last'), create a
    // random third point on the segment connecting this two points, and
    // write it to `first2'. Repeat this n times, thus writing n points to
    // `first2' that are collinear with points in the range [`first',
    // `last'). Three random numbers are needed from `rnd' for each point.
    // Returns the value of `first2' after inserting the n points.
    // Precondition: The expression `to_double((*first).x()
    // )' and `to_double((*first).y())' must be legal.
{
    return  _random_collinear_points_2( first, last, n, first2, rnd,
                                        ch_value_type(first));
}

template <class RandomAccessIterator, class OutputIterator>
OutputIterator random_collinear_points_2(
                       RandomAccessIterator first,
                       RandomAccessIterator last,
                       std::size_t n,
                       OutputIterator first2)
{
    return  random_collinear_points_2( first, last, n, first2,
                                       default_random);
}



CGAL_END_NAMESPACE
#endif // CGAL_POINT_GENERATORS_2_H //
// EOF //

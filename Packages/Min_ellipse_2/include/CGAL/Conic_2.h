// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-I $
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Conic_2.h
// package       : $CGAL_Package: Min_ellipse_2 $
// chapter       : Geometric Optimisation
//
// source        : web/Conic_2.aw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Bernd Gärtner, Sven Schönherr <sven@inf.ethz.ch>
// coordinator   : ETH Zürich (Bernd Gärtner <gaertner@inf.ethz.ch>)
//
// implementation: 2D Conic
// ============================================================================

#ifndef CGAL_CONIC_2_H
#define CGAL_CONIC_2_H

#ifndef CGAL_REP_CLASS_DEFINED
#  error  no representation class defined
#endif // CGAL_REP_CLASS_DEFINED

#ifdef CGAL_HOMOGENEOUS_H
#  include <CGAL/ConicHPA2.h>
#endif

#ifdef CGAL_CARTESIAN_H
#  include <CGAL/ConicCPA2.h>
#endif

#ifndef CGAL_POINT_2_H
#  include <CGAL/Point_2.h>
#endif

#ifndef CGAL_IO_FORWARD_DECL_WINDOW_STREAM_H
#include <CGAL/IO/forward_decl_window_stream.h>
#endif

CGAL_BEGIN_NAMESPACE

template < class _R >
class Optimisation_ellipse_2;

CGAL_END_NAMESPACE

template < class _R >
CGAL::Window_stream&
operator << ( CGAL::Window_stream&,
              const CGAL::Optimisation_ellipse_2<_R>&);

CGAL_BEGIN_NAMESPACE

template < class _R>
class Conic_2 : public _R::Conic_2 {

    friend  class Optimisation_ellipse_2<_R>;

  public:

    // types
    typedef  _R                    R;
    typedef  typename _R::RT       RT;
    typedef  typename _R::FT       FT;
    typedef  typename _R::Conic_2  _Conic_2;

    // construction
    Conic_2 ()
    {}
    
    Conic_2 (RT r, RT s, RT t, RT u, RT v, RT w)
        : R::Conic_2 (r, s, t, u, v, w)
    {}
    
    

    // general access
    RT r () const
    {
        return _Conic_2::r();
    }
    
    RT s () const
    {
        return _Conic_2::s();
    }
    
    RT t () const
    {
        return _Conic_2::t();
    }
    
    RT u () const
    {
        return _Conic_2::u();
    }
    
    RT v () const
    {
        return _Conic_2::v();
    }
    
    RT w () const
    {
        return _Conic_2::w();
    }
    
    CGAL::Point_2<R> center () const
    {
        return _Conic_2::center();
    }
    
    

    // type related access
    Conic_type conic_type () const
    {
        return _Conic_2::conic_type();
    }
    
    bool is_hyperbola () const
    {
        return _Conic_2::is_hyperbola();
    }
    
    bool is_parabola () const
    {
        return _Conic_2::is_parabola();
    }
    
    bool is_ellipse () const
    {
        return _Conic_2::is_ellipse();
    }
    
    bool is_empty () const
    {
        return _Conic_2::is_empty();
    }
    
    bool is_trivial () const
    {
        return _Conic_2::is_trivial();
    }
    
    bool is_degenerate () const
    {
        return _Conic_2::is_degenerate();
    }
    
    

    // orientation related access
    CGAL::Orientation orientation () const
    {
        return _Conic_2::orientation ();
    }
    
    CGAL::Oriented_side oriented_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::oriented_side (p);
    }
    
    bool has_on_positive_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_positive_side (p);
    }
    
    bool has_on_negative_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_negative_side (p);
    }
    
    bool has_on_boundary (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_boundary (p);
    }
    
    bool has_on (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on (p);
    }
    
    Convex_side convex_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::convex_side (p);
    }
    
    bool has_on_convex_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_convex_side (p);
    }
    
    bool has_on_nonconvex_side (const CGAL::Point_2<R>& p) const
    {
        return _Conic_2::has_on_nonconvex_side (p);
    }
    
    

    // comparisons
    bool operator == ( const Conic_2<_R>& c) const
    {
        return _Conic_2::operator == ( (Conic_2)c);
    }
    
    bool operator != ( const Conic_2<_R>& c) const
    {
        return( ! operator == ( c));
    }

    // set methods
    void set (RT r, RT s, RT t,
              RT u, RT v, RT w)
    {
        _Conic_2::set (r, s, t, u, v, w);
    }
    
    void set_opposite ()
    {
        _Conic_2::set_opposite();
    }
    
    void set_linepair (const CGAL::Point_2<R>& p1, const CGAL::Point_2<R>& p2,
                       const CGAL::Point_2<R>& p3, const CGAL::Point_2<R>& p4)
    {
        _Conic_2::set_linepair (p1, p2, p3, p4);
    }
    
    void set_ellipse (const CGAL::Point_2<R>& p1, const CGAL::Point_2<R>& p2,
                      const CGAL::Point_2<R>& p3)
    {
        _Conic_2::set_ellipse (p1, p2, p3);
    }
    
    void set_ellipse (const CGAL::Point_2<R>& p1, const CGAL::Point_2<R>& p2,
                      const CGAL::Point_2<R>& p3, const CGAL::Point_2<R>& p4,
                      CGAL::Orientation o = POSITIVE)
    {
        _Conic_2::set_ellipse (p1, p2, p3, p4, o);
    }
    
    void set (const CGAL::Point_2<R>& p1, const CGAL::Point_2<R>& p2,
              const CGAL::Point_2<R>& p3, const CGAL::Point_2<R>& p4,
              const CGAL::Point_2<R>& p5,
              CGAL::Orientation o = POSITIVE)
    {
        _Conic_2::set (p1, p2, p3, p4, p5, o);
    }
    
    

  private:
    void set_linear_combination (
        const RT& a1, const Conic_2<R>& c1,
        const RT& a2, const Conic_2<R>& c2)
    {
        _Conic_2::set_linear_combination (a1, c1, a2, c2);
    }
    
    static void set_two_linepairs (const CGAL::Point_2<R>& p1,
                                   const CGAL::Point_2<R>& p2,
                                   const CGAL::Point_2<R>& p3,
                                   const CGAL::Point_2<R>& p4,
                                   Conic_2<R>& pair1,
                                   Conic_2<R>& pair2)
    {
        _Conic_2::set_two_linepairs (p1, p2, p3, p4, pair1, pair2);
    }
    
    void set_ellipse (const Conic_2<R>& pair1,
                      const Conic_2<R>& pair2)
    {
        _Conic_2::set_ellipse (pair1, pair2);
    }
    
    void set (const Conic_2<R>& c1, const Conic_2<R>& c2,
              const CGAL::Point_2<R>& p)
    {
        _Conic_2::set( c1, c2, p);  analyse();
    }
    
    CGAL::Sign vol_derivative (RT dr, RT ds,
                               RT dt, RT du,
                               RT dv, RT dw) const
    {
        return _Conic_2::vol_derivative (dr, ds, dt, du, dv, dw);
    }
    
    double vol_minimum (RT dr, RT ds,
                        RT dt, RT du,
                        RT dv, RT dw) const
    {
        return _Conic_2::vol_minimum (dr, ds, dt, du, dv, dw);
    }
    
    
};



#ifndef CGAL_NO_OSTREAM_INSERT_CONIC_2
template< class _R>
std::ostream& operator << ( std::ostream& os, const Conic_2<_R>& c)
{
    return( os << c.r() << ' ' << c.s() << ' ' << c.t() << ' '
               << c.u() << ' ' << c.v() << ' ' << c.w());
}
#endif // CGAL_NO_OSTREAM_INSERT_CONIC_2

CGAL_END_NAMESPACE

#endif // CGAL_CONIC_2_H

// ===== EOF ==================================================================

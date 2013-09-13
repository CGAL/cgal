// Copyright (c) 2000,2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Bernd Gaertner, Sven Schoenherr <sven@inf.ethz.ch>


#ifndef CGAL_CONIC_2_H
#define CGAL_CONIC_2_H

#include <CGAL/Point_2.h>
#include <CGAL/Conic_misc.h>
#include <CGAL/Dimension.h>

namespace CGAL {

template < class R_ >
class Optimisation_ellipse_2;

template < class R_>
class Conic_2 : public R_::Kernel_base::Conic_2 {

    friend  class Optimisation_ellipse_2<R_>;

  public:

    typedef Dimension_tag<2>  Ambient_dimension;
    typedef Dimension_tag<1>  Feature_dimension;

    // types
    typedef  R_                    R;
    typedef  typename R_::RT       RT;
    typedef  typename R_::FT       FT;
    typedef  typename R_::Kernel_base::Conic_2  _Conic_2;

    // construction
    Conic_2 ()
    {}
    
    Conic_2 (RT r, RT s, RT t, RT u, RT v, RT w)
        : _Conic_2 (r, s, t, u, v, w)
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

    bool is_circle () const
    {
        return _Conic_2::is_circle();
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
    bool operator == ( const Conic_2<R_>& c) const
    {
        return _Conic_2::operator == ( (Conic_2)c);
    }
    
    bool operator != ( const Conic_2<R_>& c) const
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
    
    void set_circle (const CGAL::Point_2<R>& p1, const CGAL::Point_2<R>& p2,
		     const CGAL::Point_2<R>& p3)
    {
        // the unique circle through the three points
        _Conic_2::set_circle(p1, p2, p3);
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
        _Conic_2::set( c1, c2, p);  this->analyse();
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
template< class R_>
std::ostream& operator << ( std::ostream& os, const Conic_2<R_>& c)
{
    return( os << c.r() << ' ' << c.s() << ' ' << c.t() << ' '
               << c.u() << ' ' << c.v() << ' ' << c.w());
}
#endif // CGAL_NO_OSTREAM_INSERT_CONIC_2

} //namespace CGAL

#endif // CGAL_CONIC_2_H

// ===== EOF ==================================================================

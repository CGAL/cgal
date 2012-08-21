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

#ifndef CGAL_CONICHPA2_H
#define CGAL_CONICHPA2_H

// includes
#include <CGAL/Conic_misc.h>
#include <CGAL/kernel_assertions.h>

namespace CGAL {

template < class PT, class DA>
class ConicHPA2;

template < class PT, class DA>
class _Min_ellipse_2_adapterH2__Ellipse;


template < class _PT, class _DA>
class ConicHPA2
{
  public:
    // types
    typedef           _PT      PT;
    typedef           _DA      DA;
    typedef  typename _DA::RT  RT;

  //private:
    //friend class Conic_2< CGAL::Homogeneous<RT> >;
    friend class _Min_ellipse_2_adapterH2__Ellipse<PT,DA>;

    DA                  dao;
    RT                  _r, _s, _t, _u, _v, _w;
    Conic_type          type;
    CGAL::Orientation   o;
    bool                empty, trivial, degenerate;
    
    
    void
    set_linear_combination (const RT& a1, const ConicHPA2<PT,DA>& c1,
                            const RT& a2, const ConicHPA2<PT,DA>& c2)
    {
        _r = a1 * c1.r() + a2 * c2.r();
        _s = a1 * c1.s() + a2 * c2.s();
        _t = a1 * c1.t() + a2 * c2.t();
        _u = a1 * c1.u() + a2 * c2.u();
        _v = a1 * c1.v() + a2 * c2.v();
        _w = a1 * c1.w() + a2 * c2.w();
    }
    
    static void set_two_linepairs (const PT& p1,
                                   const PT& p2,
                                   const PT& p3,
                                   const PT& p4,
                                   ConicHPA2<PT,DA>& pair1,
                                   ConicHPA2<PT,DA>& pair2)
    {
        RT x1, y1, h1, x2, y2, h2, x3, y3, h3, x4, y4, h4;
        const DA& da = pair1.da();
        da.get (p1, x1, y1, h1);
        da.get (p2, x2, y2, h2);
        da.get (p3, x3, y3, h3);
        da.get (p4, x4, y4, h4);
    
        CGAL::Orientation side1_24 = (CGAL::Orientation)(CGAL_NTS sign
                                       (-h2*x1*y4+h1*x2*y4
                                        +h2*x4*y1-h4*x2*y1
                                        +h4*x1*y2-h1*x4*y2)),
                         side3_24 = (CGAL::Orientation)(CGAL_NTS sign
                                      (-h2*x3*y4+h3*x2*y4
                                       +h2*x4*y3-h4*x2*y3
                                       +h4*x3*y2-h3*x4*y2));
        if (side1_24 != side3_24) {
            // (counter)clockwise order
            pair1.set_linepair (p1, p2, p3, p4);
            pair2.set_linepair (p2, p3, p4, p1);
        } else {
            CGAL::Orientation side1_32 = (CGAL::Orientation)(CGAL_NTS sign
                                           (-h3*x1*y2+h1*x3*y2
                                            +h3*x2*y1-h2*x3*y1
                                            +h2*x1*y3-h1*x2*y3));
            if (side1_32 != side3_24) {
                // p1, p2 need to be swapped
                pair1.set_linepair (p2, p1, p3, p4);
                pair2.set_linepair (p1, p3, p4, p2);
            } else {
                // p2, p3 need to be swapped
                pair1.set_linepair (p1, p3, p2, p4);
                pair2.set_linepair (p3, p2, p4, p1);
            }
        }
    }
    
    void set_ellipse (const ConicHPA2<PT,DA>& pair1,
                      const ConicHPA2<PT,DA>& pair2)
    {
        RT b = RT(2) * (pair1.r() * pair2.s() + pair1.s() * pair2.r()) -
               pair1.t() * pair2.t();
        set_linear_combination (pair2.det()-b, pair1,
                                pair1.det()-b, pair2);
    }
    
    void set (const ConicHPA2<PT,DA>& c1,
              const ConicHPA2<PT,DA>& c2,
              const PT& p)
    {
        set_linear_combination (c2.evaluate(p), c1, -c1.evaluate(p), c2);
    }
    
    CGAL::Sign vol_derivative (RT dr, RT ds, RT dt,
                              RT du, RT dv, RT dw) const
    {
        RT a1 = RT(4)*r()*ds+RT(4)*dr*s()-RT(2)*t()*dt,
           a0 = RT(4)*r()*s()-t()*t(),
           b1 = (RT(4)*r()*s()-t()*t())*dw+(RT(4)*r()*ds+RT(4)*dr*s()-
                RT(2)*t()*dt)*w()-u()*u()*ds -
                RT(2)*u()*du*s()-v()*v()*dr-RT(2)*v()*dv*r()+u()*v()*dt+
                (u()*dv+du*v())*t(),
           b0 = (RT(4)*r()*s()-t()*t())*w()
                -u()*u()*s()-v()*v()*r()+u()*v()*t(),
           c0 = -RT(2)*a0*b1 + RT(3)*a1*b0;
    
        return CGAL_NTS sign ((int)-CGAL_NTS sign (c0)*o);
    }
    
    double vol_minimum (RT dr, RT ds, RT dt, RT du, RT dv, RT dw) const
    {
        RT a2 = RT(4)*dr*ds-dt*dt,
           a1 = RT(4)*r()*ds+RT(4)*dr*s()-RT(2)*t()*dt,
           a0 = RT(4)*r()*s()-t()*t(),
           b3 = (RT(4)*dr*ds-dt*dt)*dw-du*du*ds-dv*dv*dr+du*dv*dt,
           b2 = (RT(4)*r()*ds+RT(4)*dr*s()-RT(2)*t()*dt)*dw+
                (RT(4)*dr*ds-dt*dt)*w()-RT(2)*u()*du*ds-du*du*s()-
                RT(2)*v()*dv*dr-dv*dv*r()+(u()*dv+du*v())*dt+du*dv*t(),
           b1 = (RT(4)*r()*s()-t()*t())*dw+(RT(4)*r()*ds+RT(4)*dr*s()-
                RT(2)*t()*dt)*w()-u()*u()*ds -
                RT(2)*u()*du*s()-v()*v()*dr-RT(2)*v()*dv*r()+u()*v()*dt+
                (u()*dv+du*v())*t(),
           b0 = (RT(4)*r()*s()-t()*t())*w()
                -u()*u()*s()-v()*v()*r()+u()*v()*t(),
           c3 = -RT(3)*a1*b3 + RT(2)*a2*b2,
           c2 = -RT(6)*a0*b3 - a1*b2 + RT(4)*a2*b1,
           c1 = -RT(4)*a0*b2 + a1*b1 + RT(6)*a2*b0,
           c0 = -RT(2)*a0*b1 + RT(3)*a1*b0;
    
           double roots[3];
           int nr_roots = solve_cubic
                                (CGAL::to_double(c3), CGAL::to_double(c2),
                                 CGAL::to_double(c1), CGAL::to_double(c0),
                                 roots[0], roots[1], roots[2]);
           CGAL_kernel_precondition (nr_roots > 0); // minimum exists
           return best_value (roots, nr_roots,
                                 CGAL::to_double(a2), CGAL::to_double(a1),
                                 CGAL::to_double(a0), CGAL::to_double(b3),
                                 CGAL::to_double(b2), CGAL::to_double(b1),
                                 CGAL::to_double(b0));
    }
    
    

  protected:
    RT det () const
    {
        return RT(4)*s()*r() - t()*t();
    }
    
    void analyse( )
    {
        RT d = det();
        type = (Conic_type)(CGAL_NTS sign(d));
        switch (type) {
        case HYPERBOLA:
            {
                trivial = empty = false;
                RT z_prime = d*w() - u()*u()*s() - v()*v()*r() + u()*v()*t();
                o = (CGAL::Orientation)(CGAL_NTS sign (z_prime));
                degenerate = (o == CGAL::ZERO);
                
                
            }
            break;
        case PARABOLA:
            {
                if (!CGAL_NTS is_zero (r())) {
                    trivial         = false;
                    degenerate      = (t()*u() == RT(2)*r()*v());
                    if (degenerate) {
                        CGAL::Sign discr = (CGAL::Sign)
                                        CGAL_NTS sign(u()*u()-RT(4)*r()*w());
                        switch (discr) {
                            case CGAL::NEGATIVE:
                                empty = true;
                                o = (CGAL::Orientation)(CGAL_NTS sign (w()));
                                break;
                            case CGAL::ZERO:
                                empty = false;
                                o = (CGAL::Orientation)(CGAL_NTS sign (r()));
                                break;
                            case CGAL::POSITIVE:
                                empty = false;
                                o = CGAL::ZERO;
                                break;
                        }
                    } else {
                        empty = false;
                        o = (CGAL::Orientation)(-CGAL_NTS sign (r()));
                    }
                } else if (!CGAL_NTS is_zero (s())) {
                    trivial         = false;
                    degenerate      = (t()*v() == RT(2)*s()*u());
                    if (degenerate) {
                        CGAL::Sign discr = (CGAL::Sign)
                                        CGAL_NTS sign(v()*v()-RT(4)*s()*w());
                        switch (discr) {
                            case CGAL::NEGATIVE:
                                empty = true;
                                o = (CGAL::Orientation)(CGAL_NTS sign (w()));
                                break;
                            case CGAL::ZERO:
                                empty = false;
                                o = (CGAL::Orientation)(CGAL_NTS sign (s()));
                                break;
                            case CGAL::POSITIVE:
                                empty = false;
                                o = CGAL::ZERO;
                                break;
                        }
                    } else {
                        empty = false;
                        o = (CGAL::Orientation)(-CGAL_NTS sign (s()));
                    }
                } else { // r=0, s=0
                    degenerate      = true;
                    bool uv_zero    =    CGAL_NTS is_zero (u())
                                      && CGAL_NTS is_zero (v());
                    trivial         = uv_zero && CGAL_NTS is_zero (w());
                    empty           = uv_zero && !trivial;
                    if (empty)
                        o = (CGAL::Orientation)(CGAL_NTS sign (w()));
                    else if (trivial)
                        o = CGAL::POSITIVE;
                    else
                        o = CGAL::ZERO;
                }
                
                
            }
            break;
        case ELLIPSE:
            {
                trivial = false;
                RT z_prime = d*w() - u()*u()*s() - v()*v()*r() + u()*v()*t();
                if (CGAL_NTS is_positive (r())) {
                    empty = (CGAL::POSITIVE == CGAL_NTS sign (z_prime));
                    empty ? o = CGAL::POSITIVE : o = CGAL::NEGATIVE;
                } else {
                    empty = (CGAL::NEGATIVE == CGAL_NTS sign (z_prime));
                    empty ? o = CGAL::NEGATIVE : o = CGAL::POSITIVE;
                }
                degenerate = empty || CGAL_NTS is_zero (z_prime);
            }
            break;
        }
    }
    
    RT evaluate (const PT& p) const
    {
        RT x, y, h;
        dao.get (p, x, y, h);
        return  r()*x*x + s()*y*y + t()*x*y + u()*x*h + v()*y*h + w()*h*h;
    }
    
    

  public:
    ConicHPA2 ( const DA& da = DA()) : dao( da) { }
    
    ConicHPA2 (RT r, RT s, RT t, RT u, RT v, RT w, const DA& da = DA())
        : dao( da), _r(r), _s(s), _t(t), _u(u), _v(v), _w(w)
    {
        analyse();
    }
    
    const DA&  da() const
    {
        return dao;
    }
    
    RT r() const { return _r;}
    RT s() const { return _s;}
    RT t() const { return _t;}
    RT u() const { return _u;}
    RT v() const { return _v;}
    RT w() const { return _w;}
    
    PT center () const
    {
        CGAL_kernel_precondition (type != PARABOLA);
        // PT p;
	// replaced previous line by following hack (no idea
	// why original version doesn't work)
        typename DA::Point_2 p;
        RT two = RT(2);
        dao.set( p, two*s()*u() - t()*v(), two*r()*v() - t()*u(), -det());
        return p;
    }
    
    Conic_type conic_type () const
    {
        return type;
    }
    
    bool is_hyperbola () const
    {
        return (type == HYPERBOLA);
    }
    
    bool is_parabola () const
    {
        return (type == PARABOLA);
    }
    
    bool is_ellipse () const
    {
        return (type == ELLIPSE);
    }
    
    bool is_circle () const
    {
        return (type == ELLIPSE && (r()==s()) && CGAL_NTS is_zero (t()));
    }
   
    bool is_empty () const
    {
        return empty;
    }
    
    bool is_trivial () const
    {
        return trivial;
    }
    
    bool is_degenerate () const
    {
        return degenerate;
    }
    
    CGAL::Orientation orientation () const
    {
        return o;
    }
    
    CGAL::Oriented_side oriented_side (const PT& p) const
    {
        return (CGAL::Oriented_side)(CGAL_NTS sign (evaluate (p)));
    }
    
    bool has_on_positive_side (const PT& p) const
    {
        return (CGAL_NTS is_positive (evaluate(p)));
    }
    
    bool has_on_negative_side (const PT& p) const
    {
        return (CGAL_NTS is_negative (evaluate(p)));
    }
    
    bool has_on_boundary (const PT& p) const
    {
       return (CGAL_NTS is_zero (evaluate(p)));
    }
    
    bool has_on (const PT& p) const
    {
       return (CGAL_NTS is_zero (evaluate(p)));
    }
    
    Convex_side convex_side (const PT& p) const
    {
        switch (o) {
        case CGAL::POSITIVE:
            return (Convex_side)( CGAL_NTS sign (evaluate (p)));
        case CGAL::NEGATIVE:
            return (Convex_side)(-CGAL_NTS sign (evaluate (p)));
        case CGAL::ZERO:
            return (Convex_side)( CGAL_NTS sign (CGAL_NTS abs (evaluate(p))));
        }
        // keeps g++ happy
        return( Convex_side( 0));
    }
    
    bool has_on_convex_side (const PT& p) const
    {
        return (convex_side (p) == ON_CONVEX_SIDE);
    }
    
    bool has_on_nonconvex_side (const PT& p) const
    {
        return (convex_side (p) == ON_NONCONVEX_SIDE);
    }
    
    bool operator == ( const ConicHPA2<_PT,_DA>& c) const
    {
        // find coefficient != 0
        RT  factor1(0);
        if ( ! CGAL_NTS is_zero( r())) factor1 = r(); else
        if ( ! CGAL_NTS is_zero( s())) factor1 = s(); else
        if ( ! CGAL_NTS is_zero( t())) factor1 = t(); else
        if ( ! CGAL_NTS is_zero( u())) factor1 = u(); else
        if ( ! CGAL_NTS is_zero( v())) factor1 = v(); else
        if ( ! CGAL_NTS is_zero( w())) factor1 = w(); else
        CGAL_kernel_assertion_msg( false, "all coefficients zero");
    
        // find coefficient != 0
        RT  factor2(0);
        if ( ! CGAL_NTS is_zero( c.r())) factor2 = c.r(); else
        if ( ! CGAL_NTS is_zero( c.s())) factor2 = c.s(); else
        if ( ! CGAL_NTS is_zero( c.t())) factor2 = c.t(); else
        if ( ! CGAL_NTS is_zero( c.u())) factor2 = c.u(); else
        if ( ! CGAL_NTS is_zero( c.v())) factor2 = c.v(); else
        if ( ! CGAL_NTS is_zero( c.w())) factor2 = c.w(); else
        CGAL_kernel_assertion_msg( false, "all coefficients zero");
    
        return(    ( r()*factor2 == c.r()*factor1)
                && ( s()*factor2 == c.s()*factor1)
                && ( t()*factor2 == c.t()*factor1)
                && ( u()*factor2 == c.u()*factor1)
                && ( v()*factor2 == c.v()*factor1)
                && ( w()*factor2 == c.w()*factor1));
    }
    
    void set (RT r_, RT s_, RT t_, RT u_, RT v_, RT w_)
    {
        _r = r_; _s = s_; _t = t_; _u = u_; _v = v_; _w = w_;
        analyse();
     }
    
    void set_opposite ()
    {
        _r = -r(); _s = -s(); _t = -t(); _u = -u(); _v = -v(); _w = -w();
        o = CGAL::opposite(orientation());
    }
     
  void set_circle (const PT& p1, const PT& p2, const PT& p3) 
  { 
     // the circle will have r = s = det*h1*h2*h3, t=0
     RT x1, y1, h1, x2, y2, h2, x3, y3, h3;
     dao.get (p1, x1, y1, h1);
     dao.get (p2, x2, y2, h2);
     dao.get (p3, x3, y3, h3);
    
     // precondition: p1, p2, p3 not collinear
     RT det = -h1*x3*y2+h3*x1*y2+h1*x2*y3-h2*x1*y3+h2*x3*y1-h3*x2*y1;
     CGAL_kernel_precondition (!CGAL_NTS is_zero (det));
     
     // Cramer's rule
     RT sqr1 = (-x1*x1 - y1*y1)*h2*h3;
     RT sqr2 = (-x2*x2 - y2*y2)*h1*h3;
     RT sqr3 = (-x3*x3 - y3*y3)*h1*h2;

     _u = -h1*sqr3*y2+h3*sqr1*y2+h1*sqr2*y3-h2*sqr1*y3+h2*sqr3*y1-h3*sqr2*y1;
     _v = -h1*x3*sqr2+h3*x1*sqr2+h1*x2*sqr3-h2*x1*sqr3+h2*x3*sqr1-h3*x2*sqr1;
     _w = -sqr1*x3*y2+sqr3*x1*y2+sqr1*x2*y3-sqr2*x1*y3+sqr2*x3*y1-sqr3*x2*y1;
     _r = det*h1*h2*h3;
     _s = _r;
     _t = RT(0);

     analyse();
     CGAL_kernel_postcondition(is_circle());
     CGAL_kernel_postcondition(has_on_boundary(p1));
     CGAL_kernel_postcondition(has_on_boundary(p2));
     CGAL_kernel_postcondition(has_on_boundary(p3));
  }
 
    void set_linepair (const PT& p1, const PT& p2, const PT& p3,
                       const PT& p4, const DA& da = DA())
    {
        RT x1, y1, h1, x2, y2, h2, x3, y3, h3, x4, y4, h4;
        da.get (p1, x1, y1, h1);
        da.get (p2, x2, y2, h2);
        da.get (p3, x3, y3, h3);
        da.get (p4, x4, y4, h4);
    
        // precondition: p1 != p2, p3 != p4
        CGAL_kernel_precondition
            ( ((x1*h2 != x2*h1) || (y1*h2 != y2*h1)) &&
              ((x3*h4 != x4*h3) || (y3*h4 != y4*h3)) );
    
        RT h1x2_x1h2 = h1*x2-x1*h2;
        RT h3x4_x3h4 = h3*x4-x3*h4;
        RT y1h2_h1y2 = y1*h2-h1*y2;
        RT y3h4_h3y4 = y3*h4-h3*y4;
        RT x1y2_y1x2 = x1*y2-y1*x2;
        RT x3y4_y3x4 = x3*y4-y3*x4;
    
        _r = y1h2_h1y2 * y3h4_h3y4;
        _s = h1x2_x1h2 * h3x4_x3h4;
        _t = h1x2_x1h2 * y3h4_h3y4 + y1h2_h1y2 * h3x4_x3h4;
        _u = x1y2_y1x2 * y3h4_h3y4 + y1h2_h1y2 * x3y4_y3x4;
        _v = x1y2_y1x2 * h3x4_x3h4 + h1x2_x1h2 * x3y4_y3x4;
        _w = x1y2_y1x2 * x3y4_y3x4;
    
        analyse();
    }
    
    void set_ellipse (const PT& p1, const PT& p2, const PT& p3)
    {
        RT x1, y1, h1, x2, y2, h2, x3, y3, h3;
        dao.get (p1, x1, y1, h1);
        dao.get (p2, x2, y2, h2);
        dao.get (p3, x3, y3, h3);
    
        // precondition: p1, p2, p3 not collinear
        RT det = -h1*x3*y2+h3*x1*y2+h1*x2*y3-h2*x1*y3+h2*x3*y1-h3*x2*y1;
        CGAL_kernel_precondition (!CGAL_NTS is_zero (det));
    
        RT x1x1 = x1*x1, y1y1 = y1*y1,
           x2x2 = x2*x2, y2y2 = y2*y2,
           x3x3 = x3*x3, y3y3 = y3*y3,  // x_i^2, y_i^2
           x1h1 = x1*h1, y1h1 = y1*h1,
           x2h2 = x2*h2, y2h2 = y2*h2,
           x3h3 = x3*h3, y3h3 = y3*h3,  // x_i h_i, y_i h_i
           h1h1 = h1*h1,
           h2h2 = h2*h2,
           h3h3 = h3*h3,                // h_i^2
           two = RT(2);                 // 2
    
        _r = y1y1*h2h2*h3h3 - y1h1*y2h2*h3h3 - y1h1*h2h2*y3h3 +
             h1h1*y2y2*h3h3 - h1h1*y2h2*y3h3 + h1h1*h2h2*y3y3;
    
        _s = x1x1*h2h2*h3h3 - x1h1*x2h2*h3h3 - x1h1*h2h2*x3h3 +
             h1h1*x2x2*h3h3 - h1h1*x2h2*x3h3 + h1h1*h2h2*x3x3;
    
        _t = -two*x1*y1*h2h2*h3h3 + x1h1*y2h2*h3h3 + x1h1*h2h2*y3h3 +
                 y1h1*x2h2*h3h3 -two*h1h1*x2*y2*h3h3 + h1h1*x2h2*y3h3 +
                 y1h1*h2h2*x3h3 + h1h1*y2h2*x3h3 -two*h1h1*h2h2*x3*y3;
    
        _u = -(h1h1*y2y2*x3h3 - h1h1*x2*y2*y3h3 - h1h1*y2h2*x3*y3 +
                   x1h1*h2h2*y3y3 + h1h1*x2h2*y3y3 +y1y1*x2h2*h3h3 +
                   y1y1*h2h2*x3h3 - x1*y1*y2h2*h3h3 - y1h1*x2*y2*h3h3 -
                   x1*y1*h2h2*y3h3 - y1h1*h2h2*x3*y3 + x1h1*y2y2*h3h3);
    
        _v = -(h1h1*x2x2*y3h3 - h1h1*x2*y2*x3h3 - h1h1*x2h2*x3*y3 +
                   y1h1*h2h2*x3x3 + h1h1*y2h2*x3x3 +x1x1*y2h2*h3h3 +
                   x1x1*h2h2*y3h3 - x1*y1*x2h2*h3h3 - x1h1*x2*y2*h3h3 -
                   x1*y1*h2h2*x3h3 - x1h1*h2h2*x3*y3 + y1h1*x2x2*h3h3);
    
        _w = y1y1*x2h2*x3h3 - x1*y1*y2h2*x3h3 - y1h1*x2*y2*x3h3 +
             y1h1*y2h2*x3x3 - x1*y1*x2h2*y3h3 + y1h1*x2x2*y3h3 -
             y1h1*x2h2*x3*y3 + x1h1*y2y2*x3h3 + x1x1*y2h2*y3h3 -
             x1h1*x2*y2*y3h3 - x1h1*y2h2*x3*y3 + x1h1*x2h2*y3y3;
    
        type = ELLIPSE;
        degenerate = trivial = empty = false;
        o = CGAL::NEGATIVE;
        if (CGAL_NTS is_positive (det)) set_opposite ();
    
    }
    
    void set_ellipse (const PT& p1, const PT& p2,
                      const PT& p3, const PT& p4,
                      CGAL::Orientation _o = POSITIVE)
    {
        ConicHPA2<PT,DA> pair1, pair2;
        set_two_linepairs (p1, p2, p3, p4, pair1, pair2);
        set_ellipse (pair1, pair2);
        analyse();
        if (o != _o) set_opposite();
    }
    
    void set (const PT& p1, const PT& p2, const PT& p3, const PT& p4,
              const PT& p5, CGAL::Orientation _o = POSITIVE)
    {
        ConicHPA2<PT,DA> c1; c1.set_linepair (p1, p2, p3, p4);
        ConicHPA2<PT,DA> c2; c2.set_linepair (p1, p4, p2, p3);
        set_linear_combination (c2.evaluate (p5), c1,
                               -c1.evaluate (p5), c2);
        analyse();
        // precondition: all points distinct <=> conic nontrivial
        CGAL_kernel_precondition (!is_trivial());
        if (o != _o) set_opposite();
    }
    
    
};

#ifndef CGAL_NO_OSTREAM_INSERT_CONICHPA2
template< class _PT, class _DA>
std::ostream& operator << ( std::ostream& os, const ConicHPA2<_PT,_DA>& c)
{
    return( os << c.r() << ' ' << c.s() << ' ' << c.t() << ' '
               << c.u() << ' ' << c.v() << ' ' << c.w());
}

template< class _PT, class _DA>
std::istream& operator >> ( std::istream& is, ConicHPA2<_PT,_DA>& c)
{
    typedef  typename _DA::RT                  RT;

    RT  r, s, t, u, v, w;
    is >> r >> s >> t >> u >> v >> w;
    c.set( r, s, t, u, v, w);

    return( is);
}
#endif // CGAL_NO_OSTREAM_INSERT_CONICHPA2

} //namespace CGAL

#endif // CGAL_CONICHPA2_H

// ===== EOF ==================================================================

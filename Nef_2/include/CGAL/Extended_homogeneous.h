// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_EXTENDED_HOMOGENEOUS_H
#define CGAL_EXTENDED_HOMOGENEOUS_H

#include <CGAL/license/Nef_2.h>


#include <CGAL/basic.h>
#include <CGAL/Homogeneous.h> 
#include <CGAL/Point_2.h> 
#include <CGAL/Line_2_Line_2_intersection.h> 
#include <CGAL/squared_distance_2.h> 
#include <CGAL/number_utils.h>
#include <CGAL/Nef_polynomial.h>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 5
#include <CGAL/Nef_2/debug.h>
#include <CGAL/Nef_2/Line_to_epoint.h>
#include <CGAL/Is_extended_kernel.h>

namespace CGAL {


template <class T> class Extended_homogeneous;

template<class T>
struct Is_extended_kernel<Extended_homogeneous<T> > {
       typedef Tag_true value_type;
};

/*{\Moptions outfile=ExtendedKernelTraits_2.man}*/
/*{\Moptions print_title=yes }*/ 
/*{\Msubst Extended_homogeneous ExtendedKernelTraits_2}*/ 
/*{\Manpage {ExtendedKernelTraits_2}{}{Extended Kernel Traits}{K}}*/

template <class RT_>
class Extended_homogeneous : public 
  CGAL::Homogeneous< CGAL::Nef_polynomial<RT_> > { public:

/*{\Mdefinition |\Mname| is a kernel concept providing extended
geometry\cgalFootnote{It is called extended geometry for simplicity,
though it is not a real geometry in the classical sense.}. Let |\Mvar|
be an instance of the data type |\Mname|.  The central notion of
extended geometry are extended points. An extended point represents
either a standard affine point of the Cartesian plane or a
non-standard point representing the equivalence class of rays where
two rays are equivalent if one is contained in the other.

Let $R$ be an infinimaximal number\cgalFootnote{A finite but very large
number.}, $F$ be the square box with corners $NW(-R,R)$, $NE(R,R)$,
$SE(R,-R)$, and $SW(-R,-R)$. Let $p$ be a non-standard point and let
$r$ be a ray defining it. If the frame $F$ contains the source point
of $r$ then let $p(R)$ be the intersection of $r$ with the frame $F$,
if $F$ does not contain the source of $r$ then $p(R)$ is undefined.
For a standard point let $p(R)$ be equal to $p$ if $p$ is contained in
the frame $F$ and let $p(R)$ be undefined otherwise. Clearly, for any
standard or non-standard point $p$, $p(R)$ is defined for any
sufficiently large $R$. Let $f$ be any function on standard points,
say with $k$ arguments. We call $f$ {\em extensible} if for any $k$
points $p_1$, \ldots, $p_k$ the function value
$f(p_1(R),\ldots,p_k(R))$ is constant for all sufficiently large
$R$. We define this value as $f(p_1,\ldots,p_k)$.  Predicates like
lexicographic order of points, orientation, and incircle tests are
extensible.

An extended segment is defined by two extended points such that it
is either an affine segment, an affine ray, an affine line, or a
segment that is part of the square box. Extended directions extend
the affine notion of direction to extended objects.

This extended geometry concept serves two purposes. It offers
functionality for changing between standard affine and extended
geometry. At the same time it provides extensible geometric primitives
on the extended geometric objects.}*/

  typedef CGAL::Homogeneous< CGAL::Nef_polynomial<RT_> >  Base;
  typedef Extended_homogeneous<RT_> Self;

  /*{\Mtypes 8.5}*/
  /*{\Mtext \headerline{Affine kernel types}}*/

  typedef CGAL::Homogeneous<RT_> Standard_kernel;
  /*{\Mtypemember the standard affine kernel.}*/

  typedef RT_ Standard_RT;
  /*{\Mtypemember the standard ring type.}*/

  typedef typename Standard_kernel::FT Standard_FT;
  /*{\Xtypemember the field type.}*/

  typedef typename Standard_kernel::Point_2     Standard_point_2;
  /*{\Mtypemember standard points.}*/

  typedef typename Standard_kernel::Segment_2   Standard_segment_2;
  /*{\Mtypemember standard segments.}*/

  typedef typename Standard_kernel::Line_2      Standard_line_2;
  /*{\Mtypemember standard oriented lines.}*/

  typedef typename Standard_kernel::Direction_2 Standard_direction_2;
  /*{\Mtypemember standard directions.}*/

  typedef typename Standard_kernel::Ray_2       Standard_ray_2;
  /*{\Mtypemember standard rays.}*/

  typedef typename Standard_kernel::Aff_transformation_2 
    Standard_aff_transformation_2;
  /*{\Mtypemember standard affine transformations.}*/

  /*{\Mtext \headerline{Extended kernel types}}*/

  typedef typename Base::RT RT;
  /*{\Mtypemember the ring type of our extended kernel.}*/

  typedef typename Base::Point_2      Point_2;
  /*{\Mtypemember extended points.}*/

  typedef typename Base::Segment_2    Segment_2;
  /*{\Mtypemember extended segments.}*/

  typedef typename Base::Direction_2  Direction_2;
  /*{\Mtypemember extended directions.}*/

  typedef typename Base::Line_2       Line_2;
  // used only internally

  enum Point_type { SWCORNER=1, LEFTFRAME, NWCORNER, 
                    BOTTOMFRAME, STANDARD, TOPFRAME,
                    SECORNER, RIGHTFRAME, NECORNER };
  /*{\Menum a type descriptor for extended points.}*/


  public:
  Point_2 epoint(const Standard_RT& m1, const Standard_RT& n1, 
                 const Standard_RT& m2, const Standard_RT& n2, 
                                const Standard_RT& n3) const
  { return Point_2(RT(n1,m1),RT(n2,m2),RT(n3)); }

  Point_2 construct_point(const Standard_line_2& l, Point_type& t) const
  { 
    t = (Point_type)Line_to_epoint<Standard_kernel>::determine_type(l);
    Point_2 res;
    switch (t) {
      case SWCORNER:   res = epoint(-1, 0, -1, 0, 1); break;
      case NWCORNER:   res = epoint(-1, 0,  1, 0, 1); break;
      case SECORNER:   res = epoint( 1, 0, -1, 0, 1); break; 
      case NECORNER:   res = epoint( 1, 0,  1, 0, 1); break;  
      case LEFTFRAME:  
        res = epoint(-l.b(), 0,  l.a(), -l.c(), l.b()); break; 
      case RIGHTFRAME: 
        res = epoint( l.b(), 0, -l.a(), -l.c(), l.b()); break; 
      case BOTTOMFRAME: 
        res = epoint( l.b(), -l.c(), -l.a(), 0, l.a()); break; 
      case TOPFRAME:
        res = epoint(-l.b(), -l.c(),  l.a(), 0, l.a()); break; 
      default: CGAL_error_msg("EPoint type not correct!");
    }
    return res;
  }

  template <class Forward_iterator>
  void determine_frame_radius(Forward_iterator start, Forward_iterator end,
                              Standard_RT& R0) const
  { Standard_RT R, mx, nx, my, ny;
    while ( start != end ) {
      Point_2 p = *start++;
      if ( is_standard(p) ) {
        R = (CGAL::max)(CGAL_NTS abs(p.hx()[0])/p.hw()[0], 
			CGAL_NTS abs(p.hy()[0])/p.hw()[0]);
      } else {
        RT rx = CGAL_NTS abs(p.hx()), ry = CGAL_NTS abs(p.hy());
        mx = ( rx.degree()>0 ? rx[1] : Standard_RT(0) ); nx = rx[0];
        my = ( ry.degree()>0 ? ry[1] : Standard_RT(0) ); ny = ry[0];
        if ( mx > my )      R = CGAL_NTS abs((ny-nx)/(mx-my));
        else if ( mx < my ) R = CGAL_NTS abs((nx-ny)/(my-mx));
        else /* mx == my */ R = CGAL_NTS abs(nx-ny)/(2*p.hw()[0]);
      }
      R0 = (CGAL::max)(R+1,R0);
    }
  }


  /*{\Moperations 2}*/
  /*{\Mtext \headerline{Interfacing the affine kernel types}}*/

  Point_2 construct_point(const Standard_point_2& p) const
  /*{\Mop creates an extended point and initializes it to the 
  standard point |p|.}*/
  { return Point_2(p.hx(), p.hy(), p.hw()); }

  Point_2 construct_point(const Standard_point_2& p1, 
                          const Standard_point_2& p2, 
                          Point_type& t) const
  /*{\Xop creates an extended point and initializes it to the equivalence
  class of all the rays underlying the oriented line |l(p1,p2)|. 
  |t| returns the type of the new extended point.}*/
  { return construct_point(Standard_line_2(p1,p2),t); }

  Point_2 construct_point(const Standard_line_2& l) const
  /*{\Mop creates an extended point and initializes it to the equivalence
  class of all the rays underlying the oriented line |l|. }*/
  { Point_type dummy; return construct_point(l,dummy); }

  Point_2 construct_point(const Standard_point_2& p1, 
                          const Standard_point_2& p2) const
  /*{\Mop creates an extended point and initializes it to the equivalence
  class of all the rays underlying the oriented line |l(p1,p2)|.}*/
  { return construct_point(Standard_line_2(p1,p2)); }

  Point_2 construct_point(const Standard_point_2& p, 
                          const Standard_direction_2& d) const
  /*{\Mop creates an extended point and initializes it to the equivalence
  class of all the rays underlying the ray starting in |p| in direction |d|.}*/
  { return construct_point(Standard_line_2(p,d)); }

  Point_2 construct_opposite_point(const Standard_line_2& l) const
  /*{\Mop creates an extended point and initializes it to the equivalence
  class of all the rays underlying the oriented line opposite to |l|. }*/
  { Point_type dummy; return construct_point(l.opposite(),dummy); }


  Point_type type(const Point_2& p) const
  /*{\Mop determines the type of |p| and returns it.}*/
  {
    CGAL_assertion(p.hx().degree()>=0 && p.hy().degree()>=0 );
    CGAL_assertion(p.hw().degree()==0);
    if (p.hx().degree() == 0 && p.hy().degree() == 0) 
      return STANDARD;
    // now we are on the square frame
    RT rx = p.hx();
    RT ry = p.hy();
    int sx = CGAL_NTS sign(rx);
    int sy = CGAL_NTS sign(ry);
    if (sx < 0) rx = -rx;
    if (sy < 0) ry = -ry;
    if (rx>ry) {
      if (sx > 0) return RIGHTFRAME;
      else        return LEFTFRAME;
    }
    if (rx<ry) {
      if (sy > 0) return TOPFRAME;
      else        return BOTTOMFRAME;
    }
    // now (rx == ry) 
    if (sx==sy) {
      if (sx < 0) return SWCORNER;
      else        return NECORNER;
    } else { CGAL_assertion(sx==-sy);
      if (sx < 0) return NWCORNER;
      else        return SECORNER;
    }
  }


  bool is_standard(const Point_2& p) const
  /*{\Mop returns |true| iff |p| is a standard point.}*/
  { return (type(p)==STANDARD);  }

  Standard_point_2 standard_point(const Point_2& p) const
  /*{\Mop returns the standard point represented by |p|.
  \precond |\Mvar.is_standard(p)|.}*/
  { CGAL_assertion(type(p)==STANDARD);
    CGAL_assertion(p.hw() > RT(0));
    return Standard_point_2(p.hx()[0],p.hy()[0],p.hw()[0]);
  }

  Standard_line_2 standard_line(const Point_2& p) const
  /*{\Mop returns the oriented line representing the 
  bundle of rays defining |p|.
  \precond |!\Mvar.is_standard(p)|.}*/
  { CGAL_assertion(type(p)!=STANDARD);
    RT hx = p.hx(), hy = p.hy(), hw = p.hw();
    Standard_RT dx,dy;
    if (hx.degree()>0) dx=hx[1]; else dx=0;
    if (hy.degree()>0) dy=hy[1]; else dy=0;
    Standard_point_2 p0(hx[0],hy[0],hw[0]);
    Standard_point_2 p1(hx[0]+dx,hy[0]+dy,hw[0]);
    return Standard_line_2(p0,p1);
  }

  Standard_ray_2 standard_ray(const Point_2& p) const
  /*{\Mop a ray defining |p|. \precond |!\Mvar.is_standard(p)|.}*/
  { CGAL_assertion(type(p)!=STANDARD);
    Standard_line_2 l = standard_line(p);
    Standard_direction_2 d = l.direction();
    Standard_point_2 q = l.point(0);
    return Standard_ray_2(q,d);
  }

  Point_2 NE() const { return construct_point(Standard_line_2(-1, 1,0)); }
  /*{\Mop returns the point on the northeast frame corner.}*/

  Point_2 SE() const { return construct_point(Standard_line_2( 1, 1,0)); }
  /*{\Mop returns the point on the southeast frame corner.}*/

  Point_2 NW() const { return construct_point(Standard_line_2(-1,-1,0)); }
  /*{\Mop returns the point on the northwest frame corner.}*/

  Point_2 SW() const { return construct_point(Standard_line_2( 1,-1,0)); }
  /*{\Mop returns the point on the southwest frame corner.}*/


  Line_2 upper() const { return construct_line(NW(),NE()); }
  /*{\Xop returns the line underlying the upper frame segment.}*/

  Line_2 lower() const { return construct_line(SW(),SE()); }
  /*{\Xop returns the line underlying the lower frame segment.}*/

  Line_2 left()  const { return construct_line(SW(),NW()); }
  /*{\Xop returns the line underlying the left frame segment.}*/

  Line_2 right() const { return construct_line(SE(),NE()); }
  /*{\Xop returns the line underlying the right frame segment.}*/


  /*{\Mtext \headerline{Geometric kernel calls}}*/

  Point_2 source(const Segment_2& s) const
  /*{\Mop returns the source point of |s|.}*/
  { typename Base::Construct_vertex_2 _source = this->construct_vertex_2_object();
    return _source(s,0); }

  Point_2 target(const Segment_2& s) const
  /*{\Mop returns the target point of |s|.}*/
  { typename Base::Construct_vertex_2 _target = this->construct_vertex_2_object();
    return _target(s,1); }

  Segment_2 construct_segment(const Point_2& p, const Point_2& q) const
  /*{\Mop constructs a segment |pq|.}*/
  { typename Base::Construct_segment_2 _segment =
      this->construct_segment_2_object();
    return _segment(p,q); }

  void simplify(Point_2& p) const
  /*{\Xop only used internally.}*/
  { CGAL_NEF_TRACEN("simplify("<<p<<")");
    RT x=p.hx(), y=p.hy(), w=p.hw();
    RT common = x.is_zero() ? y : (RT) RT::gcd(x,y);
    common = RT::gcd(common,w);
    p = Point_2(x/common,y/common,w/common);
    CGAL_NEF_TRACEN("canceled="<<p);
  }

  Line_2 construct_line(const Standard_line_2& l)  const
  /*{\Xop only used internally.}*/
  { return Line_2(l.a(),l.b(),l.c()); }

  Line_2 construct_line(const Point_2& p1, const Point_2& p2) const
  /*{\Xop only used internally.}*/
  { Line_2 l(p1,p2);
      CGAL_NEF_TRACEN("eline("<<p1<<p2<<")="<<l);
    RT a=l.a(), b=l.b(), c=l.c();
    RT common = a.is_zero() ? b : (RT) RT::gcd(a,b);
    common = RT::gcd(common,c);
    l =  Line_2(a/common,b/common,c/common);
      CGAL_NEF_TRACEN("canceled="<<l);
    return l; 
  }

  int orientation(const Segment_2& s, const Point_2& p) const
  /*{\Mop returns the orientation of |p| with respect to the line
  through |s|.}*/
  { typename Base::Orientation_2 _orientation =
      this->orientation_2_object();
    return static_cast<int> ( _orientation(source(s),target(s),p) ); 
  }

  int orientation(const Point_2& p1, const Point_2& p2, const Point_2& p3) 
  const
  /*{\Mop returns the orientation of |p3| with respect to the line
  through |p1p2|.}*/
  { typename Base::Orientation_2 _orientation =
      this->orientation_2_object();
    return static_cast<int> ( _orientation(p1,p2,p3) ); 
  }

  bool left_turn(const Point_2& p1, const Point_2& p2, const Point_2& p3) 
  const
  /*{\Mop return true iff the |p3| is left of the line through |p1p2|.}*/
  { return orientation(p1,p2,p3) > 0; }
   
  bool is_degenerate(const Segment_2& s) const
  /*{\Mop return true iff |s| is degenerate.}*/
  { typename Base::Is_degenerate_2 _is_degenerate =
      this->is_degenerate_2_object();
    return _is_degenerate(s); }

  int compare_xy(const Point_2& p1, const Point_2& p2) const
  /*{\Mop returns the lexicographic order of |p1| and |p2|.}*/
  { typename Base::Compare_xy_2 _compare_xy =
      this->compare_xy_2_object();
    return static_cast<int>( _compare_xy(p1,p2) );
  }

  int compare_x(const Point_2& p1, const Point_2& p2) const
  /*{\Mop returns the order on the $x$-coordinates of |p1| and |p2|.}*/
  { typename Base::Compare_x_2 _compare_x =
      this->compare_x_2_object();
    return static_cast<int>( _compare_x(p1,p2) );
  }

  int compare_y(const Point_2& p1, const Point_2& p2) const
  /*{\Mop returns the order on the $y$-coordinates of |p1| and |p2|.}*/
  { typename Base::Compare_y_2 _compare_y =
      this->compare_y_2_object();
    return static_cast<int>( _compare_y(p1,p2) );
  }

  Point_2 intersection(
    const Segment_2& s1, const Segment_2& s2) const
  /*{\Mop returns the point of intersection of the lines supported by 
  |s1| and |s2|. \precond the intersection point exists.}*/
  { typename Base::Intersect_2 _intersect =
      this->intersect_2_object();
    typename Base::Construct_line_2 _line =
      this->construct_line_2_object();
    Point_2 p; 
    CGAL::Object result =
      _intersect(_line(s1),_line(s2));
    if ( !CGAL::assign(p, result) )
    CGAL_error_msg("intersection: no intersection.");
    simplify(p);
    return p;
  }

  Direction_2 construct_direction(
    const Point_2& p1, const Point_2& p2) const
  /*{\Mop returns the direction of the vector |p2| - |p1|.}*/
  { typename Base::Construct_direction_2 _direction =
      this->construct_direction_2_object();
    return _direction(construct_line(p1,p2)); }

  bool strictly_ordered_ccw(const Direction_2& d1, 
    const Direction_2& d2, const Direction_2& d3) const
  /*{\Mop returns |true| iff |d2| is in the interior of the
  counterclockwise angular sector between |d1| and |d3|.}*/
  { 
    if ( d1 < d2 )  return ( d2 < d3 )||( d3 <= d1 );
    if ( d1 > d2 )  return ( d2 < d3 )&&( d3 <= d1 );
    return false;
  }

  bool strictly_ordered_along_line(
    const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  /*{\Mop returns |true| iff |p2| is in the relative interior of the
  segment |p1p3|.}*/
  { typename Base::Are_strictly_ordered_along_line_2 _ordered =
      this->are_strictly_ordered_along_line_2_object();
    return _ordered(p1,p2,p3);
  }

  bool contains(const Segment_2& s, const Point_2& p) const
  /*{\Mop returns true iff |s| contains |p|.}*/
  { typename Base::Has_on_2 _contains = this->has_on_2_object();
    return _contains(s,p);
  }

  bool first_pair_closer_than_second(
    const Point_2& p1, const Point_2& p2, 
    const Point_2& p3, const Point_2& p4) const
  /*{\Mop returns true iff $\Labs{p1-p2} < \Labs{p3-p4}$.}*/
  { return ( squared_distance(p1,p2) < squared_distance(p3,p4) ); }

  void scale_first_by_second(RT& r1, RT& r2, RT& w) const
  { CGAL_assertion(w.degree()==0&&w!=RT(0)&& r2[1]!=Standard_RT(0));
    Standard_RT w_res = w[0]*r2[1];
    int sm2 = CGAL_NTS sign(r2[1]);
    RT r2_res = RT(Standard_RT(0),sm2 * w_res); 
    RT r1_res = RT(r2[1]*r1[0]-r1[1]*r2[0], w[0]*r1[1]*sm2);
    r1 = r1_res; r2 = r2_res; w = w_res;
  }

  Point_2 transform(const Point_2& p, 
                    const Standard_aff_transformation_2& t) const
  {
    RT tpx = t.homogeneous(0,0)*p.hx() + t.homogeneous(0,1)*p.hy() +
      t.homogeneous(0,2)*p.hw();
    RT tpy = t.homogeneous(1,0)*p.hx() + t.homogeneous(1,1)*p.hy() +
      t.homogeneous(1,2)*p.hw();
    RT tpw = t.homogeneous(2,2)*p.hw();
    if ( is_standard(p) ) {
      Point_2 res(tpx,tpy,tpw); simplify(res);
      return res;
    } 
    RT tpxa = CGAL_NTS abs(tpx);
    RT tpya = CGAL_NTS abs(tpy);
    if ( tpxa > tpya ) {
      scale_first_by_second(tpy,tpx,tpw);
    } else { // tpxa <= tpya
      scale_first_by_second(tpx,tpy,tpw);
    }
    Point_2 res(tpx,tpy,tpw); simplify(res);
    return res;
  }

  const char* output_identifier() const { return "Extended_homogeneous"; }
  /*{\Mop returns a unique identifier for kernel object input/output.}*/



};

} //namespace CGAL
#endif // CGAL_EXTENDED_HOMOGENEOUS_H

// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_BOUNDED_KERNEL_H
#define CGAL_BOUNDED_KERNEL_H

#include <CGAL/license/Nef_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Intersections_2/Line_2_Line_2.h>


#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 51
#include <CGAL/Nef_2/debug.h>
#include <CGAL/Nef_2/Line_to_epoint.h>

namespace CGAL {

template <class T> class Bounded_kernel;

template<class Kernel>
struct Is_extended_kernel;

template<class T>
struct Is_extended_kernel<Bounded_kernel<T> > {
       typedef Tag_false value_type;
};

/*{\Xanpage {Bounded_kernel}{}{An extended geometric kernel model}{K}}*/

template <class T>
class Bounded_kernel
  : public T
{

public:
  typedef T Base;
  typedef Bounded_kernel<T> Self;
  typedef T Standard_kernel;


  typedef typename Standard_kernel::RT Standard_RT;
  /*{\Xtypemember the standard ring type.}*/

  typedef typename Standard_kernel::FT Standard_FT;
  /*{\Xtypemember the field type.}*/

  typedef typename Standard_kernel::Point_2     Standard_point_2;
  /*{\Xtypemember standard points.}*/

  typedef typename Standard_kernel::Segment_2   Standard_segment_2;
  /*{\Xtypemember standard segments.}*/

  typedef typename Standard_kernel::Line_2      Standard_line_2;
  /*{\Xtypemember standard oriented lines.}*/

  typedef typename Standard_kernel::Direction_2 Standard_direction_2;
  /*{\Xtypemember standard directions.}*/

  typedef typename Standard_kernel::Ray_2       Standard_ray_2;
  /*{\Xtypemember standard rays.}*/

  typedef typename Standard_kernel::Aff_transformation_2
  Standard_aff_transformation_2;
  /*{\Xtypemember standard affine transformations.}*/

  /*{\Xtext \headerline{Extended kernel types}}*/

  typedef typename Base::RT RT;
  /*{\Xtypemember the ring type of our extended kernel.}*/

  typedef typename Base::FT FT;
  /*{\Xtypemember the ring type of our extended kernel.}*/

  typedef typename Base::Point_2      Point_2;
  /*{\Xtypemember extended points.}*/

  typedef typename Base::Segment_2    Segment_2;
  /*{\Xtypemember extended segments.}*/

  typedef typename Base::Line_2       Line_2;
  /*{\Xtypemember extended lines.}*/

  typedef typename Base::Direction_2  Direction_2;
  /*{\Xtypemember extended directions.}*/

  enum Point_type { SWCORNER=1, LEFTFRAME, NWCORNER,
                    BOTTOMFRAME, STANDARD, TOPFRAME,
                    SECORNER, RIGHTFRAME, NECORNER };
  /*{\Xenum a type descriptor for extended points.}*/

  Point_2 epoint(const Standard_FT& /*m1*/, const Standard_FT& /*n1*/,
                 const Standard_FT& /*m2*/, const Standard_FT& /*n2*/) const
  {
    CGAL_error_msg( "Bounded_kernel::epoint(..) should not be called");
    return Point_2();
  }

public:

  Point_2
  construct_point(const Standard_point_2& p) const
  {
    return p;
  }

  Point_2
  construct_point(const Standard_line_2& , Point_type& ) const
  {
    CGAL_error_msg( "Bounded_kernel::construct_point(Line,Point_type) should not be called");
    return Point_2();
  }

  Point_2
  construct_point(const Standard_point_2& ,
                  const Standard_point_2& ,
                  Point_type& /*t*/) const
  {
    CGAL_error_msg( "Bounded_kernel::construct_point(Point,Point) should not be called");
    return Point_2();
  }

  Point_2
  construct_point(const Standard_line_2& ) const
  {
    CGAL_error_msg( "Bounded_kernel::construct_point(Line) should not be called");
    return Point_2();
  }

  Point_2
  construct_point(const Standard_point_2& ,
                  const Standard_point_2& ) const
  {
    CGAL_error_msg( "Bounded_kernel::construct_point(Point,Point) should not be called");
    return Point_2();
   }

  Point_2 construct_point(const Standard_point_2& ,
                          const Standard_direction_2& ) const
  {
    CGAL_error_msg( "Bounded_kernel::construct_point(Point,Direction) should not be called");
    return Point_2();
  }

  Point_2
  construct_opposite_point(const Standard_line_2& /*l*/) const
  {
    CGAL_error_msg( "Bounded_kernel::construct_opposite_point(..) should not be called");
    return Point_2();
  }

  Point_type
  type(const Point_2& /*p*/) const
  {
    return STANDARD;
  }


  bool
  is_standard(const Point_2& /*p*/) const
  {
    return true;
  }

  Standard_point_2
  standard_point(const Point_2& p) const
  {
    return p;
  }

  Standard_line_2
  standard_line(const Point_2& /*p*/) const
  {
    CGAL_error_msg( "Bounded_kernel::standard_line(..) should not be called");
    return Standard_line_2();
  }

  Standard_ray_2
  standard_ray(const Point_2& /*p*/) const
  {
    CGAL_error_msg( "Bounded_kernel::standard_ray(..) should not be called");
    return Standard_ray_2();
  }

  Point_2
  NE() const
  {
    CGAL_error_msg( "Bounded_kernel::NE(..) should not be called");
    return Point_2();
  }


  Point_2
  SE() const
  {
    CGAL_error_msg( "Bounded_kernel::SE(..) should not be called");
    return Point_2();
  }


  Point_2
  NW() const
  {
    CGAL_error_msg( "Bounded_kernel::NW(..) should not be called");
    return Point_2();
  }


  Point_2
  SW() const
  {
    CGAL_error_msg( "Bounded_kernel::SW(..) should not be called");
    return Point_2();
  }


  Line_2
  upper() const
  {
    CGAL_error_msg( "Bounded_kernel::upper(..) should not be called");
    return Line_2();
  }


  Line_2
  lower() const
  {
    CGAL_error_msg( "Bounded_kernel::lower(..) should not be called");
    return Line_2();
  }


  Line_2
  left()  const
  {
    CGAL_error_msg( "Bounded_kernel::left(..) should not be called");
    return Line_2();
  }



  Line_2
  right() const
  {
    CGAL_error_msg( "Bounded_kernel::right(..) should not be called");
    return Line_2();
  }



  Point_2
  source(const Segment_2& s) const
  {
    typename Base::Construct_vertex_2 _source =
      this->construct_vertex_2_object();
    return _source(s,0); }

  Point_2
  target(const Segment_2& s) const
  {
    typename Base::Construct_vertex_2 _target =
      this->construct_vertex_2_object();
    return _target(s,1); }

  Segment_2
  construct_segment(const Point_2& p, const Point_2& q) const
  {
    typename Base::Construct_segment_2 _segment =
      this->construct_segment_2_object();
    return _segment(p,q);
  }

  Line_2
  construct_line(const Standard_line_2& l)  const
  {
    return Line_2(l.a(),l.b(),l.c());
  }

  Line_2
  construct_line(const Point_2& p1, const Point_2& p2) const
  {
    return Line_2(p1,p2);
  }


  int
  orientation(const Segment_2& s, const Point_2& p) const
  {
    typename Base::Orientation_2 _orientation =
      this->orientation_2_object();
    return static_cast<int> ( _orientation(source(s),target(s),p) );
  }

  int
  orientation(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
    typename Base::Orientation_2 _orientation =
      this->orientation_2_object();
    return static_cast<int> ( _orientation(p1,p2,p3) );
  }

  bool
  left_turn(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
    return orientation(p1,p2,p3) > 0;
  }

  bool
  is_degenerate(const Segment_2& s) const
  {
    typename Base::Is_degenerate_2 _is_degenerate =
      this->is_degenerate_2_object();
    return _is_degenerate(s);
  }

  int
  compare_xy(const Point_2& p1, const Point_2& p2) const
  {
    typename Base::Compare_xy_2 _compare_xy =
      this->compare_xy_2_object();
    return static_cast<int>( _compare_xy(p1,p2) );
  }

  int
  compare_x(const Point_2& p1, const Point_2& p2) const
  {
    typename Base::Compare_x_2 _compare_x =
      this->compare_x_2_object();
    return static_cast<int>( _compare_x(p1,p2) );
  }

  int
  compare_y(const Point_2& p1, const Point_2& p2) const
  {
    typename Base::Compare_y_2 _compare_y =
      this->compare_y_2_object();
    return static_cast<int>( _compare_y(p1,p2) );
  }


  Point_2
  intersection(const Segment_2& s1, const Segment_2& s2) const
  {
    typename Base::Intersect_2 _intersect =
      this->intersect_2_object();
    typename Base::Construct_line_2 _line =
      this->construct_line_2_object();
    Point_2 p;
    Line_2 l1 = _line(s1);
    Line_2 l2 = _line(s2);

    CGAL::Object result =
      _intersect(l1, l2);
    if ( !CGAL::assign(p, result) )
      CGAL_error_msg("intersection: no intersection.");
    return p;
  }

  Direction_2
  construct_direction(const Point_2& p1, const Point_2& p2) const
  {
    typename Base::Construct_direction_2 _direction =
      this->construct_direction_2_object();
    return _direction(construct_line(p1,p2));
  }

  bool
  strictly_ordered_ccw(const Direction_2& d1,
                       const Direction_2& d2, const Direction_2& d3) const
  {
    if ( d1 < d2 )  return ( d2 < d3 )||( d3 <= d1 );
    if ( d1 > d2 )  return ( d2 < d3 )&&( d3 <= d1 );
    return false;
  }

  bool
  contains(const Segment_2& s, const Point_2& p) const
  {
    typename Base::Has_on_2 _contains = this->has_on_2_object();
    return _contains(s,p);
  }

  bool
  strictly_ordered_along_line(const Point_2& p1, const Point_2& p2, const Point_2& p3) const
  {
    typename Base::Are_strictly_ordered_along_line_2 _ordered =
      this->are_strictly_ordered_along_line_2_object();
    return _ordered(p1,p2,p3);
  }

  bool
  first_pair_closer_than_second(const Point_2& p1, const Point_2& p2,
                                const Point_2& p3, const Point_2& p4) const
  {
    return ( squared_distance(p1,p2) < squared_distance(p3,p4) );
  }

  template <class Forward_iterator>
  void
  determine_frame_radius(Forward_iterator start, Forward_iterator end,
                         Standard_RT& R0) const
  {
    Standard_RT R;
    while ( start != end ) {
      Point_2 p = *start++;
      if ( is_standard(p) ) {
        R = (CGAL::max)(CGAL_NTS abs(p.x()[0]), CGAL_NTS abs(p.y()[0]));
      } else {
        RT rx = CGAL_NTS abs(p.x()), ry = CGAL_NTS abs(p.y());
        if ( rx[1] > ry[1] )      R = CGAL_NTS abs(ry[0]-rx[0])/(rx[1]-ry[1]);
        else if ( rx[1] < ry[1] ) R = CGAL_NTS abs(rx[0]-ry[0])/(ry[1]-rx[1]);
        else /* rx[1] == ry[1] */ R = CGAL_NTS abs(rx[0]-ry[0])/2;
      }
      R0 = (CGAL::max)(R+1,R0);
    }
  }



  const char*
  output_identifier() const
  {
    return "Bounded_kernel";
  }

};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_BOUNDED_KERNEL_H

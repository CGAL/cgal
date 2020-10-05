// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado,
//             Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a
// STREP (FET Open) Project under Contract No  IST-006413
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_LINE_ARC_3_H
#define CGAL_LINE_ARC_3_H

#include <CGAL/license/Circular_kernel_3.h>


#include <CGAL/result_of.h>
#include <CGAL/Bbox_3.h>


namespace CGAL {
  template <class SK>
    class Line_arc_3
    : public SK::Kernel_base::Line_arc_3
  {

    typedef typename SK::RT                    RT;
    typedef typename SK::FT                    FT;
    typedef typename SK::Line_3                Line_3;
    typedef typename SK::Point_3               Point_3;
    typedef typename SK::Plane_3               Plane_3;
    typedef typename SK::Sphere_3              Sphere_3;
    typedef typename SK::Segment_3             Segment_3;
    typedef typename SK::Circular_arc_point_3  Circular_arc_point_3;
    typedef typename SK::Kernel_base::Line_arc_3 RLine_arc_3;


  public:
    typedef  RLine_arc_3 Rep;
    typedef  SK   R;

    const Rep& rep() const
      {
        return *this;
      }

    Rep& rep()
      {
        return *this;
      }

    Line_arc_3()
      : RLine_arc_3(typename R::Construct_line_arc_3()())
      {}

    Line_arc_3(const Line_3& l,
               const Circular_arc_point_3& s,
               const Circular_arc_point_3& t)
      : RLine_arc_3(typename R::Construct_line_arc_3()(l,s,t))
      {}

    Line_arc_3(const Point_3& s,
               const Point_3& t)
      : RLine_arc_3(typename R::Construct_line_arc_3()(s,t))
      {}

    Line_arc_3(const Segment_3 &s)
      : RLine_arc_3(typename R::Construct_line_arc_3()(s))
      {}

    // Not Documented
    Line_arc_3(const Line_3 &l,
               const Sphere_3 &s,
               bool less_xyz_first = true)
      : RLine_arc_3(typename R::Construct_line_arc_3()(l,s,less_xyz_first))
    {}

    // Not Documented
    Line_arc_3(const Sphere_3 &s,
               const Line_3 &l,
               bool less_xyz_first = true)
      : RLine_arc_3(typename R::Construct_line_arc_3()(s,l,less_xyz_first))
    {}

    // Not Documented
    Line_arc_3(const Line_3 &l,
               const Sphere_3 &s1, bool less_xyz_s1,
               const Sphere_3 &s2, bool less_xyz_s2)
      : RLine_arc_3(typename R::Construct_line_arc_3()(l,s1,less_xyz_s1,
                                                         s2,less_xyz_s2))
    {}

    // Not Documented
    Line_arc_3(const Sphere_3 &s1, bool less_xyz_s1,
               const Sphere_3 &s2, bool less_xyz_s2,
               const Line_3 &l)
      : RLine_arc_3(typename R::Construct_line_arc_3()(s1,less_xyz_s1,
                                                       s2,less_xyz_s2,l))
    {}

    // Not Documented
    Line_arc_3(const Line_3 &l,
               const Plane_3 &p1,
               const Plane_3 &p2)
      : RLine_arc_3(typename R::Construct_line_arc_3()(l,p1,p2))
    {}

    // Not Documented
    Line_arc_3(const Plane_3 &p1,
               const Plane_3 &p2,
               const Line_3 &l)
      : RLine_arc_3(typename R::Construct_line_arc_3()(p1,p2,l))
    {}

    Line_arc_3(const RLine_arc_3 &a)
     : RLine_arc_3(a)
      {}

    typename cpp11::result_of<typename R::Construct_circular_source_vertex_3(Line_arc_3)>::type
    source() const
    {
      return typename R::Construct_circular_source_vertex_3()(*this);
    }

    typename cpp11::result_of<typename R::Construct_circular_target_vertex_3(Line_arc_3)>::type
    target() const
    {
      return typename R::Construct_circular_target_vertex_3()(*this);
    }

    typename cpp11::result_of<typename R::Construct_circular_min_vertex_3(Line_arc_3)>::type
    lower_xyz_extremity() const
    {
      return typename R::Construct_circular_min_vertex_3()(*this);
    }

    typename cpp11::result_of<typename R::Construct_circular_max_vertex_3(Line_arc_3)>::type
    higher_xyz_extremity() const
    {
      return typename R::Construct_circular_max_vertex_3()(*this);
    }

    typename cpp11::result_of<typename R::Construct_line_3(Line_arc_3)>::type
    supporting_line() const
    {
      return typename R::Construct_line_3()(*this);
    }

    Bbox_3 bbox() const
    { return typename R::Construct_bbox_3()(*this); }

  };

  template < typename SK >
  inline
  bool
  operator==(const Line_arc_3<SK> &p,
             const Line_arc_3<SK> &q)
  {
    return SK().equal_3_object()(p, q);
  }

  template < typename SK >
  inline
  bool
  operator!=(const Line_arc_3<SK> &p,
             const Line_arc_3<SK> &q)
  {
    return ! (p == q);
  }

  template < typename SK >
  std::ostream &
  operator<<(std::ostream & os, const Line_arc_3<SK> &a)
  {

    return os << a.supporting_line() << " "
              << a.source() << " "
              << a.target() << " ";
  }

  template < typename SK >
  std::istream &
  operator>>(std::istream & is, Line_arc_3<SK> &a)
  {
    typename SK::Line_3 l;
    typename SK::Circular_arc_point_3 p1;
    typename SK::Circular_arc_point_3 p2;
    is >> l >> p1 >> p2 ;
    if (is)
      a = Line_arc_3<SK>(l, p1, p2);
    return is;
  }
}

#endif

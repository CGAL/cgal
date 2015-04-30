// Copyright (c) 2008  INRIA Sophia-Antipolis (France).
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
// Author(s) : Monique Teillaud, Sylvain Pion, Pedro Machado, 
//             Sebastien Loriot, Julien Hazebrouck, Damien Leroy

// Partially supported by the IST Programme of the EU as a 
// STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)

#ifndef CGAL_SPHERICAL_KERNEL_LINE_ARC_3_H
#define CGAL_SPHERICAL_KERNEL_LINE_ARC_3_H

#include <CGAL/Circular_kernel_3/internal_functions_on_sphere_3.h>
#include <CGAL/Circular_kernel_3/Intersection_traits.h>
#include <boost/tuple/tuple.hpp>

namespace CGAL {
  namespace internal{
    template <class SK> class Line_arc_3 {

      typedef typename SK::Plane_3              Plane_3;
      typedef typename SK::Sphere_3             Sphere_3;
      typedef typename SK::Point_3              Point_3;
      typedef typename SK::Segment_3            Segment_3;
      typedef typename SK::Circular_arc_point_3 Circular_arc_point_3;
      typedef typename SK::Line_3               Line_3;
      typedef typename SK::FT                   FT;

    private:
      typedef boost::tuple<Line_3, Circular_arc_point_3, 
                             Circular_arc_point_3>  Rep;
      typedef typename SK::template Handle<Rep>::type  Base;

      Base base;
      mutable unsigned char begin_less_xyz_than_end_flag;

      bool begin_less_xyz_than_end() const {
        if(begin_less_xyz_than_end_flag == 0) {
          if(SK().compare_xyz_3_object()(source(), target()) < 0)
            begin_less_xyz_than_end_flag = 2;
          else begin_less_xyz_than_end_flag = 1;
        } return begin_less_xyz_than_end_flag == 2;
      }

    public:
      Line_arc_3()
      : begin_less_xyz_than_end_flag(0) 
      {}

      Line_arc_3(const Line_3 &l, 
                 const Circular_arc_point_3 &s,
                 const Circular_arc_point_3 &t) 
      : begin_less_xyz_than_end_flag(0)
      {
        // l must pass through s and t, and s != t
        CGAL_kernel_precondition(SK().has_on_3_object()(l,s));
        CGAL_kernel_precondition(SK().has_on_3_object()(l,t));
        CGAL_kernel_precondition(s != t);
        base = Rep(l,s,t);
      }

      Line_arc_3(const Segment_3 &s) 
      : begin_less_xyz_than_end_flag(0)
      {
        base = Rep(s.supporting_line(),
                   s.source(),
                   s.target());
      }

      Line_arc_3(const Point_3 &s,
                 const Point_3 &t) 
      : begin_less_xyz_than_end_flag(0)
      {
        CGAL_kernel_precondition(s != t);
        base = Rep(SK().construct_line_3_object()(s,t),s,t);
      }

      Line_arc_3(const Line_3 &l, 
                 const Sphere_3 &s,
                 bool less_xyz_first = true) 
      {
        std::vector<typename SK3_Intersection_traits<SK, Line_3, Sphere_3>::type> sols;
         SK().intersect_3_object()(l, s, std::back_inserter(sols));
         // l must intersect s in 2 points 
         CGAL_kernel_precondition(sols.size() == 2);
         const std::pair<typename SK::Circular_arc_point_3, unsigned>& pair1=
           *boost::get<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[0]);
         const std::pair<typename SK::Circular_arc_point_3, unsigned>& pair2=
            *boost::get<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols[1]);
         if(less_xyz_first) {
           *this = Line_arc_3(l, pair1.first, pair2.first);
         } else {
           *this = Line_arc_3(l, pair2.first, pair1.first);
         } 
      }

      Line_arc_3(const Line_3 &l, 
                 const Sphere_3 &s1, bool less_xyz_s1,
                 const Sphere_3 &s2, bool less_xyz_s2) 
      {
        std::vector<typename SK3_Intersection_traits<SK, Line_3, Sphere_3>::type> sols1, sols2;
         SK().intersect_3_object()(l, s1, std::back_inserter(sols1));
         SK().intersect_3_object()(l, s2, std::back_inserter(sols2));
         // l must intersect s1 and s2
         CGAL_kernel_precondition(sols1.size() > 0);
         CGAL_kernel_precondition(sols2.size() > 0);
         const std::pair<typename SK::Circular_arc_point_3, unsigned>& pair1=
            *boost::get<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols1[(sols1.size()==1)?(0):(less_xyz_s1?0:1)]);
         const std::pair<typename SK::Circular_arc_point_3, unsigned>& pair2=
            *boost::get<std::pair<typename SK::Circular_arc_point_3, unsigned> >(&sols2[(sols2.size()==1)?(0):(less_xyz_s2?0:1)]);
         // the source and target must be different
         CGAL_kernel_precondition(pair1.first != pair2.first);
         *this = Line_arc_3(l, pair1.first, pair2.first);
      }

      Line_arc_3(const Line_3 &l, 
                 const Plane_3 &p1,
                 const Plane_3 &p2) 
      {
         // l must not be on p1 or p2
         CGAL_kernel_precondition(!SK().has_on_3_object()(p1,l));
         CGAL_kernel_precondition(!SK().has_on_3_object()(p2,l));
         // l must intersect p1 and p2
         typedef typename SK3_Intersection_traits<SK, Line_3, Plane_3>::type Intersection;
         Intersection i1 = SK().intersect_3_object()(l, p1);
         Intersection i2 = SK().intersect_3_object()(l, p2);
         const typename SK::Point_3* point1=boost::get<typename SK::Point_3>( & *i1 );
         const typename SK::Point_3* point2=boost::get<typename SK::Point_3>( & *i2 );
         CGAL_assertion(point1!=NULL);
         CGAL_assertion(point2!=NULL);
         // the source and target must be different
         CGAL_kernel_precondition(*point1 != *point2);
         *this = Line_arc_3(l, *point1, *point2);
      }

      const Line_3& supporting_line() const 
      {
        return get_pointee_or_identity(base).template get<0>();
      }

      const Circular_arc_point_3& source() const 
      {
        return get_pointee_or_identity(base).template get<1>();
      }

      const Circular_arc_point_3& target() const 
      {
        return get_pointee_or_identity(base).template get<2>();
      }

      const Circular_arc_point_3& lower_xyz_extremity() const
      {
        return begin_less_xyz_than_end() ? source() : target();
      }

      const Circular_arc_point_3& higher_xyz_extremity() const
      {
        return begin_less_xyz_than_end() ? target() : source();
      }

      const CGAL::Bbox_3 bbox() const {
        return source().bbox() + target().bbox();
      }

      bool operator==(const Line_arc_3 &) const;
      bool operator!=(const Line_arc_3 &) const;

    };

    template < class SK >
    CGAL_KERNEL_INLINE
    bool
    Line_arc_3<SK>::operator==(const Line_arc_3<SK> &t) const
    {
      if (CGAL::identical(base, t.base))
        return true;
      return CGAL::SphericalFunctors::non_oriented_equal<SK>(*this, t);
    }

    template < class SK >
    CGAL_KERNEL_INLINE
    bool
    Line_arc_3<SK>::operator!=(const Line_arc_3<SK> &t) const
    {
      return !(*this == t);
    }

  }
}

#endif

// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Envelope_3/include/CGAL/Env_plane_traits_3.h $
// $Id: Env_plane_traits_3.h 51989 2009-09-21 10:55:53Z efif $
//
// Author(s)     : Ophir Setter

#ifndef CGAL_L1_VORONOI_DIAGRAM_TRAITS_2_H
#define CGAL_L1_VORONOI_DIAGRAM_TRAITS_2_H

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/representation_tags.h>
#include <CGAL/enum.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/number_utils.h> 
#include <CGAL/Envelope_3/Envelope_base.h>
#include <CGAL/Envelope_3/Env_plane_traits_3_functions.h>

namespace CGAL{

template <class Kernel_>
class L1_voronoi_traits_2 : public Arr_linear_traits_2<Kernel_>
{
public:
  typedef Kernel_                              Kernel;
  typedef typename Kernel::FT                  FT;
  typedef Arr_linear_traits_2<Kernel>          Base;
  typedef L1_voronoi_traits_2<Kernel>        Self;
  typedef std::size_t                          Multiplicity;

  typedef typename Base::Point_2               Point_2;
  typedef typename Base::Curve_2               Curve_2;
  typedef typename Base::X_monotone_curve_2    X_monotone_curve_2;
  typedef typename Kernel::Segment_2           Segment_2;
  typedef typename Kernel::Ray_2               Ray_2;
  typedef typename Kernel::Line_2              Line_2;
  typedef typename Kernel::Direction_2         Direction_2;
  typedef std::pair<Curve_2, Multiplicity>     Intersection_curve;

  typedef typename Base::Left_side_category    Left_side_category;
  typedef typename Base::Bottom_side_category  Bottom_side_category;
  typedef typename Base::Top_side_category     Top_side_category;
  typedef typename Base::Right_side_category   Right_side_category;
  
  typedef Point_2                  Xy_monotone_surface_3;
  typedef Point_2                  Surface_3;

  // Returns the distance between a two points in L1 metric.
  static FT distance(const Point_2& p1, const Point_2& p2) {
    FT dist = CGAL::abs(p1.x() - p2.x()) + CGAL::abs(p1.y() - p2.y());
    return dist;
  }

  // Returns the midpoint (under the L1 metric) that is on the rectangle
  // defined by the two points (the rectangle can be degenerate).
  // As there are to enpoints, the index determines which is returned
  static Point_2 midpoint(const Point_2& p1, const Point_2& p2, std::size_t index) {
    const Point_2 *pp1;
    const Point_2 *pp2;
    
    if (index % 2 == 0) {
      pp1 = &p1;
      pp2 = &p2;
    } else {
      pp1 = &p2;
      pp2 = &p1;
    }


    FT delta_x = pp2->x() - pp1->x();
    FT delta_y = pp2->y() - pp1->y();
    
    FT sign_x = CGAL::sign(delta_x);
    FT sign_y = CGAL::sign(delta_y);

    FT abs_x = CGAL::abs(delta_x);
    FT abs_y = CGAL::abs(delta_y);

    FT dist = (abs_x + abs_y) / 2;
    FT mid_x = pp1->x();
    FT mid_y = pp1->y();
    
    // Walk on the horizontal edge of the rectangle and then on the vertical.

    // There is a chance that the width of the rectangle is smaller then the mid-dist.
    FT walk_x = (CGAL::min)(abs_x, dist);
    mid_x += sign_x * walk_x;
    dist -= walk_x;
    
    CGAL_assertion(abs_y > dist);
    mid_y += sign_y * dist;
    
    return Point_2(mid_x, mid_y);
  }

  static Comparison_result compare_z_at_xy (const X_monotone_curve_2& cv,
                                            const Xy_monotone_surface_3& h1,
                                            const Xy_monotone_surface_3& h2,
                                            bool above) {
    CGAL::Comparison_result side = (above == true) ? CGAL::LARGER : CGAL::SMALLER;

    Line_2 l;
    if (cv.is_segment())
      l = cv.segment().supporting_line();
    else if (cv.is_ray())
      l = cv.ray().supporting_line();
    else
      l = cv.line();

    if (l.is_vertical()) {
      // Could be a tie.
      // To be "above" the curve, we acutually need to have smaller x coordinate,
      // the order of the comparison function here is opposite to the none vertical
      // case.
      side = CGAL::opposite(side);
      CGAL::Comparison_result res = CGAL::compare_x_at_y(h1, l);
      if (res == side)
        return CGAL::SMALLER;

      res = CGAL::compare_x_at_y(h2, l);
      if (res == side)
        return CGAL::LARGER;

      return CGAL::EQUAL;
    } else {
      // One of the points in indeed closer (tie can only happen on vertical lines).
      CGAL::Comparison_result res = CGAL::compare_y_at_x(h1, l);
      if (l.is_horizontal()) {
        CGAL_assertion(CGAL::compare_y_at_x(h2, l) != res);
        if (res == side)
          return CGAL::SMALLER;

        res = CGAL::compare_y_at_x(h2, l);
        if (res == side)
          return CGAL::LARGER;
        
        return CGAL::EQUAL;
      }
      
      CGAL_assertion(CGAL::compare_y_at_x(h2, l) == CGAL::opposite(res));
      if (res == side) 
        return CGAL::SMALLER;
      else
        return CGAL::LARGER;
    }
  }
  class Make_xy_monotone_3 {
  public:
    template <class OutputIterator>
      OutputIterator operator()(const Surface_3& s,
                                bool /* is_lower */,
                                OutputIterator o) const {
      *o++ = s;
      return o;
    }
  };
  
  Make_xy_monotone_3 make_xy_monotone_3_object() const {
    return Make_xy_monotone_3();
  }
  
  class Compare_z_at_xy_3 {
  public:
    Comparison_result operator()(const Point_2& p,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return CGAL::compare(distance(p, h1), distance(p, h2));
    }

    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      Kernel k;
      Point_2 p;
      if(cv.is_segment())
        p = k.construct_midpoint_2_object()(cv.left(), cv.right());
      else
        if(cv.is_ray())
          p = k.construct_point_on_2_object()(cv.ray(), 1);
        else {
          CGAL_assertion(cv.is_line());
          p = k.construct_point_on_2_object()(cv.line(), 1);
        }
      return this->operator()(p, h1, h2); 
    }
    
    Comparison_result operator()(const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      // should happen only if the points are equal.
      CGAL_assertion(h1 == h2);
      return EQUAL;
    }
  };

  Compare_z_at_xy_3 compare_z_at_xy_3_object() const
  {
    return Compare_z_at_xy_3();
  }

  class Compare_z_at_xy_above_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return compare_z_at_xy (cv, h1, h2, true);
    }

  };

  Compare_z_at_xy_above_3 compare_z_at_xy_above_3_object() const
  {
    return Compare_z_at_xy_above_3();
  }

  class Compare_z_at_xy_below_3
  {
  public:
    Comparison_result operator()(const X_monotone_curve_2& cv,
                                 const Xy_monotone_surface_3& h1,
                                 const Xy_monotone_surface_3& h2) const {
      return compare_z_at_xy (cv, h1, h2, false);
    }
  };

  Compare_z_at_xy_below_3 compare_z_at_xy_below_3_object() const {
    return Compare_z_at_xy_below_3();
  }
  
  class Construct_projected_boundary_2 {
  public:
    template <class OutputIterator>
    OutputIterator operator()(const Xy_monotone_surface_3& s,
                              OutputIterator o) const {
      return o;
    }
  };

  Construct_projected_boundary_2 
  construct_projected_boundary_2_object() const
  {
    return Construct_projected_boundary_2();
  }


  class Construct_projected_intersections_2 {
  public:
    template <class OutputIterator>
      OutputIterator operator()(const Xy_monotone_surface_3& s1,
                                const Xy_monotone_surface_3& s2,
                                OutputIterator o) const {
      // The bisector construction is based on the description of the bisector
      // from Handbook of Computational Geometry - the chapter on Voronoi diagrams
      // by F. Aurenhammer and R. Klein.

      FT delta_x = s2.x() - s1.x();
      FT delta_y = s2.y() - s1.y();

      CGAL::Sign s_x = CGAL::sign(delta_x);
      CGAL::Sign s_y = CGAL::sign(delta_y);

      if (s_x == CGAL::ZERO && s_y == CGAL::ZERO)
        return o;
      
      // The sites have the same x/y coordinate.
      if (s_x == CGAL::ZERO || s_y == CGAL::ZERO) {
        *o++ = CGAL::make_object(Intersection_curve(CGAL::bisector(s1, s2), 1));
        return o;
      }

      Point_2 p1 = midpoint(s1, s2, 0);
      Point_2 p2 = midpoint(s1, s2, 1);
      CGAL_assertion(p1 != p2);

      Point_2 *top, *bottom, *left, *right;
      if (CGAL::sign(p2.x() - p1.x()) == CGAL::POSITIVE) {
        left = &p1;
        right = &p2;
      }
      else {
        CGAL_assertion(CGAL::sign(p2.x() - p1.x()) == CGAL::NEGATIVE);
        left = &p2;
        right = &p1;
      }
      
      if (CGAL::sign(p2.y() - p1.y()) == CGAL::POSITIVE) {
        bottom = &p1;
        top = &p2;
      }
      else {
        CGAL_assertion(CGAL::sign(p2.y() - p1.y()) == CGAL::NEGATIVE);
        bottom = &p2;
        top = &p1;
      }
      
      // We construct the diagonal line either way.
      *o++ = CGAL::make_object(Intersection_curve(Segment_2(p1, p2), 1));

      // Now construct vertical rays. Two or four rays. If it is only two rays,
      // then the multiplicity of all the curves is 1.
      CGAL::Sign s_d = CGAL::sign(CGAL::abs(delta_y) - CGAL::abs(delta_x));
      std::size_t mult = (s_d == CGAL::ZERO) ? 0 : 1;

      if (s_d != CGAL::POSITIVE) {
        // horizontal rectangle or square = vertical rays.
        *o++ = CGAL::make_object(Intersection_curve(Ray_2(*top, Direction_2(0, 1)), mult));
        *o++ = CGAL::make_object(Intersection_curve(Ray_2(*bottom, Direction_2(0, -1)), mult));
      }

      if (s_d != CGAL::NEGATIVE) {
        // vertical rectangle or square = horizontal rays.
        *o++ = CGAL::make_object(Intersection_curve(Ray_2(*right, Direction_2(1, 0)), mult));
        *o++ = CGAL::make_object(Intersection_curve(Ray_2(*left, Direction_2(-1, 0)), mult));
      }
      
      return o;
    }
  };

  Construct_projected_intersections_2
  construct_projected_intersections_2_object() const {
    return Construct_projected_intersections_2();
  }

};

} //namespace CGAL

#endif // CGAL_L1_VORONOI_DIAGRAM_TRAITS_2_H

 // Copyright (c) 2009  GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Laurent Rineau


#ifndef CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H
#define CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H

#include <CGAL/Profile_timer.h>
#include <CGAL/intersections.h>

namespace CGAL {

namespace TriangulationProjectionTraitsCartesianFunctors {

template <class Traits>
class Projected_orientation_with_normal_3
{
  // private members
  const Traits& traits;

  // private type aliases
  typedef typename Traits::K K;
  typedef typename Traits::Point_2 Point;
public:
  typedef typename K::Orientation Orientation;
  typedef Orientation result_type;

  Projected_orientation_with_normal_3(const Traits& traits_)
    : traits(traits_)
  {
  }

  Orientation operator()(const Point& p,
			 const Point& q,
			 const Point& r) const
  {    
    CGAL_TIME_PROFILER("Projected_orientation_with_normal_3");
    return orientation(q-p, r-p, traits.normal());
  }
}; // end class Projected_orientation_with_normal_3<Traits>

template <class Traits>
class Projected_side_of_oriented_circle_with_normal_3
{
  // private members
  const Traits& traits;

  // private types aliases
  typedef typename Traits::K K;
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::FT FT;

public:
  typedef typename K::Oriented_side Oriented_side;
  typedef Oriented_side result_type;

  Projected_side_of_oriented_circle_with_normal_3(const Traits& traits_)
    : traits(traits_)
  {
  }

  Oriented_side operator()(const Point& p,
			   const Point& q,
			   const Point& r,
			   const Point& t) const
  {
    CGAL_TIME_PROFILER("Projected_side_of_oriented_circle_with_normal_3")
    const Vector_3& u = traits.normal();

    const Vector_3 tp = p - t;
    const Vector_3 tq = q - t;
    const Vector_3 tr = r - t;

    const FT tp2 = tp * tp;
    const FT tq2 = tq * tq;
    const FT tr2 = tr * tr;
    const FT u2  =  u * u;

    const FT k_p = tp * u;
    const FT k_q = tq * u;
    const FT k_r = tr * u;

    return 
       sign_of_determinant(tp.x(), tp.y(), tp.z(), (tp2 + k_p) * u2 - k_p * k_p,
			   tr.x(), tr.y(), tr.z(), (tr2 + k_r) * u2 - k_r * k_r,
			   tq.x(), tq.y(), tq.z(), (tq2 + k_q) * u2 - k_q * k_q,
                           u.x(),  u.y(),  u.z(), u2 * u2);
    // Note that q and r have been swapped in the determinant above, to
    // inverse its sign.
  }
}; // end class Projected_side_of_oriented_circle_with_normal_3

template <class Traits>
class Projected_square_distance_to_line_with_normal_3
{
  // private members
  const Traits& traits;
  
  // private types aliases
  typedef typename Traits::K K;
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Line_2 Line;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::FT FT;

public:
  Projected_square_distance_to_line_with_normal_3(const Traits& traits_)
    : traits(traits_)
  {}

  FT operator()(const Line& line, const Point& p)
  {
    const Vector_3& vl = line.to_vector();
    const Point& pl = line.point();
    const Vector_3 v = cross_product(traits.normal(),
                                     vl);
    if(v == NULL_VECTOR) {
      // den == 0 if the line is vertical
      // In that case, the distance is the distance to the line
      const Vector_3 w = cross_product(pl - p,
                                       vl);
      return (w * w) / (vl * vl);
    }
    else {
      const FT det = determinant(traits.normal(), 
                                 vl, 
                                 pl - p);
      return (det * det) / ( v * v );
    }
  }
}; // end class Projected_square_distance_to_line_with_normal_3

template <class Traits>
class Projected_intersect_3
{
  // private members
  const Traits& traits;
  
  // private types aliases
  typedef typename Traits::K K;
  typedef typename Traits::Point_2 Point;
  typedef typename Traits::Line_2 Line;
  typedef typename Traits::Segment_2 Segment;
  typedef typename K::Plane_3 Plane_3;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::FT FT;
public:
  Projected_intersect_3(const Traits& traits_)
    : traits(traits_)
  {}

  Object operator()(const Segment& s1, const Segment& s2)
  {
    CGAL_TIME_PROFILER("Projected_intersect_3")
    const Vector_3 u1 = cross_product(s1.to_vector(), traits.normal());
    if(u1 == NULL_VECTOR)
      return K().intersect_3_object()(s1.supporting_line(), s2);

    const Vector_3 u2 = cross_product(s2.to_vector(), traits.normal());
    if(u2 == NULL_VECTOR)
      return K().intersect_3_object()(s1, s2.supporting_line());
    
    const Plane_3 plane_1(s1.source(), u1);
    const Plane_3 plane_2(s2.source(), u2);
    
    Object planes_intersection = intersection(plane_1, plane_2);
    if(planes_intersection.empty()) {
      std::cerr << "planes_intersection is empty\n";
      return planes_intersection;
    }
    if(const Line* line = object_cast<Line>(&planes_intersection))
    {
      const Point& pi = line->point(0);
      if(cross_product(traits.normal(), pi - s1.source())
         * cross_product(traits.normal(), pi - s1.target()) > FT(0)
         ||
         cross_product(traits.normal(), pi - s2.source())
         * cross_product(traits.normal(), pi - s2.target()) > FT(0) )
      {
        // the intersection of the lines is not inside the segments
        std::cerr << "intersection not inside\n";
        return Object();
      }
      else
      {
        // Let the plane passing through s1.source() and with normal
        // the cross product of s1.to_vector() and s2.to_vector(). That
        // plane should intersect *l, now.
        return intersection(*line, Plane_3(s1.source(), 
                                           cross_product(s1.to_vector(),
                                                         s2.to_vector())));
      }
    }
    if(object_cast<Plane_3>(&planes_intersection))
    {
      std::cerr << "coplanar lines\n";
      CGAL_error();
      return Object();
    }
    return Object();
  }
}; // end class Projected_intersect_3

} // end namespace TriangulationProjectionTraitsCartesianFunctors


template < class Kernel >
class Triangulation_2_projection_traits_3
{
  typedef Triangulation_2_projection_traits_3<Kernel> Self;

  typename Kernel::Vector_3 n;

public:
  typedef typename Kernel::Vector_3 Vector_3;


  Triangulation_2_projection_traits_3(const Vector_3& n_)
    : n(n_)
  {}

  const Vector_3& normal() const
  {
    return n;
  }

  Triangulation_2_projection_traits_3(const Self& other)
    : n(other.n)
  {}

  Self& operator=(const Self& other)
  {
    if(this != &other) {
      n = other.n;
    }
    return *this;
  }

  typedef Kernel K;
  typedef typename K::FT          FT;
  typedef typename K::Point_3     Point_2;
  typedef typename K::Segment_3   Segment_2;
  typedef typename K::Triangle_3  Triangle_2;
  typedef typename K::Line_3      Line_2;

  // Maybe not a good choice
  typedef typename K::Less_xy_3            Less_x_2;
  typedef typename K::Less_z_3             Less_y_2;

  typedef typename K::Compare_xy_3                           Compare_x_2;
  typedef typename K::Compare_z_3                            Compare_y_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
    Projected_orientation_with_normal_3<Self>                Orientation_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
    Projected_side_of_oriented_circle_with_normal_3<Self>    Side_of_oriented_circle_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
  Projected_square_distance_to_line_with_normal_3<Self>      Compute_squared_distance_2;

  typedef TriangulationProjectionTraitsCartesianFunctors::
  Projected_intersect_3<Self>                                Intersect_2;

  typedef typename K::Construct_segment_3  Construct_segment_2;
  typedef typename K::Construct_line_3     Construct_line_2;
  typedef typename K::Construct_triangle_3 Construct_triangle_2;
  
  Less_x_2
  less_x_2_object() const
    { return Less_x_2();}

  Less_y_2
  less_y_2_object() const
    { return Less_y_2();}

  Compare_x_2
  compare_x_2_object() const
  {
    return Compare_x_2();
  }

  Compare_y_2
  compare_y_2_object() const
  { 
    return Compare_y_2();
  }

  Orientation_2 
  orientation_2_object() const
  {
    return Orientation_2(*this); 
  }

  Side_of_oriented_circle_2
  side_of_oriented_circle_2_object() const
  {
    return Side_of_oriented_circle_2(*this);
  }
  
  Compute_squared_distance_2
  compute_squared_distance_2_object() const
  {
    return Compute_squared_distance_2(*this);
  }
  
  Intersect_2
  intersect_2_object () const
  {
    return Intersect_2(*this);
  }

  Construct_segment_2  construct_segment_2_object() const
    {return Construct_segment_2();}
  
  Construct_line_2  construct_line_2_object() const
    {return Construct_line_2();}
  
  Construct_triangle_2  construct_triangle_2_object() const
    {return Construct_triangle_2();}

}; // end class Triangulation_2_projection_traits_3<Kernel>
  

} // end namespace CGAL

#endif // CGAL_TRIANGULATION_2_PROJECTION_TRAITS_3_H

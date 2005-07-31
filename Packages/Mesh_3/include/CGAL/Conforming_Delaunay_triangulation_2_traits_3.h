// Copyright (c) 2004  INRIA Sophia-Antipolis (France).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_CONFORMING_DELAUNAY_TRIANGULATION_2_TRAITS_3_H
#define CGAL_CONFORMING_DELAUNAY_TRIANGULATION_2_TRAITS_3_H

#include <CGAL/Triangulation_2_traits_3.h>
namespace CGAL {

template<class K>
struct Conforming_Delaunay_triangulation_2_traits_3
  : public CGAL::Triangulation_2_traits_3<K>
{
  typedef CGAL::Triangulation_2_traits_3<K> Base;
  typedef typename Base::Rep Rep;
  typedef typename Base::Point_2 Point_2;
  typedef typename Base::Segment_2 Segment_2;
  typedef typename Base::Triangle_2 Triangle_2;
 
  typedef typename Base::Compare_x_2 Compare_x_2;
  typedef typename Base::Compare_y_2 Compare_y_2;
  typedef typename Base::Orientation_2 Orientation_2;
  typedef typename Base::Side_of_oriented_circle_2 Side_of_oriented_circle_2;  
  typedef typename Base::Construct_segment_2 Construct_segment_2;
  typedef typename Base::Construct_triangle_2 Construct_triangle_2;

  // for ConformingDelaunayTriangulationTraits_2<Tr>
  typedef typename Rep::FT FT;
  typedef typename Rep::Vector_3 Vector_2;
  typedef typename Rep::Construct_vector_3 Construct_vector_2;
  typedef typename Rep::Construct_scaled_vector_3 Construct_scaled_vector_2;
  typedef typename Rep::Construct_translated_point_3 
                                   Construct_translated_point_2;
  typedef typename Rep::Construct_midpoint_3 Construct_midpoint_2;
  typedef typename Rep::Compute_squared_distance_3 Compute_squared_distance_2;
  typedef typename Rep::Angle_3 Angle_2;

  // for DelaunayMeshTraits_2<Tr>
  typedef typename Rep::Construct_circumcenter_3 Construct_circumcenter_2;

  Construct_vector_2
  construct_vector_2_object() const
  { return Construct_vector_2(); }

  Construct_scaled_vector_2
  construct_scaled_vector_2_object() const
  { return Construct_scaled_vector_2(); }

  Construct_translated_point_2
  construct_translated_point_2_object() const
  { return Construct_translated_point_2(); }

  Construct_midpoint_2
  construct_midpoint_2_object() const
  { return Construct_midpoint_2(); }

  Compute_squared_distance_2
  compute_squared_distance_2_object() const
  { return Compute_squared_distance_2(); }

  Angle_2
  angle_2_object() const
  { return Angle_2(); }

  Construct_circumcenter_2
  construct_circumcenter_2_object() const
  { return Construct_circumcenter_2(); }
};

} // end namespace CGAL

#endif //CGAL_CONFORMING_DELAUNAY_TRIANGULATION_2_TRAITS_3_H

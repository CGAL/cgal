// Copyright (c) 2006-2009   INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_TRIANGULATION_TRAITS_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Periodic_3_offset_3.h>
#include <CGAL/Periodic_3_construct_point_3.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Traits_with_offsets_adaptor.h>


namespace CGAL { 

template < class Kernel, class Off = typename CGAL::Periodic_3_offset_3 >
class Periodic_3_triangulation_traits_base_3
  : public Kernel
{
public:
  typedef Kernel                                                 K;
  typedef Off                                                    Offset;
  typedef Periodic_3_triangulation_traits_base_3< K, Offset >    Self;  

  typedef typename K::RT                RT;
  typedef typename K::FT                FT;
  typedef typename K::Point_3           Point_3;
  typedef typename K::Vector_3          Vector_3;
  typedef Offset                        Periodic_3_offset_3;
  typedef typename K::Iso_cuboid_3      Iso_cuboid_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3 Point;

  typedef typename K::Segment_3         Segment_3;
  typedef typename K::Triangle_3        Triangle_3;
  typedef typename K::Tetrahedron_3     Tetrahedron_3;

  // Triangulation predicates
  typedef Traits_with_offsets_adaptor<Self, typename K::Compare_xyz_3>
      Compare_xyz_3;
  typedef Traits_with_offsets_adaptor<Self, typename K::Orientation_3>
      Orientation_3;
  
  // Triangulation constructions
  typedef Periodic_3_construct_point_3<Self, typename K::Construct_point_3>
      Construct_point_3;
  typedef Traits_with_offsets_adaptor<Self, typename K::Construct_segment_3>
      Construct_segment_3;
  typedef Traits_with_offsets_adaptor<Self, typename K::Construct_triangle_3>
      Construct_triangle_3;
  typedef Traits_with_offsets_adaptor<Self, typename K::Construct_tetrahedron_3>
      Construct_tetrahedron_3;

  // Access
  void set_domain(const Iso_cuboid_3& domain) {
    _domain = domain;
  }
  
  Iso_cuboid_3 get_domain() const {
    return _domain;
  }

  // Operations
  Compare_xyz_3
  compare_xyz_3_object() const {
    return Compare_xyz_3(&_domain);
  }
  Orientation_3
  orientation_3_object() const {
    return Orientation_3(&_domain);
  }
  Construct_point_3
  construct_point_3_object() const {
    return Construct_point_3(_domain);
  }
  Construct_segment_3
  construct_segment_3_object() const {
    return Construct_segment_3(&_domain);
  }
  Construct_triangle_3
  construct_triangle_3_object() const {
    return Construct_triangle_3(&_domain);
  }
  Construct_tetrahedron_3
  construct_tetrahedron_3_object() const {
    return Construct_tetrahedron_3(&_domain);
  }

protected:
  Iso_cuboid_3 _domain;
};

namespace future_release
{
template < typename K, typename Off = CGAL::Periodic_3_offset_3 >
class Periodic_3_triangulation_traits_3;
}

} //namespace CGAL

// Partial specialization for Filtered_kernel<CK>.
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_filtered_traits_3.h>

namespace CGAL {
// This declaration is needed to break the cyclic dependency.
template < typename K, typename Off >
class Periodic_3_triangulation_filtered_traits_3;

namespace future_release {
template < class K, class Off>
class Periodic_3_triangulation_traits_3
  : public Periodic_3_triangulation_traits_base_3<K, Off>
{
};

template < typename CK, typename Off >
class Periodic_3_triangulation_traits_3 < Filtered_kernel<CK>, Off>
  : public Periodic_3_triangulation_filtered_traits_3 <
  Filtered_kernel<CK>, Off >
{
public:
  typedef Filtered_kernel<CK>  Kernel;
};

template < class Off >
class Periodic_3_triangulation_traits_3<CGAL::Epick, Off>
  : public Periodic_3_triangulation_filtered_traits_3<CGAL::Epick, Off>
{
  typedef CGAL::Epick Kernel;
};

}
} //namespace CGAL

#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>


namespace CGAL
{
template < typename K, typename Off >
class Periodic_3_Delaunay_triangulation_traits_3;

// Periodic_3_triangulation_traits_3 should not be used as traits for Periodic_3_Delaunay_triangulation_3 anymore.
template < class Kernel, class Off = typename CGAL::Periodic_3_offset_3 >
class CGAL_DEPRECATED Periodic_3_triangulation_traits_3 : public Periodic_3_Delaunay_triangulation_traits_3<Kernel, Off>
{
};
}

#endif // CGAL_PERIODIC_3_TRIANGULATION_TRAITS_3_H

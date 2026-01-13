// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Léo Valque

#ifndef CGAL_THREE_POINT_CUT_PLANE_TRAITS_H
#define CGAL_THREE_POINT_CUT_PLANE_TRAITS_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
\ingroup PMP_corefinement_grp

The class `Three_point_cut_plane_traits<Kernel>` is a model of the
TODO concept. This traits class represents planes by three defining points and supports efficient point–plane orientation and plane–line intersection operations.

TODO documented or remains internal?

*/
template <class Kernel>
struct Three_point_cut_plane_traits
{
  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  struct Plane_3: public std::array<typename Kernel::Point_3, 3>{
    using Base = std::array<typename Kernel::Point_3, 3>;
    using Explicit_plane = typename Kernel::Plane_3;

    Plane_3(const Point_3 &a, const Point_3 &b, const Point_3 &c): Base({a, b, c}){}
    Plane_3(const std::array<Point_3, 3> &arr): Base(arr){}

    // Warning: it is slow (Planes are constructed each time)
    bool operator<(const Plane_3 &b) const{
      Explicit_plane pa = explicit_plane();
      Explicit_plane pb = b.explicit_plane();
      Comparison_result res = compare(pa.a(), pb.a());
      if(res == EQUAL)
        res = compare(pa.b(), pb.b());
      if(res == EQUAL)
        res = compare(pa.c(), pb.c());
      if(res == EQUAL)
        res = compare(pa.d(), pb.d());
      return res == SMALLER;
    };

    bool operator==(const Plane_3 &b) const{
      Explicit_plane pa = explicit_plane();
      Explicit_plane pb = b.explicit_plane();
      return pa==pb;
    }

    Explicit_plane explicit_plane() const{
      return  Explicit_plane((*this)[0], (*this)[1], (*this)[2]);
    }

  };
  using Vector_3 = typename Kernel::Vector_3;

  struct Does_not_support_CDT2{};

  struct Oriented_side_3
  {
    Oriented_side operator()(const Plane_3& plane, const Point_3& p)  const
    {
      return orientation(plane[0], plane[1], plane[2], p);
    }
  };

  struct Construct_plane_line_intersection_point_3
  {
    Point_3 operator()(const Plane_3& plane, const Point_3& p, const Point_3& q)
    {
      typename Kernel::Construct_plane_line_intersection_point_3 construction;
      return construction(plane[0], plane[1], plane[2], p, q);
    }
  };

  struct Construct_orthogonal_vector_3{
    Vector_3 operator()(const Plane_3& plane)
    {
      return typename Kernel::Plane_3(plane[0], plane[1], plane[2]).orthogonal_vector();
    }
  };

  struct Compute_squared_distance_3
  {
    using Compute_scalar_product_3 = typename Kernel::Compute_scalar_product_3;
    FT operator()(const Plane_3& plane, const Point_3& p)
    {
      typename Kernel::Plane_3 pl(plane[0], plane[1], plane[2]);
      return Compute_scalar_product_3()(Vector_3(ORIGIN, p), pl.orthogonal_vector())+pl.d();
    }
  };

  Oriented_side_3 oriented_side_3_object() const
  {
    return Oriented_side_3();
  }

  Construct_plane_line_intersection_point_3 construct_plane_line_intersection_point_3_object() const
  {
    return Construct_plane_line_intersection_point_3();
  }

  Construct_orthogonal_vector_3 construct_orthogonal_vector_3_object() const
  {
    return Construct_orthogonal_vector_3();
  }

  Compute_squared_distance_3 compute_squared_distance_3_object() const { return Compute_squared_distance_3(); }

};

} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_THREE_POINT_CUT_PLANE_TRAITS_H
// Copyright (c) 2004-2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_POISSON_MESH_CRITERIA_3_H
#define CGAL_POISSON_MESH_CRITERIA_3_H

#include <CGAL/license/Poisson_surface_reconstruction_3.h>


#include <iostream>
#include <CGAL/number_utils.h>
#include <CGAL/utils.h> // for min

namespace CGAL {

namespace internal {
namespace Poisson {

template <class Tr>
class Constant_sizing_field {
  typedef typename Tr::FT FT;
  FT sq_radius_bound;
public:
  FT cell_radius_bound() const { return CGAL::approximate_sqrt(sq_radius_bound); }

  Constant_sizing_field(FT sq_radius_bound = 0.)
    : sq_radius_bound(sq_radius_bound) {}

  template <typename Point>
  FT operator()(const Point&) const { return sq_radius_bound; }
}; // end class Constant_sizing_field

} // end namespace Poisson
} // end namespace internal

template <
  class Tr,
  typename Sizing_field = internal::Poisson::Constant_sizing_field<Tr>,
  typename Second_sizing_field = internal::Poisson::Constant_sizing_field<Tr>
  >
class Poisson_mesh_cell_criteria_3
{
  const Sizing_field& sq_sizing_field;
  const Second_sizing_field& sq_second_sizing_field;
  double sq_radius_edge_bound_;

  typedef Poisson_mesh_cell_criteria_3<Tr,
                                       Sizing_field,
                                       Second_sizing_field> Self;
public:
  struct Cell_quality
  {
    double aspect_ratio_; // squared radius edge ratio (square(3))
    double sq_size_;

    Cell_quality() {}
    Cell_quality(double _aspect, double _sq_size)
    : aspect_ratio_(_aspect)
    , sq_size_(_sq_size)
    {}

    double sq_size() const { return sq_size_; }
    double aspect() const { return aspect_ratio_; }

    // q1<q2 means q1 is prioritized over q2
    // ( q1 == *this, q2 == q )
    bool operator<(const Cell_quality& q) const
    {
      if( sq_size() > 1 )
        if( q.sq_size() > 1 )
          return ( sq_size() > q.sq_size() );
        else
          return true; // *this is big but not q
      else
        if( q.sq_size() >  1 )
          return false; // q is big but not *this
      return( aspect() > q.aspect() );
    }
  };

  typedef typename Tr::Cell_handle Cell_handle;

  Poisson_mesh_cell_criteria_3(const double sq_radius_edge_bound = 2, //< radius edge ratio bound (ignored if zero)
                  const double sq_radius_bound = 0) //< squared cell radius bound (ignored if zero)
    : sq_sizing_field(sq_radius_bound),
      sq_second_sizing_field(),
      sq_radius_edge_bound_(sq_radius_edge_bound)
  {}

  Poisson_mesh_cell_criteria_3(const double sq_radius_edge_bound,   //< radius edge ratio bound (ignored if zero)
                               const Sizing_field& _sq_sizing_field, ///< squared sizing field for cell radius bound
                               const Second_sizing_field& _sq_second_sizing_field)
    : sq_sizing_field(_sq_sizing_field),
      sq_second_sizing_field(_sq_second_sizing_field),
      sq_radius_edge_bound_(sq_radius_edge_bound)
  {}

  inline
  double sq_radius_edge_bound() const
  {
    return sq_radius_edge_bound_;
  }

  inline
  void set_sq_radius_edge_bound(const double sq_radius_edge_bound)
  {
    sq_radius_edge_bound_ = sq_radius_edge_bound;
  }

  class Is_bad
  {
  protected:
    const Self* cell_criteria;
  public:
    typedef typename Tr::Point Point_3;

    Is_bad(const Self* cell_criteria)
      : cell_criteria(cell_criteria) {}

    bool operator()(const Cell_handle& c,
                    Cell_quality& qual) const
    {
      const Point_3& p = c->vertex(0)->point();
      const Point_3& q = c->vertex(1)->point();
      const Point_3& r = c->vertex(2)->point();
      const Point_3& s = c->vertex(3)->point();

      typedef typename Tr::Geom_traits Geom_traits;
      typedef typename Geom_traits::Compute_squared_radius_3 Sq_Radius;
      typedef typename Geom_traits::Compute_squared_distance_3 Sq_Distance;
      typedef typename Geom_traits::Compute_scalar_product_3 Scalar_product;

      Sq_Radius sq_radius = Geom_traits().compute_squared_radius_3_object();
      Sq_Distance sqd = Geom_traits().compute_squared_distance_3_object();

      double sq_size = to_double(sq_radius(p, q, r, s));

      // If the tetrahedron is small enough according to the second sizing
      // field, then let's say the tetrahedron is OK according to the
      // sizing fields.
      if( sq_size > cell_criteria->sq_second_sizing_field(p) )
      {
        // first sizing field
        const double sq_size_bound = cell_criteria->sq_sizing_field(p);
        if(sq_size_bound > 0.) {
          qual.sq_size_ = sq_size / sq_size_bound;
          // normalized by size bound to deal
          // with size field
          if( qual.sq_size() > 1 )
          {
            qual.aspect_ratio_ = 1; // (do not compute aspect)
            return true;
          }
        }
      }
      if(cell_criteria->sq_radius_edge_bound() == 0)
      {
        qual = Cell_quality{0,1};
        return false;
      }

      // Aspect ratio
      double min_sq_length = CGAL::to_double(sqd(p, q));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(sqd(p, r)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(sqd(p, s)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(sqd(q, r)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(sqd(q, s)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(sqd(r, s)));

      qual.aspect_ratio_ = sq_size / min_sq_length;

      if(qual.aspect_ratio_ > cell_criteria->sq_radius_edge_bound_)
        return true; // bad aspect ratio

      // Check normals
      Scalar_product scalar_product = Geom_traits().compute_scalar_product_3_object();

      const std::array<std::pair<int, int>, 6> vs = {{{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}};
      for(const auto& vv : vs)
      {
        auto v0 = c->vertex(vv.first);
        auto v1 = c->vertex(vv.second);
        // skip if one of the vertices is not an input point
        if(v0->type() != 0 || v1->type() != 0)
          continue;
        if(scalar_product(v0->normal(), v1->normal()) <= 0.)
          return true; // normals are not oriented in the same direction
      }

      return false;
    }

  }; // end Is_bad


  Is_bad is_bad_object() const
  { return Is_bad(this); }

}; // end Poisson_mesh_cell_criteria_3

  template <typename Tr>
  std::ostream& operator<<(std::ostream& os,
                           const typename Poisson_mesh_cell_criteria_3<Tr>::Cell_quality& q)
  {
    return os << q.sq_size() << ", " << q.aspect();
  }

} // end namespace CGAL

#endif // CGAL_POISSON_MESH_CRITERIA_3_H

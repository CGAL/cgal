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
  double sq_radius_bound;
public:
  double cell_radius_bound() const { return CGAL::sqrt(sq_radius_bound); }

  Constant_sizing_field(double sq_radius_bound = 0.)
    : sq_radius_bound(sq_radius_bound) {}

  template <typename Point>
  double operator()(const Point&) const { return sq_radius_bound; }
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
  Sizing_field sizing_field;
  Second_sizing_field second_sizing_field;
  double radius_edge_bound_;

  typedef Poisson_mesh_cell_criteria_3<Tr,
                                       Sizing_field,
                                       Second_sizing_field> Self;
public:
  struct Cell_quality : public std::pair<double, double>
  {
    typedef std::pair<double, double> Base;

    Cell_quality() : Base() {}
    Cell_quality(double _aspect, double _sq_size) : Base(_aspect, _sq_size) {};

    double sq_size() const { return second; }
    double aspect() const { return first; }

    // q1<q2 means q1 is prioritised over q2
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

  // inline
  // double squared_radius_bound() const
  // {
  //   return squared_radius_bound_;
  // }

  typedef typename Tr::Cell_handle Cell_handle;

  Poisson_mesh_cell_criteria_3(const double radius_edge_bound = 2, //< radius edge ratio bound (ignored if zero)
                  const double radius_bound = 0) //< cell radius bound (ignored if zero)
    : sizing_field(radius_bound*radius_bound),
      second_sizing_field(),
      radius_edge_bound_(radius_edge_bound)
  {}

  Poisson_mesh_cell_criteria_3(const double radius_edge_bound, //< radius edge ratio bound (ignored if zero)
                               const Sizing_field& sizing_field,
                               const Second_sizing_field second_sizing_field = Second_sizing_field())
    : sizing_field(sizing_field),
      second_sizing_field(second_sizing_field),
      radius_edge_bound_(radius_edge_bound)
  {}

  // inline
  // void set_squared_radius_bound(const double squared_radius_bound)
  // {
  //   squared_radius_bound_ = squared_radius_bound;
  // }

  inline
  double radius_edge_bound() const
  {
    return radius_edge_bound_;
  }

  inline
  void set_radius_edge_bound(const double radius_edge_bound)
  {
    radius_edge_bound_ = radius_edge_bound;
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
      typedef typename Geom_traits::Compute_squared_radius_3 Radius;
      typedef typename Geom_traits::Compute_squared_distance_3 Distance;

      Radius sq_radius = Geom_traits().compute_squared_radius_3_object();
      Distance distance = Geom_traits().compute_squared_distance_3_object();

      double size = to_double(sq_radius(p, q, r, s));

      // If the tetrahedron is small enough according to the second sizing
      // field, then let's say the tetrahedron is OK according to the
      // sizing fields.
      if( size > cell_criteria->second_sizing_field(p) )
      {
        // first sizing field
        const double size_bound = cell_criteria->sizing_field(p);
        if(size_bound > 0.) {
          qual.second = size / size_bound;
          // normalized by size bound to deal
          // with size field
          if( qual.sq_size() > 1 )
          {
            qual.first = 1; // (do not compute aspect)
            return true;
          }
        }
      }
      if( cell_criteria->radius_edge_bound_ == 0 )
        {
          qual = Cell_quality(0,1);
          return false;
        }

      double min_sq_length = CGAL::to_double(distance(p, q));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(p, r)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(p, s)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(q, r)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(q, s)));
      min_sq_length = (CGAL::min)(min_sq_length, to_double(distance(r, s)));

      qual.first = size / min_sq_length;

      return (qual.first > cell_criteria->radius_edge_bound_);
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

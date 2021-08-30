// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_SEARCH_TRAITS_VERTEX_HANDLE_3_H
#define CGAL_SEARCH_TRAITS_VERTEX_HANDLE_3_H

#include <CGAL/license/Point_set_processing_3.h>


#include <CGAL/Kernel_traits.h>
#include <CGAL/Search_traits.h>
#include <CGAL/Kd_tree_rectangle.h>

#include <cmath>

namespace CGAL {

  namespace internal {

/// \cond SKIP_IN_MANUAL

/// A Point_vertex_handle_3 objects wraps either
/// a Vertex_handle or a 3D point.
template <class Vertex_handle>
struct Point_vertex_handle_3
{
private:

  double m_coord[3];
  Vertex_handle m_vertex_handle;

public:

  Point_vertex_handle_3()
  {
    m_coord[0] =
    m_coord[1] =
    m_coord[2] = 0;
    //m_vertex_handle = nullptr;
  }

  Point_vertex_handle_3(double x, double y, double z,
                        Vertex_handle vertex_handle) ///< nullptr for query points
  {
    m_coord[0] = x;
    m_coord[1] = y;
    m_coord[2] = z;
    m_vertex_handle = vertex_handle;
  }

  /// Default copy constructor, operator =() and destructor are fine.

  operator Vertex_handle() const { return m_vertex_handle; }

  double x() const { return m_coord[0]; }
  double y() const { return m_coord[1]; }
  double z() const { return m_coord[2]; }

  double& x() { return m_coord[0]; }
  double& y() { return m_coord[1]; }
  double& z() { return m_coord[2]; }

  bool operator==(const Point_vertex_handle_3& p) const
  {
    return (x() == p.x()) &&
           (y() == p.y()) &&
           (z() == p.z());
  }

  bool  operator!=(const Point_vertex_handle_3& p) const { return ! (*this == p); }

  template <class V>
  friend struct Construct_cartesian_const_iterator_vertex_handle_3;

}; // end of class Point_vertex_handle_3

} // namespace internal

/// Kernel traits specialization for Point_vertex_handle_3.
template <class Vertex_handle>
struct Kernel_traits< internal::Point_vertex_handle_3<Vertex_handle> > {
  struct Kernel {
    typedef double FT;
    typedef double RT;
  };
};

  namespace internal {

/// Functor with two function operators, which return the begin and
/// past the end iterator for the Cartesian coordinates.
template <class Vertex_handle>
struct Construct_cartesian_const_iterator_vertex_handle_3
{
  typedef CGAL::internal::Point_vertex_handle_3<Vertex_handle> Point_vertex_handle_3;
  typedef const double* result_type;
  const double* operator()(const Point_vertex_handle_3& p) const
  { return static_cast<const double*>(p.m_coord); }

  const double* operator()(const Point_vertex_handle_3& p, int)  const
  { return static_cast<const double*>(p.m_coord+3); }
};

/// The class Euclidean_distance_vertex_handle_3 implements the Euclidean distance
/// for Vertex_handles.
/// To optimize distance computations squared distances are used.
///
/// @heading Is Model for the Concepts: Model of the OrthogonalDistance concept.
template <class Vertex_handle>
struct Euclidean_distance_vertex_handle_3
{
  typedef double FT;
  typedef CGAL::internal::Point_vertex_handle_3<Vertex_handle> Point_vertex_handle_3;
  typedef Point_vertex_handle_3 Query_item;
  typedef Point_vertex_handle_3 Point_d;

  double transformed_distance(const Point_vertex_handle_3& p1, const Point_vertex_handle_3& p2) const {
    double distx = p1.x()-p2.x();
    double disty = p1.y()-p2.y();
    double distz = p1.z()-p2.z();
    return distx * distx +
           disty * disty +
           distz * distz;
  }


  double min_distance_to_rectangle(const Point_vertex_handle_3& p,
                                   const Kd_tree_rectangle<double, Dimension_tag<3> >& b) const {
    double distance(0.0), h = p.x();
    if (h < b.min_coord(0)) distance += (b.min_coord(0)-h)*(b.min_coord(0)-h);
    if (h > b.max_coord(0)) distance += (h-b.max_coord(0))*(h-b.max_coord(0));
    h=p.y();
    if (h < b.min_coord(1)) distance += (b.min_coord(1)-h)*(b.min_coord(1)-h);
    if (h > b.max_coord(1)) distance += (h-b.max_coord(1))*(h-b.min_coord(1));
    h=p.z();
    if (h < b.min_coord(2)) distance += (b.min_coord(2)-h)*(b.min_coord(2)-h);
    if (h > b.max_coord(2)) distance += (h-b.max_coord(2))*(h-b.max_coord(2));
    return distance;
  }


  double min_distance_to_rectangle(const Point_vertex_handle_3& p,
                                   const Kd_tree_rectangle<double,Dimension_tag<3> >& b,
                                   std::vector<double>& dists) const {
    double distance(0.0), h = p.x();
    if (h < b.min_coord(0)){
      dists[0] = (b.min_coord(0)-h);
      distance += dists[0] * dists[0];
    }
    else if (h > b.max_coord(0)){
      dists[0] = (h-b.max_coord(0));
      distance += dists[0] * dists[0];
    }

    h=p.y();
    if (h < b.min_coord(1)){
      dists[1] = (b.min_coord(1)-h);
      distance += dists[1] * dists[1];
    }
    else if (h > b.max_coord(1)){
      dists[1] = (h-b.max_coord(1));
      distance += dists[1] * dists[1];
    }

    h=p.z();
    if (h < b.min_coord(2)){
      dists[2] = (b.min_coord(2)-h);
      distance += dists[2] * dists[2];
    }
    else if (h > b.max_coord(2)){
      dists[2] = (h-b.max_coord(2));
      distance += dists[2] * dists[2];
    }
    return distance;
  }


  double max_distance_to_rectangle(const Point_vertex_handle_3& p,
                                   const Kd_tree_rectangle<double,Dimension_tag<3> >& b) const {
    double h = p.x();

    double d0 = (h >= (b.min_coord(0)+b.max_coord(0))/2.0) ?
                (h-b.min_coord(0))*(h-b.min_coord(0)) : (b.max_coord(0)-h)*(b.max_coord(0)-h);

    h=p.y();
    double d1 = (h >= (b.min_coord(1)+b.max_coord(1))/2.0) ?
                (h-b.min_coord(1))*(h-b.min_coord(1)) : (b.max_coord(1)-h)*(b.max_coord(1)-h);
    h=p.z();
    double d2 = (h >= (b.min_coord(2)+b.max_coord(2))/2.0) ?
                (h-b.min_coord(2))*(h-b.min_coord(2)) : (b.max_coord(2)-h)*(b.max_coord(2)-h);
    return d0 + d1 + d2;
  }



  double max_distance_to_rectangle(const Point_vertex_handle_3& p,
                                   const Kd_tree_rectangle<double, Dimension_tag<3> >& b,
                                   std::vector<double>& dists) const {
    double h = p.x();
    dists[0] = (h >= (b.min_coord(0)+b.max_coord(0))/2.0)? h-b.min_coord(0) : b.max_coord(0)-h;
    double d0 = dists[0]*dists[0];

    h=p.y();
    dists[1] = (h >= (b.min_coord(1)+b.max_coord(1))/2.0)? h-b.min_coord(1): b.max_coord(1)-h;
    double d1 = dists[1]*dists[1];

    h=p.z();
    dists[2] = (h >= (b.min_coord(2)+b.max_coord(2))/2.0)? h-b.min_coord(2) : b.max_coord(2)-h;

    double d2 = dists[2]*dists[2];

    return d0 + d1 + d2;
  }


  double new_distance(double& dist, double old_off, double new_off, int /*cutting_dimension*/) const
  {
    return dist + new_off*new_off - old_off*old_off;
  }

  double transformed_distance(double d) const { return d*d; }

  double inverse_of_transformed_distance(double d) { return CGAL::sqrt(d); }

}; // end of struct Euclidean_distance_vertex_handle_3



/// Traits class for the Spatial_searching package
/// for Vertex_handles.
///
/// @heading Is Model for the Concepts: Model of the SearchTraits concept.
template <class Vertex_handle>
class Search_traits_vertex_handle_3
  : public Search_traits< double,
                          Point_vertex_handle_3<Vertex_handle>,
                          const double*,
                          Construct_cartesian_const_iterator_vertex_handle_3<Vertex_handle> ,
                          Dimension_tag<3>
                          >
{};

  } // namespace internal
/// \endcond

} //namespace CGAL

#endif // CGAL_SEARCH_TRAITS_VERTEX_HANDLE_3_H

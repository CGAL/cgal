// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Clement Jamin

#ifndef CGAL_REGULAR_TRIANGULATION_TRAITS_ADAPTER_H
#define CGAL_REGULAR_TRIANGULATION_TRAITS_ADAPTER_H

#include <CGAL/license/Triangulation.h>

#include <CGAL/basic.h>

#include <CGAL/boost/iterator/transform_iterator.hpp>

namespace CGAL {

// Wrapper class to make a model of `RegularTriangulationTraits` easily usable
// by the `Regular_triangulation` class. By using this class:
// - Point_d (used by `Triangulation` and the TDS) becomes a weighted point
// - Predicates and functors such as Less_coordinate_d or Orientation_d
//   can be called using weighted points instead of bare points (this is
//   needed because `Weighted_point_d` is not convertible to `Point_d`)
// This way, `Triangulation` works perfectly well with weighted points.

template <class RTTraits>
class Regular_triangulation_traits_adapter
  : public RTTraits
{
public:
  typedef RTTraits                                     Base;

  // Required by TriangulationTraits
  typedef typename Base::Dimension                     Dimension;
  typedef typename Base::FT                            FT;
  typedef typename Base::Flat_orientation_d            Flat_orientation_d;
  typedef typename Base::Weighted_point_d              Point_d;

  // Needed by Regular_triangulation (picked from RegularTriangulationTraits)
  typedef typename Base::Weighted_point_d              Weighted_point_d;
  typedef typename Base::Construct_point_d             Construct_point_d;
  typedef typename Base::Compute_weight_d              Compute_weight_d;
  typedef typename Base::Power_side_of_power_sphere_d  Power_side_of_power_sphere_d;
  typedef typename Base::In_flat_power_side_of_power_sphere_d 
                                                       In_flat_power_side_of_power_sphere_d;

  //===========================================================================
  // Custom types
  //===========================================================================

  // Required by SpatialSortingTraits_d
  class Less_coordinate_d
  {
    const RTTraits &m_traits;

  public:
    typedef bool result_type;

    Less_coordinate_d(const RTTraits &kernel)
      : m_traits(kernel) {}

    result_type operator()(
      Weighted_point_d const& p, Weighted_point_d const& q, int i) const
    {
      Construct_point_d cp = m_traits.construct_point_d_object();
      return m_traits.less_coordinate_d_object() (cp(p), cp(q), i);
    }
  };

  //===========================================================================

  // Required by TriangulationTraits
  class Orientation_d
  {
    const RTTraits &m_traits;

  public:
    typedef Orientation result_type;

    Orientation_d(const RTTraits &kernel)
      : m_traits(kernel) {}

    template <typename ForwardIterator> 
    result_type operator()(ForwardIterator start, ForwardIterator end) const
    {
      Construct_point_d cp = m_traits.construct_point_d_object();
      return m_traits.orientation_d_object() (
        boost::make_transform_iterator(start, cp),
        boost::make_transform_iterator(end, cp)
      );
    }
  };

  //===========================================================================

  // Required by TriangulationTraits
  class Construct_flat_orientation_d
  {
    const RTTraits &m_traits;

  public:
    typedef Flat_orientation_d result_type;
    
    Construct_flat_orientation_d(const RTTraits &kernel)
      : m_traits(kernel) {}

    template <typename ForwardIterator> 
    result_type operator()(ForwardIterator start, ForwardIterator end) const
    {
      Construct_point_d cp = m_traits.construct_point_d_object();
      return m_traits.construct_flat_orientation_d_object() (
        boost::make_transform_iterator(start, cp),
        boost::make_transform_iterator(end, cp)
      );
    }
  };


  //===========================================================================

  // Required by TriangulationTraits
  class In_flat_orientation_d
  {
    const RTTraits &m_traits;

  public:
    typedef Orientation result_type;
    
    In_flat_orientation_d(const RTTraits &kernel)
      : m_traits(kernel) {}

    template <typename ForwardIterator> 
    result_type operator()(Flat_orientation_d orient, 
      ForwardIterator start, ForwardIterator end) const
    {
      Construct_point_d cp = m_traits.construct_point_d_object();
      return m_traits.in_flat_orientation_d_object() (
        orient,
        boost::make_transform_iterator(start, cp),
        boost::make_transform_iterator(end, cp)
      );
    }
  };

  //===========================================================================
  
  // Required by TriangulationTraits
  class Contained_in_affine_hull_d
  {
    const RTTraits &m_traits;

  public:
    typedef bool result_type;
    
    Contained_in_affine_hull_d(const RTTraits &kernel)
      : m_traits(kernel) {}

    template <typename ForwardIterator> 
    result_type operator()(ForwardIterator start, ForwardIterator end, 
                           const Weighted_point_d & p) const
    {
      Construct_point_d cp = m_traits.construct_point_d_object();
      return m_traits.contained_in_affine_hull_d_object() (
        boost::make_transform_iterator(start, cp),
        boost::make_transform_iterator(end, cp),
        cp(p)
      );
    }
  };

  //===========================================================================

  // Required by TriangulationTraits
  class Compare_lexicographically_d
  {
    const RTTraits &m_traits;

  public:
    typedef Comparison_result result_type;
    
    Compare_lexicographically_d(const RTTraits &kernel)
      : m_traits(kernel) {}

    result_type operator()(
      const Weighted_point_d & p, const Weighted_point_d & q) const
    {
      Construct_point_d cp = m_traits.construct_point_d_object();
      return m_traits.compare_lexicographically_d_object()(cp(p), cp(q));
    }
  };
  
  //===========================================================================

  // Only for Triangulation_off_ostream.h (undocumented)
  class Compute_coordinate_d
  {
    const RTTraits &m_traits;

  public:
    typedef FT result_type;
    
    Compute_coordinate_d(const RTTraits &kernel)
      : m_traits(kernel) {}

    result_type operator()(
      const Weighted_point_d & p, const int i) const
    {
      Construct_point_d cp = m_traits.construct_point_d_object();
      return m_traits.compute_coordinate_d_object()(cp(p), i);
    }
  };

  //===========================================================================

  // To satisfy SpatialSortingTraits_d
  // and also for Triangulation_off_ostream.h (undocumented)
  class Point_dimension_d
  {
    const RTTraits &m_traits;

  public:
    typedef int result_type;
    
    Point_dimension_d(const RTTraits &kernel)
      : m_traits(kernel) {}

    result_type operator()(
      const Weighted_point_d & p) const
    {
      Construct_point_d cp = m_traits.construct_point_d_object();
      return m_traits.point_dimension_d_object()(cp(p));
    }
  };
  
  //===========================================================================
  // Object creation
  //===========================================================================

  Less_coordinate_d less_coordinate_d_object() const
  {
    return Less_coordinate_d(*this);
  }
  Contained_in_affine_hull_d contained_in_affine_hull_d_object() const
  { 
    return Contained_in_affine_hull_d(*this); 
  }
  Orientation_d orientation_d_object() const
  {
    return Orientation_d(*this); 
  }
  Construct_flat_orientation_d construct_flat_orientation_d_object() const
  { 
    return Construct_flat_orientation_d(*this);
  }
  In_flat_orientation_d in_flat_orientation_d_object() const
  { 
    return In_flat_orientation_d(*this);
  }
  Compare_lexicographically_d compare_lexicographically_d_object() const
  { 
    return Compare_lexicographically_d(*this);
  }
  Compute_coordinate_d compute_coordinate_d_object() const
  { 
    return Compute_coordinate_d(*this);
  }
  Point_dimension_d point_dimension_d_object() const
  { 
    return Point_dimension_d(*this);
  }
};


} //namespace CGAL

#endif // CGAL_REGULAR_TRIANGULATION_TRAITS_ADAPTER_H

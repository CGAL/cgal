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
//
// Author(s)     : Clement Jamin

#ifndef CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_H
#define CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/representation_tags.h>
#include <CGAL/Kernel_traits.h>

#include <boost/iterator/transform_iterator.hpp>

namespace CGAL {

template < class K, class Weight = typename K::RT >
class Regular_triangulation_euclidean_traits
  : public K
{
public:
  typedef K                                                 Base;

  // Types from K
  typedef typename K::Dimension                     Dimension;
  typedef typename K::FT                            FT;
  typedef typename K::Point_d                       Bare_point;
  typedef typename K::Weighted_point_d              Weighted_point;
  typedef Weighted_point                            Point_d;

  typedef typename K::Construct_weighted_point_d    Construct_weighted_point_d;
  typedef typename K::Power_test_d                  Power_test_d;
  typedef typename K::In_flat_power_test_d          In_flat_power_test_d;
  typedef typename K::Flat_orientation_d            Flat_orientation_d;
  typedef typename K::Point_drop_weight_d           Point_drop_weight_d;
  typedef typename K::Point_weight_d                Point_weight_d;
  
  //=============================================================================
  // Custom types
  //=============================================================================
  
  class Orientation_d
  {
    const K &m_kernel;

  public:
    typedef Orientation result_type;

    Orientation_d(const K &kernel)
      : m_kernel(kernel) {}

    template <typename ForwardIterator> 
    result_type operator()(ForwardIterator start, ForwardIterator end) const
    {
      Point_drop_weight_d pdw = m_kernel.point_drop_weight_d_object();
      return m_kernel.orientation_d_object() (
        boost::make_transform_iterator(start, pdw),
        boost::make_transform_iterator(end, pdw)
      );
    }
  };

  //=============================================================================

  class Construct_flat_orientation_d
  {
    const K &m_kernel;

  public:
    typedef Flat_orientation_d result_type;
    
    Construct_flat_orientation_d(const K &kernel)
      : m_kernel(kernel) {}

    template <typename ForwardIterator> 
    result_type operator()(ForwardIterator start, ForwardIterator end) const
    {
      Point_drop_weight_d pdw = m_kernel.point_drop_weight_d_object();
      return m_kernel.construct_flat_orientation_d_object() (
        boost::make_transform_iterator(start, pdw),
        boost::make_transform_iterator(end, pdw)
      );
    }
  };


  //=============================================================================

  class In_flat_orientation_d
  {
    const K &m_kernel;

  public:
    typedef Orientation result_type;
    
    In_flat_orientation_d(const K &kernel)
      : m_kernel(kernel) {}

    template <typename ForwardIterator> 
    result_type operator()(Flat_orientation_d orient, 
      ForwardIterator start, ForwardIterator end) const
    {
      Point_drop_weight_d pdw = m_kernel.point_drop_weight_d_object();
      return m_kernel.in_flat_orientation_d_object() (
        orient,
        boost::make_transform_iterator(start, pdw),
        boost::make_transform_iterator(end, pdw)
      );
    }
  };

  //=============================================================================

  class Contained_in_affine_hull_d
  {
    const K &m_kernel;

  public:
    typedef bool result_type;
    
    Contained_in_affine_hull_d(const K &kernel)
      : m_kernel(kernel) {}

    template <typename ForwardIterator> 
    result_type operator()(ForwardIterator start, ForwardIterator end, 
                           const Weighted_point & p) const
    {
      Point_drop_weight_d pdw = m_kernel.point_drop_weight_d_object();
      return m_kernel.contained_in_affine_hull_d_object() (
        boost::make_transform_iterator(start, pdw),
        boost::make_transform_iterator(end, pdw),
        pdw(p)
      );
    }
  };

  //=============================================================================

  class Compare_lexicographically_d
  {
    const K &m_kernel;

  public:
    typedef Comparison_result result_type;
    
    Compare_lexicographically_d(const K &kernel)
      : m_kernel(kernel) {}

    result_type operator()(
      const Weighted_point & p, const Weighted_point & q) const
    {
      Point_drop_weight_d pdw = m_kernel.point_drop_weight_d_object();
      return m_kernel.compare_lexicographically_d_object()(pdw(p), pdw(q));
    }
  };
  
  //=============================================================================

  class Compute_coordinate_d
  {
    const K &m_kernel;

  public:
    typedef FT result_type;
    
    Compute_coordinate_d(const K &kernel)
      : m_kernel(kernel) {}

    result_type operator()(
      const Weighted_point & p, const int i) const
    {
      Point_drop_weight_d pdw = m_kernel.point_drop_weight_d_object();
      m_kernel.compute_coordinate_d_object()(pdw(p), i);
      return m_kernel.compute_coordinate_d_object()(pdw(p), i);
    }
  };

  //=============================================================================

  class Point_dimension_d
  {
    const K &m_kernel;

  public:
    typedef int result_type;
    
    Point_dimension_d(const K &kernel)
      : m_kernel(kernel) {}

    result_type operator()(
      const Weighted_point & p) const
    {
      Point_drop_weight_d pdw = m_kernel.point_drop_weight_d_object();
      return m_kernel.point_dimension_d_object()(pdw(p));
    }
  };
  
  //=============================================================================
  // Object creation
  //=============================================================================

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

#endif // CGAL_REGULAR_TRIANGULATION_EUCLIDEAN_TRAITS_H

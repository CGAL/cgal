// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
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
// Author(s)     : 
//

#ifndef CGAL_STRUCTURE_POINT_SET_3_H
#define CGAL_STRUCTURE_POINT_SET_3_H

#include <CGAL/property_map.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/assertions.h>

#include <iterator>
#include <list>


namespace CGAL {


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
/// \cond SKIP_IN_MANUAL
namespace internal {

  template <typename Traits>
  class Point_set_structuring
  {
  public:

    typedef Point_set_structuring<Traits> Self;

    typedef typename Traits::FT FT;
    typedef typename Traits::Point_3 Point;
    typedef typename Traits::Vector_3 Vector;
    typedef typename Traits::Line_3 Line;

    typedef typename Traits::Plane_3 Plane;

    typedef typename Traits::Point_map Point_map;
    typedef typename Traits::Normal_map Normal_map;
    typedef typename Traits::Input_range Input_range;

    typedef typename Input_range::iterator Input_iterator;

    typedef Shape_detection_3::Shape_base<Traits> Shape; 
    typedef Shape_detection_3::Plane<Traits> Plane_shape;

  private:

    Traits m_traits;

    Input_iterator m_input_begin;
    Input_iterator m_input_end;
    Point_map m_point_pmap;
    Normal_map m_normal_pmap;
    
    std::vector<boost::shared_ptr<Plane_shape> > m_planes;
    
  public:

    Point_set_structuring (Traits t = Traits ())
      : m_traits (t)
    {

    }

    
    Point_set_structuring (Input_iterator begin, Input_iterator end,
                           const Shape_detection_3::Efficient_RANSAC<Traits>& shape_detection)
      : m_traits (shape_detection.traits())
    {
      m_input_begin = begin;
      m_input_end = end;

      BOOST_FOREACH (boost::shared_ptr<Shape> shape, shape_detection.shapes())
        {
          boost::shared_ptr<Plane_shape> pshape
            = boost::dynamic_pointer_cast<Plane_shape>(shape);
        
          // Ignore all shapes other than plane
          if (pshape == boost::shared_ptr<Plane_shape>())
            continue;
          m_planes.push_back (pshape);
        }

    }

    
    virtual ~Point_set_structuring ()
    {
      clear ();
    }

    void clear ()
    {

    }

    void run (double radius)
    {

    }
    
  };
  
} /* namespace internal */
/// \endcond



// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------

/// \ingroup PkgPointSetProcessing
/// TODO documentation

// This variant requires the kernel.
template <typename InputIterator,
          typename PointPMap,
          typename EfficientRANSACTraits,
          typename Kernel
>
void
structure_point_set (InputIterator first,  ///< iterator over the first input point.
                     InputIterator beyond, ///< past-the-end iterator over the input points.
                     PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3
                     Shape_detection_3::Efficient_RANSAC<EfficientRANSACTraits>&
                     shape_detection, ///< shape detection engine
                     double radius, ///< attraction radius
                     const Kernel& /*kernel*/) ///< geometric traits.
{
  internal::Point_set_structuring<EfficientRANSACTraits> pss
    (first, beyond, shape_detection);
  pss.run (radius);
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename InputIterator,
          typename PointPMap,
          typename EfficientRANSACTraits
>
void
structure_point_set (InputIterator first,    ///< iterator over the first input point.
                     InputIterator beyond,   ///< past-the-end iterator over the input points.
                     PointPMap point_pmap, ///< property map: value_type of InputIterator -> Point_3
                     Shape_detection_3::Efficient_RANSAC<EfficientRANSACTraits>&
                     shape_detection, ///< shape detection engine
                     double radius) ///< attraction radius
{
  typedef typename boost::property_traits<PointPMap>::value_type Point;
  typedef typename Kernel_traits<Point>::Kernel Kernel;
  return structure_point_set (
    first,beyond,
    point_pmap,
    shape_detection,
    radius,
    Kernel());
}
/// @endcond

/// @cond SKIP_IN_MANUAL
// This variant creates a default point property map = Identity_property_map.
template < typename InputIterator, typename EfficientRANSACTraits >
void
structure_point_set (InputIterator first,    ///< iterator over the first input point.
                     InputIterator beyond,   ///< past-the-end iterator over the input points.
                     Shape_detection_3::Efficient_RANSAC<EfficientRANSACTraits>&
                     shape_detection, ///< shape detection engine
                     double radius) ///< attraction radius
{
  return structure_point_set (
    first,beyond,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
    make_dereference_property_map(first),
#else
    make_identity_property_map(
    typename std::iterator_traits<InputIterator>::value_type()),
#endif
    shape_detection,
    radius);
}
/// @endcond


} //namespace CGAL

#endif // CGAL_STRUCTURE_POINT_SET_3_H


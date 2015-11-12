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

#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_d.h>

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

    class My_point_property_map{
      const std::vector<Point>& points;
    public:
      typedef Point value_type;
      typedef const value_type& reference;
      typedef std::size_t key_type;
      typedef boost::lvalue_property_map_tag category;  
      My_point_property_map (const std::vector<Point>& pts) : points (pts) {}
      reference operator[] (key_type k) const { return points[k]; }
      friend inline reference get (const My_point_property_map& ppmap, key_type i) 
      { return ppmap[i]; }
    };
      

    Traits m_traits;

    std::vector<Point> m_points;
    std::vector<int> m_indices;
    Point_map m_point_pmap;
    Normal_map m_normal_pmap;
    
    std::vector<boost::shared_ptr<Plane_shape> > m_planes;
    std::vector<std::pair<std::size_t, std::size_t> > m_pairs;
    
  public:

    Point_set_structuring (Traits t = Traits ())
      : m_traits (t)
    {

    }

    
    Point_set_structuring (Input_iterator begin, Input_iterator end,
                           const Shape_detection_3::Efficient_RANSAC<Traits>& shape_detection)
      : m_traits (shape_detection.traits())
    {

      for (Input_iterator it = begin; it != end; ++ it)
        m_points.push_back (get(m_point_pmap, *it));

      m_indices = std::vector<int> (m_points.size (), -1);

      BOOST_FOREACH (boost::shared_ptr<Shape> shape, shape_detection.shapes())
        {
          boost::shared_ptr<Plane_shape> pshape
            = boost::dynamic_pointer_cast<Plane_shape>(shape);
        
          // Ignore all shapes other than plane
          if (pshape == boost::shared_ptr<Plane_shape>())
            continue;
          m_planes.push_back (pshape);

          for (std::size_t i = 0; i < pshape->indices_of_assigned_points().size (); ++ i)
            m_indices[pshape->indices_of_assigned_points()[i]] = m_planes.size () - 1;
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

      std::cerr << "Finding adjacent primitives... " << std::endl;
      find_pairs_of_adjacent_primitives (radius);
      std::cerr << "Found " << m_pairs.size () << " pair(s) of adjacent primitives." << std::endl;

    }

  private:

    void find_pairs_of_adjacent_primitives (double radius)
    {
      typedef typename Traits::Search_traits Search_traits_base;
      typedef Search_traits_adapter <std::size_t, My_point_property_map, Search_traits_base> Search_traits;
      typedef CGAL::Kd_tree<Search_traits> Tree;
      typedef CGAL::Fuzzy_sphere<Search_traits> Fuzzy_sphere;

      My_point_property_map pmap (m_points);

      Tree tree (boost::counting_iterator<std::size_t> (0),
                 boost::counting_iterator<std::size_t> (m_points.size()),
                 typename Tree::Splitter(),
                 Search_traits (pmap));

      std::vector<std::vector<bool> > adjacency_table (m_planes.size (),
                                                       std::vector<bool> (m_planes.size (), false));

      //compute a basic adjacency relation (two primitives are neighbors
      //if at least one point of the primitive 1 is a k-nearest neighbor
      //of a point of the primitive 2 and vice versa)
      for (std::size_t i = 0; i < m_points.size (); ++ i)
        {
          int ind_i = m_indices[i];

          if (ind_i == -1)
            continue;

          Fuzzy_sphere query (i, radius, 0., tree.traits());
          
          std::vector<std::size_t> neighbors;
          tree.search (std::back_inserter (neighbors), query); // WIP: SegFaults so far...

          
          for (std::size_t k = 0; k < neighbors.size(); ++ k)
            {
              int ind_k = m_indices[neighbors[k]];
              if (ind_k != -1 && ind_k != ind_i)
                adjacency_table[ind_i][ind_k] = true;
            }
        }

      //verify the symmetry and store the pairs of primitives in
      //primitive_pairs
      for (std::size_t i = 0; i < adjacency_table.size() - 1; ++ i)
        for (std::size_t j = i + 1; j < adjacency_table[i].size(); ++ j)
          if ((adjacency_table[i][j]) && (adjacency_table[j][i]))
            m_pairs.push_back (std::make_pair (i, j));

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


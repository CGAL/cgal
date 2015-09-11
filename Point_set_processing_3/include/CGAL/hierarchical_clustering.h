// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
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
// $URL:
// $Id: 
//
// Author(s) : Simon Giraudot, Pierre Alliez


#ifndef HIERARCHICAL_CLUSTERING_H
#define HIERARCHICAL_CLUSTERING_H

#include <cmath>
#include <stack>

#include <CGAL/property_map.h>
#include <CGAL/basic.h>
#include <CGAL/Dimension.h>
#include <CGAL/Object.h>
#include <CGAL/centroid.h>
#include <CGAL/point_set_processing_assertions.h>
#include <CGAL/Default_diagonalize_traits.h>

namespace CGAL {


  /// \ingroup PkgPointSetProcessing
  
/// Recursively split the point set in smaller clusters until the
/// clusters have less than `size` elements or until their variation
/// factor is below `var_max`.
///
/// This method does not change the input point set: the output is not
/// a subset of the input and is stored in a different container.
///
/// \pre `1/3 > var_max > 0`
/// \pre 'size > 0`
///
/// @tparam InputIterator iterator over input points.
/// @tparam PointPMap is a model of `ReadablePropertyMap` with value type `Point_3<Kernel>`.
///        It can be omitted if the value type of `InputIterator` is convertible to `Point_3<Kernel>`.
/// @tparam OuputIterator back inserter on a container with value type `Point_3<Kernel>`.
/// @tparam DiagonalizeTraits is a model of `DiagonalizeTraits`. It
///        can be omitted: if Eigen 3 (or greater) is available and
///        `CGAL_EIGEN3_ENABLED` is defined then an overload using
///        `Eigen_diagonalize_traits` is provided. Otherwise, the internal
///        implementation `Internal_diagonalize_traits` is used.
/// @tparam Kernel Geometric traits class.
///        It can be omitted and deduced automatically from the value type of `PointPMap`.
///

// This variant requires all parameters.

  template <typename InputIterator,
	    typename PointPMap,
	    typename OutputIterator,
	    typename DiagonalizeTraits,
	    typename Kernel>
  void hierarchical_clustering (InputIterator begin,
				InputIterator end,
				PointPMap point_pmap,
				OutputIterator out,
				const unsigned int size,
				const double var_max,
				const DiagonalizeTraits&,
				const Kernel&)
  {
    typedef typename Kernel::Plane_3  Plane;
    typedef typename Kernel::Point_3  Point;
    typedef typename Kernel::Vector_3 Vector;

    // We define a cluster as a point set + its centroid (useful for
    // faster computations of centroids - to be implemented)
    typedef std::pair< std::list<Point>, Point > cluster;

    std::list<cluster> clusters_stack;
    typedef typename std::list<cluster>::iterator cluster_iterator;

    CGAL_precondition (begin != end);
    CGAL_point_set_processing_precondition
      (var_max >= 0.0 && var_max <= 1./3.);

    // The first cluster is the whole input point set
    clusters_stack.push_front (cluster (std::list<Point>(), Point (0., 0., 0.)));
    for(InputIterator it = begin; it != end; it++)
      {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
	Point point = get(point_pmap, it);
#else
	Point point = get(point_pmap, *it);
#endif
	clusters_stack.front ().first.push_back (point);
      }
    clusters_stack.front ().second = centroid (clusters_stack.front ().first.begin (),
					       clusters_stack.front ().first.end ());
    
    while (!(clusters_stack.empty ()))
      {
	cluster& current_cluster = clusters_stack.back ();

	// If the cluster only has 1 element, we add it to the list of
	// output points
	if (current_cluster.first.size () == 1)
	  {
	    *(out ++) = current_cluster.second;
	    clusters_stack.pop_back ();
	    continue;
	  }

	// Compute the covariance matrix of the set
	cpp11::array<double, 6> covariance = {{ 0., 0., 0., 0., 0., 0. }};

	for (typename std::list<Point>::iterator it = current_cluster.first.begin ();
	     it != current_cluster.first.end (); ++ it)
	  {
	    const Point& p = *it;
	    Vector d = p - current_cluster.second;
	    covariance[0] += d.x () * d.x ();
	    covariance[1] += d.x () * d.y ();
	    covariance[2] += d.x () * d.z ();
	    covariance[3] += d.y () * d.y ();
	    covariance[4] += d.y () * d.z ();
	    covariance[5] += d.z () * d.z ();
	  }

	cpp11::array<double, 3> eigenvalues = {{ 0., 0., 0. }};
	cpp11::array<double, 9> eigenvectors = {{ 0., 0., 0.,
						  0., 0., 0.,
						  0., 0., 0. }};
	// Linear algebra = get eigenvalues and eigenvectors for
	// PCA-like analysis
	DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
	  (covariance, eigenvalues, eigenvectors);
	
	// Variation of the set defined as lambda_min / (lambda_0 + lambda_1 + lambda_2)
	double var = 0.;
	for (int i = 0; i < 3; ++ i)
	  var += eigenvalues[i];
	var = eigenvalues[0] / var;

	// Split the set if size OR variance of the cluster is too large
	if (current_cluster.first.size () > size || var > var_max)
	  {
	    clusters_stack.push_front (cluster (std::list<Point>(), Point (0., 0., 0.)));
	    cluster_iterator positive_side = clusters_stack.begin ();
	    clusters_stack.push_front (cluster (std::list<Point>(), Point (0., 0., 0.)));
	    cluster_iterator negative_side = clusters_stack.begin ();
	    
	    // Compute the plane which splits the point set into 2 point sets:
	    //  * Normal to the eigenvector with highest eigenvalue
	    //  * Passes through the centroid of the set
	    Plane plane (current_cluster.second, Vector (eigenvectors[6], eigenvectors[7], eigenvectors[8]));

	    std::size_t current_cluster_size = 0;
	    typename std::list<Point>::iterator it = current_cluster.first.begin ();
	    while (it != current_cluster.first.end ())
	      {
		typename std::list<Point>::iterator current = it ++;

		std::list<Point>& side = (plane.has_on_positive_side (*current)
					  ? positive_side->first : negative_side->first);
		side.splice (side.end (), current_cluster.first, current);
		++ current_cluster_size;
	      }

	    if (positive_side->first.empty () || negative_side->first.empty ())
	      {
		cluster_iterator empty, nonempty;
		if (positive_side->first.empty ())
		  {
		    empty = positive_side;
		    nonempty = negative_side;
		  }
		else
		  {
		    empty = negative_side;
		    nonempty = positive_side;
		  }

		nonempty->second = centroid (nonempty->first.begin (), nonempty->first.end ());

		clusters_stack.erase (empty);
	      }
	    else
	      {
		// Compute the centroids
		positive_side->second = centroid (positive_side->first.begin (), positive_side->first.end ());

		// The second centroid can be computed with the first and
		// the previous ones :
		// centroid_neg = (n_total * centroid - n_pos * centroid_pos)
		//                 / n_neg;
		negative_side->second = Point ((current_cluster_size * current_cluster.second.x ()
						- positive_side->first.size () * positive_side->second.x ())
					       / negative_side->first.size (),
					       (current_cluster_size * current_cluster.second.y ()
						- positive_side->first.size () * positive_side->second.y ())
					       / negative_side->first.size (),
					       (current_cluster_size * current_cluster.second.z ()
						- positive_side->first.size () * positive_side->second.z ())
					       / negative_side->first.size ());

	      }
	    clusters_stack.pop_back ();
	  }
	// If the size/variance are small enough, add the centroid as
	// and output point
	else
	  {
	    *(out ++) = current_cluster.second;
	    clusters_stack.pop_back ();
	  }
      }
  }
/// @endcond


/// @cond SKIP_IN_MANUAL
  // This variant deduces the kernel from the iterator type.
  template <typename InputIterator,
	    typename PointPMap,
	    typename OutputIterator,
	    typename DiagonalizeTraits>
  void hierarchical_clustering (InputIterator begin,
				InputIterator end,
				PointPMap point_pmap,
				OutputIterator out,
				const unsigned int size,
				const double var_max,
				const DiagonalizeTraits& diagonalize_traits)
  {
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    hierarchical_clustering (begin, end, point_pmap, out, size, var_max,
			     diagonalize_traits, Kernel());
  }
/// @endcond

/// @cond SKIP_IN_MANUAL  
  // This variant uses default diagonalize traits
  template <typename InputIterator,
	    typename PointPMap,
	    typename OutputIterator>
  void hierarchical_clustering (InputIterator begin,
				InputIterator end,
				PointPMap point_pmap,
				OutputIterator out,
				const unsigned int size,
				const double var_max)
  {
    typedef typename boost::property_traits<PointPMap>::value_type Point;
    typedef typename Kernel_traits<Point>::Kernel Kernel;
    hierarchical_clustering (begin, end, point_pmap, out, size, var_max,
			     Default_diagonalize_traits<double, 3> (), Kernel());
  }
/// @endcond  

/// @cond SKIP_IN_MANUAL
  // This variant creates a default point property map = Identity_property_map.
  template <typename InputIterator,
	    typename OutputIterator>
  void hierarchical_clustering (InputIterator begin,
				InputIterator end,
				OutputIterator out,
				const unsigned int size = 10,
				const double var_max = 0.333)
  {
    hierarchical_clustering
      (begin, end,
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
       make_dereference_property_map(first),
#else
       make_identity_property_map (typename std::iterator_traits<InputIterator>::value_type()),
#endif
       out, size, var_max);
  }
/// @endcond  

} // namespace CGAL

#endif // HIERARCHICAL_CLUSTERING_H

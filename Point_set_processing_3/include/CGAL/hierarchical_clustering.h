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
#include <CGAL/PCA_util.h>

namespace CGAL {


  /// \ingroup PkgPointSetProcessing
  
/// Recursively split the point set in smaller clusters until the
/// clusters have less than `size` elements or until their variation
/// factor is below `var_max`.
///
/// This method does not change the input point set: the output is not
/// a subset of the input and is stored in a different container.
///
/// \pre `0 < var_max < 1/3`
/// \pre `size > 0`
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
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3  Point;
    typedef typename Kernel::Vector_3 Vector;

    // We define a cluster as a point set + its centroid (useful for
    // faster computations of centroids - to be implemented)
    typedef std::pair< std::list<Point>, Point > cluster;

    std::list<cluster> clusters_stack;
    typedef typename std::list<cluster>::iterator cluster_iterator;

    CGAL_precondition (begin != end);
    CGAL_point_set_processing_precondition (size > 0);
    CGAL_point_set_processing_precondition (var_max > 0.0);

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
	cluster_iterator current_cluster = clusters_stack.begin ();

	// If the cluster only has 1 element, we add it to the list of
	// output points
	if (current_cluster->first.size () == 1)
	  {
	    *(out ++) = current_cluster->second;
	    clusters_stack.pop_front ();
	    continue;
	  }

	// Compute the covariance matrix of the set
	cpp11::array<FT, 6> covariance = {{ 0., 0., 0., 0., 0., 0. }};

	internal::assemble_covariance_matrix_3 (current_cluster->first.begin (),
						current_cluster->first.end (),
						covariance,
						current_cluster->second, Kernel(),
						(Point*)NULL, CGAL::Dimension_tag<0>());
	// for (typename std::list<Point>::iterator it = current_cluster->first.begin ();
	//      it != current_cluster->first.end (); ++ it)
	//   {
	//     const Point& p = *it;
	//     Vector d = p - current_cluster->second;
	//     covariance[0] += d.x () * d.x ();
	//     covariance[1] += d.x () * d.y ();
	//     covariance[2] += d.x () * d.z ();
	//     covariance[3] += d.y () * d.y ();
	//     covariance[4] += d.y () * d.z ();
	//     covariance[5] += d.z () * d.z ();
	//   }

	cpp11::array<FT, 3> eigenvalues = {{ 0., 0., 0. }};
	cpp11::array<FT, 9> eigenvectors = {{ 0., 0., 0.,
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
	if (current_cluster->first.size () > size || var > var_max)
	  {
	    clusters_stack.push_front (cluster (std::list<Point>(), Point (0., 0., 0.)));
	    cluster_iterator negative_side = clusters_stack.begin ();
	    // positive_side is built directly from current_cluster
	    
	    // The plane which splits the point set into 2 point sets:
	    //  * Normal to the eigenvector with highest eigenvalue
	    //  * Passes through the centroid of the set
	    Vector v (eigenvectors[6], eigenvectors[7], eigenvectors[8]);

	    std::size_t current_cluster_size = 0;
	    typename std::list<Point>::iterator it = current_cluster->first.begin ();
	    while (it != current_cluster->first.end ())
	      {
		typename std::list<Point>::iterator current = it ++;

		// Test if point is on negative side of plane and
		// transfer it to the negative_side cluster if it is
		if (Vector (current_cluster->second, *current) * v < 0)
		  negative_side->first.splice (negative_side->first.end (),
					       current_cluster->first, current);
		++ current_cluster_size;
	      }

	    // If one of the clusters is empty, only keep the non-empty one
	    if (current_cluster->first.empty () || negative_side->first.empty ())
	      {
		cluster_iterator empty, nonempty;
		if (current_cluster->first.empty ())
		  {
		    empty = current_cluster;
		    nonempty = negative_side;
		  }
		else
		  {
		    empty = negative_side;
		    nonempty = current_cluster;
		  }

		nonempty->second = centroid (nonempty->first.begin (), nonempty->first.end ());

		clusters_stack.erase (empty);
	      }
	    else
	      {
		// Save old centroid for faster computation
		Point old_centroid = current_cluster->second;
		
		// Compute the first centroid
		current_cluster->second = centroid (current_cluster->first.begin (), current_cluster->first.end ());

		// The second centroid can be computed with the first and
		// the old ones :
		// centroid_neg = (n_total * old_centroid - n_pos * first_centroid)
		//                 / n_neg;
		negative_side->second = Point ((current_cluster_size * old_centroid.x ()
						- current_cluster->first.size () * current_cluster->second.x ())
					       / negative_side->first.size (),
					       (current_cluster_size * old_centroid.y ()
						- current_cluster->first.size () * current_cluster->second.y ())
					       / negative_side->first.size (),
					       (current_cluster_size * old_centroid.z ()
						- current_cluster->first.size () * current_cluster->second.z ())
					       / negative_side->first.size ());
	      }
	  }
	// If the size/variance are small enough, add the centroid as
	// and output point
	else
	  {
	    *(out ++) = current_cluster->second;
	    clusters_stack.pop_front ();
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

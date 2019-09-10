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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s) : Simon Giraudot, Pierre Alliez


#ifndef HIERARCHY_SIMPLIFY_POINT_SET_H
#define HIERARCHY_SIMPLIFY_POINT_SET_H

#include <CGAL/license/Point_set_processing_3.h>

#include <CGAL/disable_warnings.h>

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
#include <CGAL/squared_distance_3.h>
#include <CGAL/Iterator_range.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

namespace CGAL {


  /// \cond SKIP_IN_MANUAL
  namespace internal {

    template < typename InputIterator,
	       typename PointMap,
	       typename K >
    typename K::Point_3
    hsps_centroid(InputIterator begin, 
		 InputIterator end,
		 PointMap& point_map,	     
		 const K&)
    {
      typedef typename K::Point_3 Point;
      typedef typename K::FT FT;

      CGAL_precondition(begin != end);

      FT x = (FT)0., y = (FT)0., z = (FT)0.;
      unsigned int nb_pts = 0;
      while(begin != end) 
	{
	  typename boost::property_traits<PointMap>::reference point =
            get(point_map, *begin);
	  x += point.x ();  y += point.y ();  z += point.z ();
	  ++ nb_pts;
	  ++ begin;
	}
      return Point (x/nb_pts, y/nb_pts, z/nb_pts);
    }

    template < typename Input_type,
	       typename PointMap,
	       typename K >
    void
    hsc_terminate_cluster (std::list<Input_type>& cluster,
			   std::list<Input_type>& points_to_keep,
			   std::list<Input_type>& points_to_remove,
			   PointMap& point_map,
			   const typename K::Point_3& centroid,
			   const K&)
    {
      typedef typename std::list<Input_type>::iterator Iterator;
      typedef typename K::FT FT;

      FT dist_min = (std::numeric_limits<FT>::max)();

      typename std::list<Input_type>::iterator point_min;
      for (Iterator it = cluster.begin (); it != cluster.end (); ++ it)
	{
	  FT dist = CGAL::squared_distance (get(point_map, *it), centroid);
	  if (dist < dist_min)
	    {
	      dist_min = dist;
	      point_min = it;
	    }
	}

      points_to_keep.splice (points_to_keep.end (), cluster, point_min);
      points_to_remove.splice (points_to_remove.end (), cluster, cluster.begin (), cluster.end ());
    }			       
			    



  } // namespace internal
  /// \endcond

  /**
     \ingroup PkgPointSetProcessingAlgorithms
  
     Recursively split the point set in smaller clusters until the
     clusters have less than `size` elements or until their variation
     factor is below `var_max`. 
  
     This method modifies the order of input points so as to pack all remaining points first,
     and returns an iterator over the first point to remove (see erase-remove idiom).
     For this reason it should not be called on sorted containers.
  
     \pre `0 < maximum_variation < 1/3`
     \pre `size > 0`
  
     \tparam PointRange is a model of `Range`. The value type of
     its iterator is the key type of the named parameter `point_map`.

     \param points input point range.
     \param np optional sequence of \ref psp_namedparameters "Named Parameters" among the ones listed below.

     \cgalNamedParamsBegin
       \cgalParamBegin{point_map} a model of `ReadWritePropertyMap` with value type `geom_traits::Point_3`.
       If this parameter is omitted, `CGAL::Identity_property_map<geom_traits::Point_3>` is used.\cgalParamEnd
       \cgalParamBegin{size} maximum cluster size.\cgalParamEnd
       \cgalParamBegin{maximum_variation} maximum cluster variation value.\cgalParamEnd
       \cgalParamBegin{diagonalize_traits} a model of `DiagonalizeTraits`. It can be omitted:
       if Eigen 3 (or greater) is available and `CGAL_EIGEN3_ENABLED` is defined then an overload
       using `Eigen_diagonalize_traits` is provided. Otherwise, the internal implementation
       `CGAL::Diagonalize_traits` is used.\cgalParamEnd
       \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `Kernel`\cgalParamEnd
     \cgalNamedParamsEnd

     \return iterator over the first point to remove.
  */
  template <typename PointRange,
	    typename NamedParameters>
  typename PointRange::iterator
  hierarchy_simplify_point_set (PointRange& points,
                                const NamedParameters& np)
  {
    using boost::choose_param;
  
    // basic geometric types
    typedef typename Point_set_processing_3::GetPointMap<PointRange, NamedParameters>::type PointMap;
    typedef typename Point_set_processing_3::GetK<PointRange, NamedParameters>::Kernel Kernel;
    typedef typename GetDiagonalizeTraits<NamedParameters, double, 3>::type DiagonalizeTraits;

    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    typedef typename Kernel::FT FT;

    PointMap point_map = choose_param(get_param(np, internal_np::point_map), PointMap());
    unsigned int size = choose_param(get_param(np, internal_np::size), 10);
    double var_max = choose_param(get_param(np, internal_np::maximum_variation), 1./3.);

    typedef typename std::iterator_traits<typename PointRange::iterator>::value_type Input_type;

    // We define a cluster as a point set + its centroid (useful for
    // faster computations of centroids - to be implemented)
    typedef std::pair< std::list<Input_type>, Point > cluster;

    std::list<cluster> clusters_stack;
    typedef typename std::list<cluster>::iterator cluster_iterator;

    CGAL_precondition (points.begin() != points.end());
    CGAL_point_set_processing_precondition (size > 0);
    CGAL_point_set_processing_precondition (var_max > 0.0);

    // The first cluster is the whole input point set
    clusters_stack.push_front (cluster (std::list<Input_type>(), Point (0., 0., 0.)));
    std::copy (points.begin(), points.end(), std::back_inserter (clusters_stack.front ().first));
    
    clusters_stack.front ().second = internal::hsps_centroid (clusters_stack.front ().first.begin (),
                                                              clusters_stack.front ().first.end (),
                                                              point_map, Kernel());

    std::list<Input_type> points_to_keep;
    std::list<Input_type> points_to_remove;
    
    while (!(clusters_stack.empty ()))
      {
	cluster_iterator current_cluster = clusters_stack.begin ();

	// If the cluster only has 1 element, we add it to the list of
	// output points
	if (current_cluster->first.size () == 1)
	  {
	    points_to_keep.splice (points_to_keep.end (), current_cluster->first,
				   current_cluster->first.begin ());
	    clusters_stack.pop_front ();
	    continue;
	  }

	// Compute the covariance matrix of the set
	cpp11::array<FT, 6> covariance = {{ 0., 0., 0., 0., 0., 0. }};

	for (typename std::list<Input_type>::iterator it = current_cluster->first.begin ();
	     it != current_cluster->first.end (); ++ it)
	  {
	    const Point& point = get(point_map, *it);
	    Vector d = point - current_cluster->second;
	    covariance[0] += d.x () * d.x ();
	    covariance[1] += d.x () * d.y ();
	    covariance[2] += d.x () * d.z ();
	    covariance[3] += d.y () * d.y ();
	    covariance[4] += d.y () * d.z ();
	    covariance[5] += d.z () * d.z ();
	  }

	cpp11::array<FT, 3> eigenvalues = {{ 0., 0., 0. }};
	cpp11::array<FT, 9> eigenvectors = {{ 0., 0., 0.,
					      0., 0., 0.,
					      0., 0., 0. }};
	// Linear algebra = get eigenvalues and eigenvectors for
	// PCA-like analysis
	DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
	  (covariance, eigenvalues, eigenvectors);
	
	// Variation of the set defined as lambda_min / (lambda_0 + lambda_1 + lambda_2)
	double var = eigenvalues[0] / (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);

	// Split the set if size OR variance of the cluster is too large
	if (current_cluster->first.size () > size || var > var_max)
	  {
	    clusters_stack.push_front (cluster (std::list<Input_type>(), Point (0., 0., 0.)));
	    cluster_iterator negative_side = clusters_stack.begin ();
	    // positive_side is built directly from current_cluster
	    
	    // The plane which splits the point set into 2 point sets:
	    //  * Normal to the eigenvector with highest eigenvalue
	    //  * Passes through the centroid of the set
	    Vector v (eigenvectors[6], eigenvectors[7], eigenvectors[8]);

	    std::size_t current_cluster_size = 0;
	    typename std::list<Input_type>::iterator it = current_cluster->first.begin ();
	    while (it != current_cluster->first.end ())
	      {
		typename std::list<Input_type>::iterator current = it ++;
                const Point& point = get(point_map, *current);

		// Test if point is on negative side of plane and
		// transfer it to the negative_side cluster if it is
		if (Vector (current_cluster->second, point) * v < 0)
		  negative_side->first.splice (negative_side->first.end (),
					       current_cluster->first, current);
		++ current_cluster_size;
	      }

	    // If one of the clusters is empty, stop to avoid infinite
	    // loop and keep the non-empty one
	    if (current_cluster->first.empty () || negative_side->first.empty ())
	      {
		cluster_iterator nonempty = (current_cluster->first.empty ()
					     ? negative_side : current_cluster);

		// Compute the centroid
		nonempty->second = internal::hsps_centroid (nonempty->first.begin (),
							    nonempty->first.end (),
							    point_map, Kernel());
		
		internal::hsc_terminate_cluster (nonempty->first,
						 points_to_keep,
						 points_to_remove,
						 point_map,
						 nonempty->second,
						 Kernel ());
		clusters_stack.pop_front ();
		clusters_stack.pop_front ();
	      }
	    else
	      {
		// Save old centroid for faster computation
		Point old_centroid = current_cluster->second;
		
		// Compute the first centroid
		current_cluster->second = internal::hsps_centroid (current_cluster->first.begin (),
								   current_cluster->first.end (),
								   point_map, Kernel());

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
	    internal::hsc_terminate_cluster (current_cluster->first,
					     points_to_keep,
					     points_to_remove,
					     point_map,
					     current_cluster->second,
					     Kernel ());
	    clusters_stack.pop_front ();
	  }
      }
    typename PointRange::iterator first_point_to_remove =
      std::copy (points_to_keep.begin(), points_to_keep.end(), points.begin());
    std::copy (points_to_remove.begin(), points_to_remove.end(), first_point_to_remove);

    return first_point_to_remove;

  }


  /// \cond SKIP_IN_MANUAL
  // variant with default NP
  template <typename PointRange>
  typename PointRange::iterator
  hierarchy_simplify_point_set (PointRange& points)
  {
    return hierarchy_simplify_point_set
      (points, CGAL::Point_set_processing_3::parameters::all_default(points));
  }

#ifndef CGAL_NO_DEPRECATED_CODE
  // deprecated API
  template <typename ForwardIterator,
	    typename PointMap,
	    typename DiagonalizeTraits,
	    typename Kernel>
  CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::hierarchy_simplify_point_set(), please update your code")
  ForwardIterator hierarchy_simplify_point_set (ForwardIterator begin,
						ForwardIterator end,
						PointMap point_map,
						const unsigned int size,
						const double var_max,
						const DiagonalizeTraits&,
						const Kernel&)
  {
    CGAL::Iterator_range<ForwardIterator> points (begin, end);
    return hierarchy_simplify_point_set
      (points,
       CGAL::parameters::point_map (point_map).
       size (size).
       maximum_variation (var_max).
       diagonalize_traits (DiagonalizeTraits()).
       geom_traits(Kernel()));
  }

  // deprecated API
  template <typename ForwardIterator,
	    typename PointMap,
	    typename DiagonalizeTraits>
  CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::hierarchy_simplify_point_set(), please update your code")
  ForwardIterator hierarchy_simplify_point_set (ForwardIterator begin,
						ForwardIterator end,
						PointMap point_map,
						const unsigned int size,
						const double var_max,
						const DiagonalizeTraits& diagonalize_traits)
  {
    CGAL::Iterator_range<ForwardIterator> points (begin, end);
    return hierarchy_simplify_point_set
      (points,
       CGAL::parameters::point_map (point_map).
       size (size).
       maximum_variation (var_max).
       diagonalize_traits (diagonalize_traits));
  }

  // deprecated API
  template <typename ForwardIterator,
	    typename PointMap >
  CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::hierarchy_simplify_point_set(), please update your code")
  ForwardIterator hierarchy_simplify_point_set (ForwardIterator begin,
						ForwardIterator end,
						PointMap point_map,
						const unsigned int size,
						const double var_max)
  {
    CGAL::Iterator_range<ForwardIterator> points (begin, end);
    return hierarchy_simplify_point_set
      (points,
       CGAL::parameters::point_map (point_map).
       size (size).
       maximum_variation (var_max));
  }
  
  // deprecated API
  template <typename ForwardIterator >
  CGAL_DEPRECATED_MSG("you are using the deprecated V1 API of CGAL::hierarchy_simplify_point_set(), please update your code")
  ForwardIterator hierarchy_simplify_point_set (ForwardIterator begin,
						ForwardIterator end,
						const unsigned int size = 10,
						const double var_max = 0.333)
  {
    CGAL::Iterator_range<ForwardIterator> points (begin, end);
    return hierarchy_simplify_point_set
      (points,
       CGAL::parameters::size (size).
       maximum_variation (var_max));
  }
#endif // CGAL_NO_DEPRECATED_CODE
  /// \endcond  

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // HIERARCHY_SIMPLIFY_POINT_SET_H

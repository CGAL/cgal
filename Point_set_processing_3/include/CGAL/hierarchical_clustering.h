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
#include <queue>

#include <CGAL/basic.h>
#include <CGAL/Dimension.h>
#include <CGAL/Object.h>
#include <CGAL/centroid.h>
#include <CGAL/eigen.h>
#include <CGAL/point_set_processing_assertions.h>


// For testing purposes, both CGAL linear algebra and Eigen library
// are used. You can switch from one to the other by commenting the
// following "define" line.
#define USE_EIGEN_LIB

#ifdef USE_EIGEN_LIB
#include <eigen3/Eigen/Dense>
#endif

#undef min
#undef max

namespace CGAL {



  // TODO: * move to PCA
  //       * adapt to the Eigen library?

  /***************************************************************************/
  template <typename InputIterator, typename K>
  typename K::FT* covariance_matrix_3 (InputIterator first,
				       InputIterator beyond,
				       typename K::Point_3 centroid)
  /***************************************************************************/
  {
    typedef typename K::FT FT;
    FT *out = new FT[6];

    for (int i = 0; i < 6; ++ i)
      out[i] = (FT)(0.0);

    for (InputIterator it = first; it != beyond; ++ it)
      {
	const Point& p = *it;
	Vector d = p - centroid;
	out[0] += d.x () * d.x ();
	out[1] += d.x () * d.y ();
	out[2] += d.y () * d.y ();
	out[3] += d.x () * d.z ();
	out[4] += d.y () * d.z ();
	out[5] += d.z () * d.z ();
      }

    return out;
  }




  /***************************************************************************/
  template <typename InputIterator, typename K>
  typename K::FT* covariance_matrix_3 (InputIterator first,
				       InputIterator beyond)
  /***************************************************************************/
  {
    return covariance_matrix_3<InputIterator, K>
      (first, beyond, centroid (first, beyond));
  }




  /***************************************************************************/
  template <typename InputIterator,
	    typename OutputIterator1,
	    typename OutputIterator2,
	    typename OutputIterator3,
	    typename K>
  void point_set_split_plane_3 (InputIterator first,
				InputIterator beyond,
				typename K::Plane_3& plane,
				OutputIterator1 points_on_positive_side,
				OutputIterator2 points_on_negative_side,
				OutputIterator3 points_on_plane)
  /***************************************************************************/
  {
    for (InputIterator it = first; it != beyond; ++ it)
      if (plane.has_on_positive_side (*it))
	*points_on_positive_side++ = *it;
      else if (plane.has_on_negative_side (*it))
	*points_on_negative_side++ = *it;
      else
	*points_on_plane++ = *it;
  }


  // TODO: move to separate file with exposed API?
  // we may need this kind of functions for future interactive editing of point sets.

  /***************************************************************************/
  template <typename InputIterator,
	    typename OutputIterator1,
	    typename OutputIterator2,
	    typename K>
  void point_set_split_plane_3 (InputIterator first,
				InputIterator beyond,
				typename K::Plane_3& plane,
				OutputIterator1 points_on_positive_side,
				OutputIterator2 points_on_negative_side,
				bool strictly_negative = true) // design to discuss
  /***************************************************************************/
  {
    for (InputIterator it = first; it != beyond; ++ it)
      if (plane.has_on_positive_side (*it))
	*points_on_positive_side++ = *it;
      else if (plane.has_on_negative_side (*it))
	*points_on_negative_side++ = *it;
      else
	{
	  if (strictly_negative)
	    *points_on_positive_side++ = *it;
	  else
	    *points_on_negative_side++ = *it;
	}
  }



  /***************************************************************************/
  template <typename InputIterator, typename K>
  InputIterator hierarchical_clustering (InputIterator first,
					 InputIterator beyond,
					 const unsigned int& size,
					 const typename K::FT& var_max)
  /***************************************************************************/
  {
    typedef typename K::FT       FT;
    typedef typename K::Plane_3  Plane;
    typedef typename K::Point_3  Point;
    typedef typename K::Vector_3 Vector;

    // We define a cluster as a point set + its centroid (useful for
    // faster computations of centroids - to be implemented)
    typedef std::pair< std::list<Point>, Point > cluster;

    CGAL_precondition (first != beyond);
    CGAL_point_set_processing_precondition
      (var_max >= (FT)0.0 && var_max <= (FT)(1./3.));

    // The first cluster is the whole input point set 
    std::list<Point> first_cluster;

    // FIXME: use std::copy instead
    for (InputIterator it = first; it != beyond; ++ it)
      first_cluster.push_back (*it);

    // Initialize a queue of clusters
    std::queue<cluster> clusters_queue;
    clusters_queue.push (cluster (first_cluster, centroid (first, beyond)));

    // The algorithm must return a set of points
    std::list<Point> points_to_keep;

    while ( !clusters_queue.empty () )
      {
	cluster current_cluster = clusters_queue.front ();
	clusters_queue.pop ();

	// If the cluster only has 1 element, we add it to the list of
	// output points
	if (current_cluster.first.size () == 1)
	  {
	    points_to_keep.push_back (current_cluster.second);
	    continue;
	  }

#ifndef USE_EIGEN_LIB

	// Compute the covariance matrix of the set
	FT* covariance
	  = covariance_matrix_3<typename std::list<Point>::iterator, K>
	  (current_cluster.first.begin (), // Beginning of the list
	   current_cluster.first.end(),    // End of the list
	   current_cluster.second);        // Centroid

	// Linear algebra = get eigenvalues and eigenvectors for
	// PCA-like analysis
	FT eigen_vectors[9];
	FT eigen_values[3];
	CGAL::internal::eigen_symmetric<FT>
	  (covariance, 3, eigen_vectors, eigen_values);

	// Variation of the set defined as lambda_2 / (lambda_0 + lambda_1 + lambda_2)
	FT var = (FT)(0.0);
	for (int i = 0; i < 3; ++ i)
	  var += eigen_values[i];
	var = eigen_values[2] / var;

	// Compute the normal to the eigenvector with highest eigenvalue
	Vector normal (eigen_vectors[0], eigen_vectors[1], eigen_vectors[2]);

#else

	// Compute the covariance matrix of the set
	Eigen::Matrix3d covariance;
	for (int i = 0; i < 3; ++ i)
	  for (int j = 0; j < 3; j ++)
	    covariance (i, j) = 0.;

	for (typename std::list<Point>::iterator it = current_cluster.first.begin ();
	     it != current_cluster.first.end (); ++ it)
	  {
	    const Point& p = *it;
	    Vector d = p - current_cluster.second;
	    covariance (0, 0) += d.x () * d.x ();
	    covariance (1, 0) += d.x () * d.y ();
	    covariance (1, 1) += d.y () * d.y ();
	    covariance (2, 0) += d.x () * d.z ();
	    covariance (2, 1) += d.y () * d.z ();
	    covariance (2, 2) += d.z () * d.z ();
	  }

	covariance (0, 1) = covariance (1, 0);
	covariance (0, 2) = covariance (2, 0);
	covariance (1, 2) = covariance (2, 1);

	// Linear algebra = get eigenvalues and eigenvectors for
	// PCA-like analysis
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>
	  eigensolver (covariance);
	if (eigensolver.info () != Eigen::Success)
	  std::abort ();
	Eigen::Vector3d eigen_values = eigensolver.eigenvalues ();
	Eigen::Vector3d eigen_vector = eigensolver.eigenvectors ().col (2);

	// Variation of the set defined as lambda_2 / (lambda_0 + lambda_1 + lambda_2)
	FT var = (FT)(0.0);
	for (int i = 0; i < 3; ++ i)
	  var += eigen_values[i];
	var = eigen_values (0) / var;

	// Compute the normal to the eigenvector with highest eigenvalue
	Vector normal (eigen_vector(0), eigen_vector(1), eigen_vector(2));

#endif


	// Split the set if size OR variance of the cluster is too large
	if (current_cluster.first.size () > size || var > var_max)
	  {
	    std::list<Point> positive_side;
	    std::list<Point> negative_side;

	    // Compute the plane which splits the point set into 2 point sets:
	    //  * Normal to the eigenvector with highest eigenvalue
	    //  * Passes through the centroid of the set
	    Plane plane (current_cluster.second, normal);

	    // Split the point sets along this plane
	    typedef typename std::list<Point>::iterator Iterator;
	    point_set_split_plane_3<Iterator,
	      std::back_insert_iterator<std::list<Point> >,
	      std::back_insert_iterator<std::list<Point> >,
	      K>
	      (current_cluster.first.begin (),
	       current_cluster.first.end (),
	       plane,
	       std::back_inserter (positive_side),
	       std::back_inserter (negative_side), true);

	    // Compute the centroids
	    Point centroid_positive
	      = centroid (positive_side.begin (), positive_side.end ());

	    // The second centroid can be computed with the first and
	    // the previous ones :
	    // centroid_neg = (n_total * centroid - n_pos * centroid_pos)
	    //                 / n_neg;
	    Point centroid_negative
	      ((current_cluster.first.size () * current_cluster.second.x ()
	      	- positive_side.size () * centroid_positive.x ())
	       / negative_side.size (),
	       (current_cluster.first.size () * current_cluster.second.y ()
	      	- positive_side.size () * centroid_positive.y ())
	       / negative_side.size (),
	       (current_cluster.first.size () * current_cluster.second.z ()
	      	- positive_side.size () * centroid_positive.z ())
	       / negative_side.size ());

	    // If the sets are non-empty, add the clusters to the queue
	    if (positive_side.size () != 0)
	      clusters_queue.push (cluster (positive_side, centroid_positive));
	    if (negative_side.size () != 0)
	      clusters_queue.push (cluster (negative_side, centroid_negative));
	  }
	// If the size/variance are small enough, add the centroid as
	// and output point
	else
	  points_to_keep.push_back (current_cluster.second);
      }

    // The output of the function is to the first point to be removed
    InputIterator first_point_to_remove =
      std::copy (points_to_keep.begin (), points_to_keep.end (), first);

    return first_point_to_remove;
  }

} // namespace CGAL

#endif // HIERARCHICAL_CLUSTERING_H

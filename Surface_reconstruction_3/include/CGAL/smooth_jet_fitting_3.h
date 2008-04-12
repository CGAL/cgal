// Copyright (c) 2007-08  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://scm.gforge.inria.fr/svn/cgal/trunk/Surface_reconstruction_3/include/CGAL/estimate_normals_jet_fitting_3.h $
// $Id: estimate_normals_jet_fitting_3.h 42587 2008-03-26 15:44:54Z lsaboret $
//
// Author(s) : Pierre Alliez and Marc Pouget

#ifndef CGAL_SMOOTH_JET_FITTING_3_H
#define CGAL_SMOOTH_JET_FITTING_3_H

#include <CGAL/basic.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Oriented_normal_3.h>

#include <iterator>
#include <list>

CGAL_BEGIN_NAMESPACE

/// Smooth one point position using jet fitting on the K
/// nearest neighbors and reprojection onto the jet.
///
/// Precondition: K >= 2.
///
/// @return point
template < typename Kernel, ///< Geometric traits class.
           typename Tree,
           typename Point>
Point
smooth_jet_fitting_3(const typename Kernel::Point_3& query, ///< 3D point to project
								     Tree& tree, ///< KD-tree
								     const unsigned int K,
                     const unsigned int degre_fitting,
										 const unsigned int degree_monge)
{
  // basic geometric types
  typedef typename Kernel::Vector_3 Vector;

	// types for K nearest neighbor search
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

	// types for jet fitting
	typedef typename CGAL::Monge_via_jet_fitting<Kernel> Monge_jet_fitting;
  typedef typename Monge_jet_fitting::Monge_form Monge_form;

	// gather set of (K+1) neighboring points
	std::vector<Point> points;

	// performs K + 1 queries (if unique the query point is
	// output first). search may be aborted when K is greater
	// than number of input points
  Neighbor_search search(tree,query,K+1);
	Search_iterator search_iterator = search.begin();
	unsigned int i;
	for(i=0;i<(K+1);i++)
	{
		if(search_iterator == search.end())
			break; // premature ending
		points.push_back(search_iterator->first);
		search_iterator++;
	}
	CGAL_precondition(points.size() >= 1);

	// performs jet fitting
	Monge_jet_fitting monge_fit;
	Monge_form monge_form = monge_fit(points.begin(), points.end(),
		                                degre_fitting, degree_monge);

	// output projection of query point onto the jet
	return monge_form.origin();
}

/// Smooth a point set using jet fitting on the K
/// nearest neighbors and reprojection onto the jet.
/// This variant requires the kernel.
///
/// Precondition: K >= 2.
template < typename InputIterator, ///< InputIterator value_type is Point_3.
           typename OutputIterator, ///< OutputIterator value_type is Point_3.
           typename Kernel ///< Geometric traits class.
>
OutputIterator ///< return past-the-end iterator of output
smooth_jet_fitting_3(InputIterator first,    ///< input points
                     InputIterator beyond,
						         OutputIterator output, ///< output points
						         const unsigned int K,   ///< number of neighbors
				      	  	 const Kernel& /*kernel*/,
										 const unsigned int degre_fitting = 2,
										 const unsigned int degree_monge = 2)
{
	// types for K-nearest neighbor search structure
	typedef typename Kernel::Point_3 Point;
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
	// but this is costly to check
  CGAL_precondition(first != beyond);

	// precondition: at least 2 nearest neighbors
  CGAL_precondition(K >= 2);

	// instanciate a KD-tree search
  Tree tree(first,beyond);

	// iterate over input points, compute and output smooth points
	InputIterator it;
	for(it = first; it != beyond; it++)
	{
		*output = smooth_jet_fitting_3<Kernel,Tree,Point>(*it,tree,K,degre_fitting,degree_monge);
		output++;
	}
	return output;
}

/// Smooth a point set using jet fitting on the K
/// nearest neighbors and reprojection onto the jet.
/// This variant requires the kernel.
///
/// Precondition: K >= 2.
template < typename InputIterator, ///< InputIterator value_type is Point_3.
           typename Kernel ///< Geometric traits class.
>
void
smooth_jet_fitting_3(InputIterator first,    ///< input points
                     InputIterator beyond,
						         const unsigned int K,   ///< number of neighbors
				      	  	 const Kernel& /*kernel*/,
										 const unsigned int degre_fitting = 2,
										 const unsigned int degree_monge = 2)
{
	// types for K-nearest neighbor search structure
	typedef typename Kernel::Point_3 Point;
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
	// but this is costly to check
  CGAL_precondition(first != beyond);

	// precondition: at least 2 nearest neighbors
  CGAL_precondition(K >= 2);

	// instanciate a KD-tree search
  Tree tree(first,beyond);

	// iterate over input points and mutate them
	InputIterator it;
	for(it = first; it != beyond; it++)
		*it = smooth_jet_fitting_3<Kernel,Tree,Point>(*it,tree,K,degre_fitting,degree_monge);
}


/// Smooth a point set using jet fitting on the K
/// nearest neighbors and reprojection onto the jet.
/// This variant deduces the kernel from iterator types.
///
/// Precondition: K >= 2.
template < typename InputIterator, ///< InputIterator value_type is Point_3
           typename OutputIterator ///< OutputIterator value_type is Point_3
>
OutputIterator ///< return past-the-end iterator of output
smooth_jet_fitting_3(InputIterator first,    ///< input points
                     InputIterator beyond,
										 OutputIterator output, ///< output points
										 const unsigned int K,   ///< number of neighbors
										 const unsigned int degre_fitting = 2,
										 const unsigned int degree_monge = 2)
{
	typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
	return smooth_jet_fitting_3(first,beyond,output,K,Kernel(),degre_fitting,degree_monge);
}

/// Smooth a point set using jet fitting on the K
/// nearest neighbors and reprojection onto the jet.
/// This function is mutating the input point set.
/// This variant deduces the kernel from iterator types.
///
/// Precondition: K >= 2.
template < typename InputIterator ///< InputIterator value_type is Point_3
>
void 
smooth_jet_fitting_3(InputIterator first,    ///< input points
                     InputIterator beyond,
										 const unsigned int K,   ///< number of neighbors
										 const unsigned int degre_fitting = 2,
										 const unsigned int degree_monge = 2)
{
	typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
	smooth_jet_fitting_3(first,beyond,K,Kernel(),degre_fitting,degree_monge);
}



CGAL_END_NAMESPACE

#endif // CGAL_SMOOTH_JET_FITTING_3_H


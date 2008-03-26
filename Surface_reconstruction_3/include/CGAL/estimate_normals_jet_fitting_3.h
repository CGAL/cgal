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
// $URL$
// $Id$
//
// Author(s) : Pierre Alliez and Laurent Saboret and Marc Pouget and Frederic Cazals

#ifndef CGAL_ESTIMATE_NORMALS_JET_FITTING_3_H
#define CGAL_ESTIMATE_NORMALS_JET_FITTING_3_H

#include <CGAL/basic.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Oriented_normal_3.h>

#include <iterator>
#include <list>

CGAL_BEGIN_NAMESPACE

/// Estimate normal direction using linear least
/// squares fitting of a plane on the K nearest neighbors.
///
/// Precondition: K >= 2.
///
/// @return Computed normal, model of OrientedNormal_3.
template < typename Kernel, ///< Geometric traits class.
           typename Tree,
           typename OrientedNormal_3
>
OrientedNormal_3
estimate_normal_jet_fitting_3(const typename Kernel::Point_3& query, ///< 3D point whose normal we want to compute
								              Tree& tree, ///< KD-tree
								              const unsigned int K,
                              const unsigned int degre_fitting = 2)
{
  // basic geometric types
  typedef typename Kernel::Point_3  Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef OrientedNormal_3 Oriented_normal;

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
	const unsigned int degree_monge = 1; // we seek for normal and not more.
	Monge_form monge_form = monge_fit(points.begin(), points.end(),
		                                degre_fitting, degree_monge);

	// output normal vector (already normalized in monge form)
	return OrientedNormal_3(monge_form.normal_direction(),
                          false /* not oriented */);
}


/// Estimate normal directions using jet fitting on the K nearest
/// neighbors.
/// This variant requires the kernel.
///
/// Precondition: K >= 2.
template < typename InputIterator, ///< InputIterator value_type is Point_3.
           typename OutputIterator, ///< OutputIterator value_type is a model of OrientedNormal_3.
           typename Kernel ///< Geometric traits class.
>
void
estimate_normals_jet_fitting_3(InputIterator first,    ///< input points
                               InputIterator beyond,
											         OutputIterator normals, ///< output normals
											         const unsigned int K,   ///< number of neighbors
									      	  	 const Kernel& /*kernel*/,
															 const unsigned int degre_fitting = 2)
{
  // Hard-code the Normal type as back_insert_iterator value_type is wrong (VC++ 2003)
	// typedef typename std::iterator_traits<OutputIterator>::value_type Normal;
  typedef CGAL::Oriented_normal_3<Kernel> Normal;

	// types for K-nearest neighbor search structure
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

	// iterate over input points, compute and output normal
	// vectors (already normalized)
	InputIterator it;
	for(it = first; it != beyond; it++)
	{
		*normals = estimate_normal_jet_fitting_3<Kernel,Tree,Normal>(*it,tree,K,degre_fitting);
		normals++;
	}
}

/// Estimate normal directions using jet fitting on the K nearest
/// neighbors.
/// This variant deduces the kernel from iterator types.
///
/// Precondition: K >= 2.
template < typename InputIterator, ///< InputIterator value_type is Point_3
           typename OutputIterator ///< OutputIterator value_type is a model of OrientedNormal_3
>
void
estimate_normals_jet_fitting_3(InputIterator first,    ///< input points
                               InputIterator beyond,
											         OutputIterator normals, ///< output normals
											         const unsigned int K,   ///< number of neighbors
															 const unsigned int degre_fitting = 2)
{
	typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
	estimate_normals_jet_fitting_3(first,beyond,normals,K,Kernel(),degre_fitting);
}


CGAL_END_NAMESPACE

#endif // CGAL_ESTIMATE_NORMALS_JET_FITTING_3_H


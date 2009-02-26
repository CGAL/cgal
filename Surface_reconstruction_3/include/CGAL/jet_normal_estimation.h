// Copyright (c) 2007-09  INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_JET_NORMAL_ESTIMATION_H
#define CGAL_JET_NORMAL_ESTIMATION_H

#include <CGAL/value_type_traits.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Orientable_normal_3.h>
#include <CGAL/surface_reconstruction_assertions.h>
#include <CGAL/Memory_sizer.h>

#include <iterator>
#include <list>

CGAL_BEGIN_NAMESPACE


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi { 


/// Estimate normal direction using jet fitting
/// on the k nearest neighbors.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param Kernel Geometric traits class.
/// @param Tree KD-tree.
/// @param OrientableNormal_3 Type of return value.
///
/// @return Computed normal, model of OrientableNormal_3.
template < typename Kernel,
           typename Tree,
           typename OrientableNormal_3
>
OrientableNormal_3
jet_normal_estimation(const typename Kernel::Point_3& query, ///< point to compute the normal at
                      Tree& tree, ///< KD-tree
                      unsigned int k, ///< number of neighbors
                      unsigned int degre_fitting)
{
  // basic geometric types
  typedef typename Kernel::Point_3  Point;
  typedef typename Kernel::Vector_3 Vector;
  typedef OrientableNormal_3 Oriented_normal;

  // types for K nearest neighbors search
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // types for jet fitting
  typedef typename CGAL::Monge_via_jet_fitting<Kernel> Monge_jet_fitting;
  typedef typename Monge_jet_fitting::Monge_form Monge_form;

  // Gather set of (k+1) neighboring points.
  // Perform k+1 queries (as in point set, the query point is
  // output first). Search may be aborted when k is greater
  // than number of input points.
  std::vector<Point> points; points.reserve(k+1);
  Neighbor_search search(tree,query,k+1);
  Search_iterator search_iterator = search.begin();
  unsigned int i;
  for(i=0;i<(k+1);i++)
  {
    if(search_iterator == search.end())
      break; // premature ending
    points.push_back(search_iterator->first);
    search_iterator++;
  }
  CGAL_surface_reconstruction_precondition(points.size() >= 1);

  // performs jet fitting
  Monge_jet_fitting monge_fit;
  const unsigned int degree_monge = 1; // we seek for normal and not more.
  Monge_form monge_form = monge_fit(points.begin(), points.end(),
                                    degre_fitting, degree_monge);

  // output normal vector (already normalized in monge form)
  return OrientableNormal_3(monge_form.normal_direction(),
                            false /* not oriented */);
}


} /* namespace CGALi */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Estimate normal directions using jet fitting on the k nearest
/// neighbors.
/// This variant requires the kernel.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param InputIterator value_type is Point_3.
/// @param OutputIterator value_type is Point_3.
/// @param Kernel Geometric traits class.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
jet_normal_estimation(InputIterator first, ///< input points
                      InputIterator beyond,
                      OutputIterator normals, ///< output normals
                      unsigned int k, ///< number of neighbors
                      const Kernel& /*kernel*/,
                      unsigned int degre_fitting)
{
  CGAL_TRACE("Call jet_normal_estimation()\n");

  // value_type_traits is a workaround as back_insert_iterator's value_type is void
  typedef typename value_type_traits<OutputIterator>::type Normal;

  // types for K nearest neighbors search structure
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_surface_reconstruction_precondition(first != beyond);

  // precondition: at least 2 nearest neighbors
  CGAL_surface_reconstruction_precondition(k >= 2);

  long memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Create KD-tree\n");

  // instanciate a KD-tree search
  Tree tree(first,beyond);

  /*long*/ memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("  Compute normals\n");

  // iterate over input points, compute and output normal
  // vectors (already normalized)
  InputIterator it;
  for(it = first; it != beyond; it++)
  {
    *normals = CGALi::jet_normal_estimation<Kernel,Tree,Normal>(*it,tree,k,degre_fitting);
    normals++;
  }

  /*long*/ memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("End of jet_normal_estimation()\n");

  return normals;
}

/// Estimate normal directions using jet fitting on the k nearest
/// neighbors.
/// This variant deduces the kernel from iterator types.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param InputIterator value_type is Point_3.
/// @param OutputIterator value_type is Point_3.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator
>
OutputIterator
jet_normal_estimation(InputIterator first, ///< input points
                      InputIterator beyond,
                      OutputIterator normals, ///< output normals
                      unsigned int k, ///< number of neighbors
                      unsigned int degre_fitting = 2)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return jet_normal_estimation(first,beyond,normals,k,Kernel(),degre_fitting);
}


CGAL_END_NAMESPACE

#endif // CGAL_JET_NORMAL_ESTIMATION_H


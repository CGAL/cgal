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
// Author(s) : Pierre Alliez and Marc Pouget

#ifndef CGAL_JET_SMOOTHING_3_H
#define CGAL_JET_SMOOTHING_3_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Monge_via_jet_fitting.h>

#include <iterator>
#include <list>

CGAL_BEGIN_NAMESPACE


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi {


/// Smooth one point position using jet fitting on the KNN
/// nearest neighbors and reprojection onto the jet.
///
/// @commentheading Precondition: KNN >= 2.
///
/// @commentheading Template Parameters:
/// @param Kernel Geometric traits class.
/// @param Tree KD-tree.
///
/// @return computed point
template <typename Kernel,
          typename Tree>
typename Kernel::Point_3
jet_smoothing_3(const typename Kernel::Point_3& query, ///< 3D point to project
                Tree& tree, ///< KD-tree
                const unsigned int KNN,
                const unsigned int degre_fitting,
                const unsigned int degree_monge)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point_3;

  // types for K nearest neighbors search
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // types for jet fitting
  typedef typename CGAL::Monge_via_jet_fitting<Kernel> Monge_jet_fitting;
  typedef typename Monge_jet_fitting::Monge_form Monge_form;

  // Gather set of (KNN+1) neighboring points.
  // Performs KNN + 1 queries (if unique the query point is
  // output first). Search may be aborted when KNN is greater
  // than number of input points.
  std::vector<Point_3> points; points.reserve(KNN+1);
  Neighbor_search search(tree,query,KNN+1);
  Search_iterator search_iterator = search.begin();
  unsigned int i;
  for(i=0;i<(KNN+1);i++)
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


} /* namespace CGALi */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Smooth a point set using jet fitting on the KNN
/// nearest neighbors and reprojection onto the jet.
/// This variant requires the kernel.
///
/// @commentheading Precondition: KNN >= 2.
///
/// @commentheading Template Parameters:
/// @param InputIterator value_type convertible to Point_3.
/// @param OutputIterator value_type convertible to Point_3.
/// @param Kernel Geometric traits class.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
jet_smoothing_3(InputIterator first,    ///< input points
                InputIterator beyond,
                OutputIterator output,  ///< output points
                const unsigned int KNN, ///< number of neighbors
                const Kernel& /*kernel*/,
                const unsigned int degre_fitting = 2,
                const unsigned int degree_monge = 2)
{
  // Point_3 types
  typedef typename std::iterator_traits<InputIterator>::value_type Input_point_3;
  typedef typename Kernel::Point_3 Point_3;

  // types for K nearest neighbors search structure
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_precondition(first != beyond);

  // precondition: at least 2 nearest neighbors
  CGAL_precondition(KNN >= 2);

  // instanciate a KD-tree search
  Tree tree(first,beyond);

  // Iterate over input points, compute and output smooth points.
  // Note: the cast to (Point_3&) ensures compatibility with classes derived from Point_3.
  for(InputIterator it = first; it != beyond; it++)
  {
    Input_point_3 point = *it;
    (Point_3&)(point) = CGALi::jet_smoothing_3<Kernel>(*it,tree,KNN,degre_fitting,degree_monge);
    *output++ = point;
  }

  return output;
}

/// Smooth a point set using jet fitting on the KNN
/// nearest neighbors and reprojection onto the jet.
/// This function is mutating the input point set.
/// This variant requires the kernel.
///
/// Warning:
/// This method moves the points, thus
/// should not be called on containers sorted wrt points position.
///
/// @commentheading Precondition: KNN >= 2.
///
/// @commentheading Template Parameters:
/// @param ForwardIterator value_type convertible to Point_3.
/// @param Kernel Geometric traits class.
template <typename ForwardIterator,
          typename Kernel>
void
jet_smoothing_3(ForwardIterator first,     ///< input/output points
                ForwardIterator beyond,
                unsigned int KNN,          ///< number of neighbors
                const Kernel& /*kernel*/,
                const unsigned int degre_fitting = 2,
                const unsigned int degree_monge = 2)
{
  // Point_3 types
  typedef typename std::iterator_traits<ForwardIterator>::value_type Input_point_3;
  typedef typename Kernel::Point_3 Point_3;

  // types for K nearest neighbors search structure
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_precondition(first != beyond);

  // precondition: at least 2 nearest neighbors
  CGAL_precondition(KNN >= 2);

  // instanciate a KD-tree search
  Tree tree(first,beyond);

  // Iterate over input points and mutate them.
  // Note: the cast to (Point_3&) ensures compatibility with classes derived from Point_3.
  ForwardIterator it;
  for(it = first; it != beyond; it++)
    (Point_3&)(*it) = CGALi::jet_smoothing_3<Kernel>(*it,tree,KNN,degre_fitting,degree_monge);
}


/// Smooth a point set using jet fitting on the KNN
/// nearest neighbors and reprojection onto the jet.
/// This variant deduces the kernel from iterator types.
///
/// @commentheading Precondition: KNN >= 2.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator
>
OutputIterator
jet_smoothing_3(InputIterator first,    ///< input points
                InputIterator beyond,
                OutputIterator output, ///< output points
                unsigned int KNN,      ///< number of neighbors
                const unsigned int degre_fitting = 2,
                const unsigned int degree_monge = 2)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Input_point_3;
  typedef typename Kernel_traits<Input_point_3>::Kernel Kernel;
  return jet_smoothing_3(first,beyond,output,KNN,Kernel(),degre_fitting,degree_monge);
}

/// Smooth a point set using jet fitting on the KNN
/// nearest neighbors and reprojection onto the jet.
/// This function is mutating the input point set.
/// This variant deduces the kernel from iterator types.
///
/// Warning:
/// This method moves the points, thus
/// should not be called on containers sorted wrt points position.
///
/// @commentheading Precondition: KNN >= 2.
///
/// @commentheading Template Parameters:
/// @param ForwardIterator value_type convertible to Point_3.
template <typename ForwardIterator>
void
jet_smoothing_3(ForwardIterator first,     ///< input/output points
                ForwardIterator beyond,
                unsigned int KNN,          ///< number of neighbors
                const unsigned int degre_fitting = 2,
                const unsigned int degree_monge = 2)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type Input_point_3;
  typedef typename Kernel_traits<Input_point_3>::Kernel Kernel;
  jet_smoothing_3(first,beyond,KNN,Kernel(),degre_fitting,degree_monge);
}


CGAL_END_NAMESPACE

#endif // CGAL_JET_SMOOTHING_3_H


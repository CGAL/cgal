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
// Author(s) : Nader Salman and Laurent Saboret

#ifndef CGAL_IMPROVED_LAPLACIAN_SMOOTH_POINT_SET_H
#define CGAL_IMPROVED_LAPLACIAN_SMOOTH_POINT_SET_H

#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <iterator>
#include <list>

CGAL_BEGIN_NAMESPACE


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace improved_laplacian_smoothing_i {


// Item in the Kd-tree: position (Point_3) + index
template <typename Kernel>
class KdTreeElement : public Kernel::Point_3
{
public:
  unsigned int index;

  // basic geometric types
  typedef typename CGAL::Origin Origin;
  typedef typename Kernel::Point_3 Point_3;

  KdTreeElement(const Origin& o = ORIGIN, unsigned int id=0)
    : Point_3(o), index(id)
  {}
  KdTreeElement(const Point_3& p, unsigned int id=0)
    : Point_3(p), index(id)
  {}
  KdTreeElement(const KdTreeElement& other)
    : Point_3(other), index(other.index)
  {}
};

// Helper class for the Kd-tree
template <typename Kernel>
class KdTreeGT : public Kernel
{
public:
  typedef KdTreeElement<Kernel> Point_3;
};

template <typename Kernel>
class KdTreeTraits : public CGAL::Search_traits_3<KdTreeGT<Kernel> >
{
public:
    typedef typename Kernel::Point_3 PointType;
};

/* Usage:
  typedef improved_laplacian_smoothing_i::KdTreeElement<Kernel> KdTreeElement;
  typedef improved_laplacian_smoothing_i::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;
*/


/// Smooth one point position using jet fitting on the k
/// nearest neighbors and reprojection onto the jet.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param Kernel Geometric traits class.
/// @param Tree KD-tree.
///
/// @return computed point
template <typename Kernel,
          typename Tree>
typename Kernel::Point_3
laplacian_smoothing_3(
                const typename Kernel::Point_3& pi, ///< 3D point to smooth
                Tree& tree, ///< KD-tree
                const unsigned int k)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;

  // types for K nearest neighbors search
  //typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef improved_laplacian_smoothing_i::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // Compute Laplacian (centroid) of k neighboring points.
  // Note: we perform k+1 queries and skip the query point which is output first.
  // TODO: we should use the functions in PCA component instead.
  Vector_3 v = CGAL::NULL_VECTOR;
  Neighbor_search search(tree,pi,k+1);
  Search_iterator search_iterator;
  for(search_iterator = search.begin(), search_iterator++; // skip pi point
      search_iterator != search.end();
      search_iterator++ )
  {
      const Point_3& p = search_iterator->first;
      v = v + (p - CGAL::ORIGIN);
  }

  Point_3 centroid = CGAL::ORIGIN + v / k;

  return centroid;
}


/// Smooth one point position using jet fitting on the k
/// nearest neighbors and reprojection onto the jet.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param Kernel Geometric traits class.
/// @param Tree KD-tree.
///
/// @return computed point
template <typename Kernel,
          typename Tree>
typename Kernel::Point_3
improved_laplacian_smooth_point_set(
                const typename Kernel::Point_3& pi, ///< 3D point to smooth
                const typename Kernel::Vector_3& bi, ///< bi movement
                Tree& tree, ///< KD-tree
                const std::vector<typename Kernel::Vector_3>& b,
                const unsigned int k,
                typename Kernel::FT beta)
{
  // basic geometric types
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;

  // types for K nearest neighbors search
  //typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef improved_laplacian_smoothing_i::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // Gather set of k neighboring points and compute the sum of their b[] values.
  // Note: we perform k+1 queries and skip the query point which is output first.
  Vector_3 bj_sum;
  Neighbor_search search(tree,pi,k+1);
  Search_iterator search_iterator;
  for(search_iterator = search.begin(), search_iterator++; // skip pi point
      search_iterator != search.end();
      search_iterator++ )
  {
    bj_sum = bj_sum + b[search_iterator->first.index];
  }

  return pi - (beta * bi + ((1-beta)/k)*bj_sum);
}


} /* namespace improved_laplacian_smoothing_i */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Improved Laplacian smoothing (Vollmer et al)
/// on the k nearest neighbors.
/// This variant requires the kernel.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param InputIterator value_type must be convertible to Point_3<Kernel>.
/// @param OutputIterator value_type must be convertible from InputIterator's value_type.
/// @param Kernel Geometric traits class.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
improved_laplacian_smooth_point_set(
                InputIterator first,    ///< iterator over the first input point.
                InputIterator beyond,   ///< past-the-end iterator over input points.
                OutputIterator output,  ///< iterator over the first output point.
                const unsigned int k,   ///< number of neighbors.
                const unsigned int iter_number, ///< number of iterations.
                const Kernel& kernel,   ///< geometric traits.
                typename Kernel::FT alpha,
                typename Kernel::FT beta)
{
  // Point_3 types
  typedef typename std::iterator_traits<InputIterator>::value_type Input_point_3;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;

  // types for K nearest neighbors search structure
  //typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef improved_laplacian_smoothing_i::KdTreeElement<Kernel> KdTreeElement;
  typedef improved_laplacian_smoothing_i::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_precondition(first != beyond);

  // precondition: at least 2 nearest neighbors
  CGAL_precondition(k >= 2);

  unsigned int i; // point index
  InputIterator it; // point iterator

  // Number of input points
  int nb_points = std::distance(first, beyond);

  // Create kd-tree
  //Tree tree(first,beyond);
  std::vector<KdTreeElement> treeElements;
  for (it = first, i=0 ; it != beyond ; ++it,++i)
  {
    treeElements.push_back(KdTreeElement(*it,i));
  }
  Tree tree(treeElements.begin(), treeElements.end());

  std::vector<Point_3>  p(nb_points); // positions at step iter_n
  std::vector<Vector_3> b(nb_points); // ...

  for(it = first, i=0; it != beyond; it++, ++i)
      p[i] = *it;

  for(int iter_n = 0; iter_n < iter_number ; ++iter_n)
  {
      // Iterate over input points, compute (original) Laplacian smooth and b[].
      for(it = first, i=0; it != beyond; it++, ++i)
      {
          Point_3 np = improved_laplacian_smoothing_i::laplacian_smoothing_3<Kernel>(*it,tree,k);
          b[i]       = alpha*(np - *it) + (1-alpha)*(np - p[i]);
          p[i]       = np;
      }

      // Iterate over input points, compute and output smooth points.
      // Note: the cast to (Point_3&) ensures compatibility with classes derived from Point_3.
      for(it = first, i=0; it != beyond; it++, ++i)
      {
          p[i] = improved_laplacian_smoothing_i::improved_laplacian_smooth_point_set<Kernel>(p[i],b[i],tree,b,k,beta);
      }
  }

  // Iterate over input points and output smooth points.
  // Note: the cast to (Point_3&) ensures compatibility with classes derived from Point_3.
  for(it = first, i=0; it != beyond; it++, ++i)
  {
    Input_point_3 point = *it;
    (Point_3&)(point) = p[i];
    *output++ = point;
  }

  return output;
}

/// Improved Laplacian smoothing (Vollmer et al)
/// on the k nearest neighbors.
/// This variant is mutating the input point set and requires the kernel.
///
/// Warning:
/// This method moves the points, thus
/// should not be called on containers sorted wrt points position.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param ForwardIterator value_type must be convertible to Point_3<Kernel>.
/// @param Kernel Geometric traits class.
template <typename ForwardIterator,
          typename Kernel>
void
improved_laplacian_smooth_point_set(
                ForwardIterator first,     ///< iterator over the first input/output point.
                ForwardIterator beyond,    ///< past-the-end iterator.
                const unsigned int k,      ///< number of neighbors.
                const unsigned int iter_number, ///< number of iterations.
                const Kernel& kernel,      ///< geometric traits.
                typename Kernel::FT alpha,
                typename Kernel::FT beta)
{
  // Point_3 types
  typedef typename std::iterator_traits<ForwardIterator>::value_type Input_point_3;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;

  // types for K nearest neighbors search structure
  //typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef improved_laplacian_smoothing_i::KdTreeElement<Kernel> KdTreeElement;
  typedef improved_laplacian_smoothing_i::KdTreeTraits<Kernel> Tree_traits;
  typedef CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::Tree Tree;
  typedef typename Neighbor_search::iterator Search_iterator;

  // precondition: at least one element in the container.
  // to fix: should have at least three distinct points
  // but this is costly to check
  CGAL_precondition(first != beyond);

  // precondition: at least 2 nearest neighbors
  CGAL_precondition(k >= 2);

  unsigned int i; // point index
  ForwardIterator it; // point iterator

  // Number of input points
  int nb_points = std::distance(first, beyond);

  // Create kd-tree
  //Tree tree(first,beyond);
  std::vector<KdTreeElement> treeElements;
  for (it = first, i=0 ; it != beyond ; ++it,++i)
  {
    treeElements.push_back(KdTreeElement(*it,i));
  }
  Tree tree(treeElements.begin(), treeElements.end());

  std::vector<Point_3>  p(nb_points); // positions at step iter_n
  std::vector<Vector_3> b(nb_points); // ...

  for(it = first, i=0; it != beyond; it++, ++i)
      p[i] = *it;

  for(int iter_n = 0; iter_n < iter_number ; ++iter_n)
  {
      // Iterate over input points, compute (original) Laplacian smooth and b[].
      for(it = first, i=0; it != beyond; it++, ++i)
      {
          Point_3 np = improved_laplacian_smoothing_i::laplacian_smoothing_3<Kernel>(*it,tree,k);
          b[i]       = alpha*(np - *it) + (1-alpha)*(np - p[i]);
          p[i]       = np;
      }

      // Iterate over input points, compute and output smooth points.
      // Note: the cast to (Point_3&) ensures compatibility with classes derived from Point_3.
      for(it = first, i=0; it != beyond; it++, ++i)
      {
          p[i] = improved_laplacian_smoothing_i::improved_laplacian_smooth_point_set<Kernel>(p[i],b[i],tree,b,k,beta);
      }
  }

  // Iterate over input points and mutate them.
  // Note: the cast to (Point_3&) ensures compatibility with classes derived from Point_3.
  for(it = first, i=0; it != beyond; it++, ++i)
    (Point_3&)(*it) = p[i];
}


/// Improved Laplacian smoothing (Vollmer et al)
/// on the k nearest neighbors.
/// This variant deduces the kernel from iterator types.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param InputIterator value_type must be convertible to Point_3<Kernel>.
/// @param OutputIterator value_type must be convertible from InputIterator's value_type.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator
>
OutputIterator
improved_laplacian_smooth_point_set(
                InputIterator first, ///< iterator over the first input point.
                InputIterator beyond, ///< past-the-end iterator over input points.
                OutputIterator output, ///< iterator over the first output point.
                unsigned int k, ///< number of neighbors.
                const unsigned int iter_number, ///< number of iterations.
                double alpha,
                double beta)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Input_point_3;
  typedef typename Kernel_traits<Input_point_3>::Kernel Kernel;
  return improved_laplacian_smooth_point_set(first,beyond,output,k,iter_number,Kernel(),alpha, beta);
}

/// Improved Laplacian smoothing (Vollmer et al)
/// on the k nearest neighbors.
/// This variant is mutating the input point set and deduces the kernel from iterator types.
///
/// Warning:
/// As this method relocates the points, it
/// should not be called on containers sorted w.r.t. point locations.
///
/// @commentheading Precondition: k >= 2.
///
/// @commentheading Template Parameters:
/// @param ForwardIterator value_type must be convertible to Point_3<Kernel>.
template <typename ForwardIterator>
void
improved_laplacian_smooth_point_set(
                ForwardIterator first, ///< iterator over the first input/output point.
                ForwardIterator beyond, ///< past-the-end iterator.
                unsigned int k, ///< number of neighbors.
                const unsigned int iter_number, ///< number of iterations.
                double alpha,
                double beta)
{
  typedef typename std::iterator_traits<ForwardIterator>::value_type Input_point_3;
  typedef typename Kernel_traits<Input_point_3>::Kernel Kernel;
  improved_laplacian_smooth_point_set(first,beyond,k,iter_number,Kernel(),alpha, beta);
}


CGAL_END_NAMESPACE

#endif // CGAL_IMPROVED_LAPLACIAN_SMOOTH_POINT_SET_H


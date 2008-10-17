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
// Author(s) : Pierre Alliez and Laurent Saboret

#ifndef CGAL_ESTIMATE_NORMALS_PCA_3_H
#define CGAL_ESTIMATE_NORMALS_PCA_3_H

#include <CGAL/value_type_traits.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Orientable_normal_3.h>
#include <CGAL/surface_reconstruction_assertions.h>
#include <CGAL/Memory_sizer.h>

#include <iterator>
#include <list>

CGAL_BEGIN_NAMESPACE


/// Estimate normal direction using linear least
/// squares fitting of a plane on the K nearest neighbors.
///
/// Precondition: KNN >= 2.
///
/// @heading Parameters:
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
estimate_normal_pca_3(const typename Kernel::Point_3& query, ///< 3D point whose normal we want to compute
                      Tree& tree, ///< KD-tree
                      const unsigned int KNN)
{
  // basic geometric types
  typedef typename Kernel::Point_3  Point;
  typedef typename Kernel::Plane_3  Plane;
  typedef typename Kernel::Vector_3 Vector;
  typedef OrientableNormal_3 Oriented_normal;

  // types for K nearest neighbors search
  typedef typename CGAL::Search_traits_3<Kernel> Tree_traits;
  typedef typename CGAL::Orthogonal_k_neighbor_search<Tree_traits> Neighbor_search;
  typedef typename Neighbor_search::iterator Search_iterator;

  // Gather set of (KNN+1) neighboring points.
  // Perform KNN+1 queries (as in point set, the query point is
  // output first). Search may be aborted when KNN is greater
  // than number of input points.
  std::vector<Point> points; points.reserve(KNN+1);
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
  CGAL_surface_reconstruction_precondition(points.size() >= 1);

  // performs plane fitting by point-based PCA
  Plane plane;
  linear_least_squares_fitting_3(points.begin(),points.end(),plane,Dimension_tag<0>());

  // output normal vector (already normalized by PCA)
  return OrientableNormal_3(plane.orthogonal_vector(),
                            false /* not oriented */);
}


/// Estimate normals direction using linear least
/// squares fitting of a plane on the K nearest neighbors.
/// This variant requires the kernel.
///
/// Precondition: KNN >= 2.
///
/// @heading Parameters:
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
estimate_normals_pca_3(InputIterator first,    ///< input points
                       InputIterator beyond,
                       OutputIterator normals, ///< output normals
                       const unsigned int KNN,   ///< number of neighbors
                       const Kernel& /*kernel*/)
{
  CGAL_TRACE("Call estimate_normals_pca_3()\n");

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
  CGAL_surface_reconstruction_precondition(KNN >= 2);

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
    *normals = estimate_normal_pca_3<Kernel,Tree,Normal>(*it,tree,KNN);
    normals++;
  }

  /*long*/ memory = CGAL::Memory_sizer().virtual_size(); CGAL_TRACE("  %ld Mb allocated\n", memory>>20);
  CGAL_TRACE("End of estimate_normals_pca_3()\n");

  return normals;
}

/// Estimate normals direction using linear least
/// squares fitting of a plane on the K nearest neighbors.
/// This variant deduces the kernel from iterator types.
///
/// Precondition: KNN >= 2.
///
/// @heading Parameters:
/// @param InputIterator value_type is Point_3.
/// @param OutputIterator value_type is Point_3.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator
>
OutputIterator
estimate_normals_pca_3(InputIterator first,    ///< input points
                       InputIterator beyond,
                       OutputIterator normals, ///< output normals
                       const unsigned int KNN)   ///< number of neighbors
{
  typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
  typedef typename Kernel_traits<Value_type>::Kernel Kernel;
  return estimate_normals_pca_3(first,beyond,normals,KNN,Kernel());
}


CGAL_END_NAMESPACE

#endif // CGAL_ESTIMATE_NORMALS_PCA_3_H


// Copyright (c) 2017 GeometryFactory Sarl (France).
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s):      Simon Giraudot


#ifndef CGAL_SCALE_SPACE_RECONSTRUCTION_3_JET_SMOOTHER_H
#define CGAL_SCALE_SPACE_RECONSTRUCTION_3_JET_SMOOTHER_H

#include <CGAL/license/Scale_space_reconstruction_3.h>

#include <CGAL/jet_smooth_point_set.h>

#ifdef CGAL_LINKED_WITH_TBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#endif // CGAL_LINKED_WITH_TBB

namespace CGAL
{

namespace Scale_space_reconstruction_3
{

/** \ingroup PkgScaleSpaceReconstruction3Classes
 *
 *  %Smoother for scale space reconstruction based on
 *  `CGAL::jet_smooth_point_set()`.
 *
 *  \cgalModels CGAL::Scale_space_reconstruction_3::Smoother
 *
 *  \tparam Geom_traits geometric traits class. It must be a
 *  model of `DelaunayTriangulationTraits_3`. It must have a
 *  `RealEmbeddable` field number type. Generally,
 *  `Exact_predicates_inexact_constructions_kernel` is preferred.
 *  \tparam ConcurrencyTag indicates whether to use concurrent
 *  processing. It can be omitted: if \ref thirdpartyTBB is available
 *  and `CGAL_LINKED_WITH_TBB` is defined then `Parallel_tag` is
 *  used. Otherwise, `Sequential_tag` is used.
 */
template <typename Geom_traits,
#ifdef DOXYGEN_RUNNING
          typename ConcurrencyTag>
#else
          typename ConcurrencyTag = CGAL::Parallel_if_available_tag>
#endif
class Jet_smoother
{
public:
  typedef typename Geom_traits::FT FT; ///< defines the point type.
  typedef typename Geom_traits::Point_3 Point; ///< defines the point typ.e
private:

  unsigned int m_k;
  unsigned int m_degree_fitting;
  unsigned int m_degree_monge;

public:

  /**
   * Constructs a jet smoother.
   *
   * \param k number of neighbors used.
   * \param degree_fitting fitting degree
   * \param degree_monge monge degree
   */
  Jet_smoother (unsigned int k = 12,
                unsigned int degree_fitting = 2,
                unsigned int degree_monge = 2)
    : m_k (k), m_degree_fitting (degree_fitting), m_degree_monge (degree_monge)
  { }

  template <typename InputIterator>
  void operator() (InputIterator begin, InputIterator end)
  {
    CGAL::Iterator_range<InputIterator> points (begin, end);
    CGAL::jet_smooth_point_set<ConcurrencyTag>
    (points, m_k, CGAL::parameters::degree_fitting(m_degree_fitting).degree_monge(m_degree_monge));
  }

};


} // namespace Scale_space_reconstruction_3

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_RECONSTRUCTION_3_JET_SMOOTHER_H

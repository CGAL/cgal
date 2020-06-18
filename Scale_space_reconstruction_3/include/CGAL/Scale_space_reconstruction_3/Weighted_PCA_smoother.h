// Copyright (C) 2013 INRIA - Sophia Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s):      Thijs van Lankveld, Simon Giraudot


#ifndef CGAL_SCALE_SPACE_RECONSTRUCTION_3_WEIGHTED_PCA_SMOOTHER_H
#define CGAL_SCALE_SPACE_RECONSTRUCTION_3_WEIGHTED_PCA_SMOOTHER_H

#include <CGAL/license/Scale_space_reconstruction_3.h>

#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_incremental_neighbor_search.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Random.h>
#include <CGAL/Default_diagonalize_traits.h>

#include <boost/function_output_iterator.hpp>

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
 *  %Smoother for scale space reconstruction based on a principal
 *  component analysis weighted by the local density of points.
 *
 *  \cgalModels CGAL::Scale_space_reconstruction_3::Smoother
 *
 *  \tparam Geom_traits geometric traits class. It must be a
 *  model of `DelaunayTriangulationTraits_3`. It must have a
 *  `RealEmbeddable` field number type. Generally,
 *  `Exact_predicates_inexact_constructions_kernel` is preferred.

 *  \tparam DiagonalizeTraits model of `DiagonalizeTraits` that
 *  determines how to diagonalize a weighted covariance matrix to
 *  approximate a weighted point set. It can be omitted if Eigen 3 (or
 *  greater) is available and `CGAL_EIGEN3_ENABLED` is defined: in
 *  that case, an overload using `Eigen_diagonalize_traits` is
 *  provided.
 *  \tparam ConcurrencyTag indicates whether to use concurrent
 *  processing. It can be omitted: if \ref thirdpartyTBB is available
 *  and `CGAL_LINKED_WITH_TBB` is defined then `Parallel_tag` is
 *  used. Otherwise, `Sequential_tag` is used.
 */
template <typename Geom_traits,
#ifdef DOXYGEN_RUNNING
          typename DiagonalizeTraits,
          typename ConcurrencyTag>
#else // DOXYGEN_RUNNING
          typename DiagonalizeTraits = CGAL::Default_diagonalize_traits<typename Geom_traits::FT, 3>,
          typename ConcurrencyTag = CGAL::Parallel_if_available_tag>
#endif // DOXYGEN_RUNNING
class Weighted_PCA_smoother
{
public:
  typedef typename Geom_traits::FT FT; ///< defines the point type.
  typedef typename Geom_traits::Point_3 Point; ///< defines the point typ.e
  typedef typename Geom_traits::Vector_3 Vector; ///< defines the vector type.
private:


  typedef boost::tuple<Point, std::size_t> Point_and_size_t;
  typedef std::vector<unsigned int> CountVec;
  typedef std::vector<Point> Pointset;

  typedef Search_traits_3<Geom_traits> Traits_base;
  typedef CGAL::Search_traits_adapter<Point_and_size_t,
                                      CGAL::Nth_of_tuple_property_map<0, Point_and_size_t>,
                                      Traits_base>                                              Search_traits;
  typedef Orthogonal_k_neighbor_search<Search_traits> Static_search;
  typedef Orthogonal_incremental_neighbor_search<Search_traits> Dynamic_search;
  typedef typename Dynamic_search::Tree               Search_tree;
  typedef Fuzzy_sphere< Search_traits >               Sphere;

  Random _generator;
  unsigned int _mean_neighbors;
  unsigned int _samples;
  FT _squared_radius;
  Search_tree _tree;  // To quickly search for nearest neighbors.
  Pointset _points;


public:

  /**
   * Constructs a weighted PCA smoother that will automatically
   * estimate the neighborhood radius.
   *
   * \param neighbors is the number of neighbors a point's neighborhood should
   * contain on average.
   * \param samples is the number of points sampled to estimate the
   * neighborhood radius.
   */
  Weighted_PCA_smoother (unsigned int neighbors = 12, unsigned int samples = 300)
    : _generator(0), _mean_neighbors (neighbors), _samples (samples),_squared_radius (-1.) { }

  /**
   * Constructs a weighted PCA smoother.
   *
   * \param squared_radius neighborhood squared radius used for principal component analysis.
   */
  Weighted_PCA_smoother (FT squared_radius)
    : _generator(0), _mean_neighbors (0), _samples (0),_squared_radius (squared_radius) { }

  /// \cond SKIP_IN_MANUAL
  template <typename InputIterator>
  void operator() (InputIterator begin, InputIterator end)
  {
    _tree.clear();
    _points.clear();

    std::size_t i = 0;
    std::size_t size = std::size_t(end - begin);
    _tree.reserve(size);
    _points.reserve(size);

    for (InputIterator it = begin; it != end; ++ it)
    {
      _points.push_back (*it);
      _tree.insert( boost::make_tuple(*it,i ++) );
    }

    _tree.build();

    if (_squared_radius == -1)
      estimate_neighborhood_squared_radius();

    // Collect the number of neighbors of each point.
    // This can be done concurrently.
    CountVec neighbors (_tree.size(), 0);
    try_parallel (ComputeNN (_points, _tree, _squared_radius, neighbors), 0, _tree.size());

    // Compute the transformed point locations.
    // This can be done concurrently.
    try_parallel (AdvanceSS (_tree, neighbors, _points), 0, _tree.size());

    i = 0;
    for (InputIterator it = begin; it != end; ++ it)
      *it = _points[i ++];
  }
  /// \endcond

  /**
   * Returns the computed (or user-specified) squared radius.
   */
  FT squared_radius()
  {
    return _squared_radius;
  }

private:

  void estimate_neighborhood_squared_radius ()
  {
    typename Geom_traits::Compute_squared_distance_3 squared_distance
      = Geom_traits().compute_squared_distance_3_object();

    unsigned int handled = 0;
    unsigned int checked = 0;
    FT radius = 0.;

    for (typename Search_tree::const_iterator it = _tree.begin(); it != _tree.end(); ++it )
    {
      unsigned int left = (unsigned int)(_tree.size() - handled);
      if (_samples >= left || _generator.get_double() < double(_samples - checked) / left)
      {
        // The neighborhood should contain the point itself as well.
        Static_search search (_tree, boost::get<0>(*it), _mean_neighbors + 1);
        radius += std::sqrt (to_double (squared_distance( boost::get<0>(*it), boost::get<0>((search.end()-1)->first))));
        ++checked;
      }
      ++handled;
    }
    radius /= double(checked);

    _squared_radius = radius * radius;
  }


  template <typename F>
  void try_parallel (const F& func, std::size_t begin, std::size_t end)
  {
#ifndef CGAL_LINKED_WITH_TBB
    CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                               "Parallel_tag is enabled but TBB is unavailable.");
#else
    if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
      tbb::parallel_for(tbb::blocked_range<std::size_t>(begin, end), func);
    else
#endif
      for (std::size_t i = begin; i < end; ++i)
        func(i);
  }

  struct Inc
  {
    unsigned int * i;

    Inc(unsigned int& i)
      : i(&i)
    {}

    template <typename T>
    void operator()(const T&) const
    {
      ++(*i);
    }

  };

  // Compute the number of neighbors of a point that lie within a fixed radius.
  class ComputeNN
  {
  private:
    typename Geom_traits::Compare_squared_distance_3 compare;

    const Pointset&     _pts;
    const Search_tree&  _tree;
    const FT           _sq_rd;
    CountVec&           _nn;

  public:
    ComputeNN(const Pointset& points, const Search_tree&  tree,
              const FT& sq_radius, CountVec& nn)
      : _pts(points), _tree(tree), _sq_rd(sqrt(sq_radius)), _nn(nn) {}

#ifdef CGAL_LINKED_WITH_TBB
    void operator()( const tbb::blocked_range< std::size_t >& range ) const {
      for( std::size_t i = range.begin(); i != range.end(); ++i )
        (*this)( i );
    }
#endif // CGAL_LINKED_WITH_TBB
    void operator()( const std::size_t& i ) const {

      Sphere sp(_pts[i], _sq_rd);

      Inc inc(_nn[i]);
      _tree.search(boost::make_function_output_iterator(inc),sp);
    }
  }; // class ComputeNN


// Advance a point to a coarser scale.
  class AdvanceSS
  {
  private:
    const Search_tree&  _tree;
    const CountVec&     _nn;
    Pointset&           _pts;

  public:
    AdvanceSS(const Search_tree& tree, const CountVec& nn, Pointset& points)
      : _tree(tree), _nn(nn),_pts(points) {}

#ifdef CGAL_LINKED_WITH_TBB
    void operator()( const tbb::blocked_range< std::size_t >& range ) const {
      for( std::size_t i = range.begin(); i != range.end(); ++i )
        (*this)( i );
    }
#endif // CGAL_LINKED_WITH_TBB
    void operator()( const std::size_t& i ) const {
      // If the neighborhood is too small, the vertex is not moved.
      if( _nn[i] < 4 )
        return;

      Static_search search(_tree, _pts[i], _nn[i]);

      Point barycenter (0., 0., 0.);
      FT weight_sum = 0.;
      unsigned int column = 0;
      // Compute total weight
      for( typename Static_search::iterator nit = search.begin();
           nit != search.end() && column < _nn[i];
           ++nit, ++column )
        weight_sum += (1.0 / _nn[boost::get<1>(nit->first)]);

      column = 0;
      // Compute barycenter
      for( typename Static_search::iterator nit = search.begin();
           nit != search.end() && column < _nn[i];
           ++nit, ++column )
      {
        Vector v (CGAL::ORIGIN, boost::get<0>(nit->first));
        barycenter = barycenter + ((1.0 / _nn[boost::get<1>(nit->first)]) / weight_sum) * v;
      }

      std::array<FT, 6> covariance = {{ 0., 0., 0., 0., 0., 0. }};
      column = 0;
      // Compute covariance matrix of Weighted PCA
      for( typename Static_search::iterator nit = search.begin();
           nit != search.end() && column < _nn[i];
           ++nit, ++column )
      {
        Vector v (barycenter, boost::get<0>(nit->first));
        FT w = (1.0 / _nn[boost::get<1>(nit->first)]);
        v = w*v;
        covariance[0] += w * v.x () * v.x ();
        covariance[1] += w * v.x () * v.y ();
        covariance[2] += w * v.x () * v.z ();
        covariance[3] += w * v.y () * v.y ();
        covariance[4] += w * v.y () * v.z ();
        covariance[5] += w * v.z () * v.z ();
      }

      // Compute the weighted least-squares planar approximation of the point set.
      std::array<FT, 9> eigenvectors = {{ 0., 0., 0.,
                                                  0., 0., 0.,
                                                  0., 0., 0. }};
      std::array<FT, 3> eigenvalues = {{ 0., 0., 0. }};
      DiagonalizeTraits::diagonalize_selfadjoint_covariance_matrix
        (covariance, eigenvalues, eigenvectors);

      // The vertex is moved by projecting it onto the plane
      // through the barycenter and orthogonal to the Eigen vector with smallest Eigen value.
      Vector norm (eigenvectors[0], eigenvectors[1], eigenvectors[2]);
      Vector b2p (barycenter, _pts[i]);

      _pts[i] = barycenter + b2p - ((norm * b2p) * norm);
    }
  }; // class AdvanceSS

};


} // namespace Scale_space_reconstruction_3

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_RECONSTRUCTION_3_WEIGHTED_PCA_SMOOTHER_H

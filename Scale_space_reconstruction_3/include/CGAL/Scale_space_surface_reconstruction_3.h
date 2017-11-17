// Copyright (C) 2013 INRIA - Sophia Antipolis (France).
// Copyright (c) 2017 GeometryFactory Sarl (France).
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s):      Thijs van Lankveld, Simon Giraudot


#ifndef CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_H
#define CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_H

#include <CGAL/license/Scale_space_reconstruction_3.h>

#include <CGAL/Scale_space_reconstruction_3/Weighted_PCA_smoother.h>
#include <CGAL/Scale_space_reconstruction_3/Alpha_shape_mesher.h>


namespace CGAL
{


/// computes a triangulated surface mesh interpolating a point set.
/** \ingroup PkgScaleSpaceReconstruction3Classes
 *
 *  Scale space reconstruction consists in combining a smoothing
 *  algorithm (see `Scale_space_reconstruction_3::Smoother`) that is
 *  used to generate a user-specified number of scales with a meshing
 *  algorithm (see `Scale_space_reconstruction_3::Mesher`).
 *
 *  The class stores the point set at the current scale and the
 *  reconstructed surface. User can either use the reconstructed
 *  surface at the current scale (smoothed version) or propagate the
 *  connectivity of the facets back to the original point set. This
 *  latest option makes it possible to reconstruct an interpolative
 *  surface on a noisy point set.
 *  
 *  \tparam Geom_traits is the geometric traits class. It must be a
 *  model of `DelaunayTriangulationTraits_3`. It must have a
 *  `RealEmbeddable` field number type. Generally,
 *  `Exact_predicates_inexact_constructions_kernel` is preferred.
 */
template <typename Geom_traits>
class Scale_space_surface_reconstruction_3
{
public:

  typedef typename Geom_traits::FT              FT;                              ///< defines the field number type.
  typedef typename Geom_traits::Point_3         Point;                           ///< defines the point type.
  typedef cpp11::array<std::size_t, 3> Facet;                           ///< defines a facet of the surface (triple of point indices).

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type                      Point_iterator;         ///< defines an iterator over the points.
  typedef const unspecified_type                Point_const_iterator;   ///< defines a constant iterator over the points.
#else
  typedef typename std::vector<Point>           Point_vector;
  typedef typename Point_vector::iterator       Point_iterator;
  typedef typename Point_vector::const_iterator Point_const_iterator;
#endif

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type                      Facet_iterator;         ///< defines an iterator over the points.
  typedef const unspecified_type                Facet_const_iterator;   ///< defines a constant iterator over the points.
#else
  typedef typename std::vector<Facet>           Facet_vector;
  typedef typename Facet_vector::iterator       Facet_iterator;
  typedef typename Facet_vector::const_iterator Facet_const_iterator;
#endif

  // Default algorithms used (same as in old API)
  typedef Scale_space_reconstruction_3::Weighted_PCA_smoother<Geom_traits> Weighted_PCA_smoother;
  typedef Scale_space_reconstruction_3::Alpha_shape_mesher<Geom_traits> Alpha_shape_mesher;
  
private:

  Point_vector m_points;
  Facet_vector m_facets;

  FT m_internal_squared_radius; // For backward compatibility

public:

  /// @{
  /// \name Initialization
  
  /**
   * Empty constructor.
   */
  Scale_space_surface_reconstruction_3 () : m_internal_squared_radius(0.) { }

  /**
   * Constructs a reconstruction object with the given point range.
   *
   * \note This constructor is equivalent to constructing an empty
   * object and then calling <code>[insert(begin, end)](\ref insert)</code>.
   *
   * \tparam InputIterator is an iterator over the point collection.
   *  The value type of the iterator must be a `Point`.
   *
   * \param begin is an iterator to the first point of the collection.
   * \param end is a past-the-end iterator for the point collection.
   */
  template <typename InputIterator>
  Scale_space_surface_reconstruction_3 (InputIterator begin, InputIterator end)
    : m_internal_squared_radius (0.)
  {
    insert (begin, end);
  }

  /// inserts a point into the scale-space at the current scale.
  /** \param p is the point to insert.
   *
   *  \note Inserting the point does not automatically construct or
   *  update the surface.
   *
   *  \note In order to construct the surface, call
   *  `#reconstruct_surface()`.
   *
   *  \sa `insert(InputIterator begin, InputIterator end)`.
   */
  void insert (const Point& p)
  {
    m_points.push_back (p);
  }
  
  /// inserts a collection of points into the scale-space at the current scale.
  /** \tparam InputIterator is an iterator over the point collection.
   *  The value type of the iterator must be a `Point`.
   *
   *  \param begin is an iterator to the first point of the collection.
   *  \param end is a past-the-end iterator for the point collection.
   *
   *  \note Inserting the points does not automatically construct or
   *  update the surface.
   *
   *  \note In order to construct the surface, call
   *  `reconstruct_surface()` after inserting the points.
   *
   *  \sa `insert(const Point& p)`.
   */
  template <typename InputIterator>
  void insert (InputIterator begin, InputIterator end)
  {
    std::copy (begin, end, std::back_inserter (m_points));
  }

  /// @}

  /// @{
  /// \name Algorithms
  
  /// increases the scale by a number of iterations.
  /** Each iteration the scale is increased, the points set at a higher scale
   *  is computed. At a higher scale, the points set is smoother.
   *
   *  \note If no smoothing algorithm is specified, the default
   *  `Weighted_PCA_smoother` is used.
   *
   *  \tparam Smoother model of `Scale_space_reconstruction_3::Smoother`.
   *
   *  \param iterations is the number of iterations to perform. If
   *  `iterations` is 0, nothing happens.
   *  \param smoother the smoother object.
   *
   *  \note This method processes the point set at the current scale. The
   *  points can be set with <code>[insert(begin, end)](\ref insert)</code>.
   *
   *  \note If the surface was already constructed, increasing the scale
   *  will not automatically adjust the surface.
   *
   *  \sa `reconstruct_surface()`.
   */
  template <typename Smoother>
  void increase_scale (std::size_t iterations = 1, const Smoother& smoother = Smoother())
  {
    for (std::size_t i = 0; i < iterations; ++ i)
      const_cast<Smoother&>(smoother) (m_points.begin(), m_points.end());
  }

  /// \cond SKIP_IN_MANUAL
  void increase_scale (std::size_t iterations = 1)
  {
    Weighted_PCA_smoother smoother;
    increase_scale (iterations, smoother);
    m_internal_squared_radius = smoother.squared_radius();
  }
  /// \endcond

  /// constructs a triangle mesh from the point set at a fixed scale.
  /** The order of the points at the current scale is the same as the order
   *  at the original scale, meaning that the surface can interpolate the
   *  point set at the original scale by applying the indices of the surface
   *  triangles to the original point set.
   *
   *  After construction, the facets of the surface can be iterated using
   *  `facets_begin()` and `facets_end()`.
   *
   *  \note If no meshing algorithm is specified, the default
   *  `Alpha_shape_mesher` is used.
   *
   *  \tparam Mesher model of `Scale_space_reconstruction_3::Mesher`.
   *
   *  \param mesher the mesher object.
   *
   *  \note This method processes the point set at the current scale. The
   *  points can be set with <code>[insert(begin, end)](\ref insert)</code>.
   *
   *  \sa `increase_scale()`.
   */
  template <typename Mesher>
  void reconstruct_surface (const Mesher& mesher = Mesher())
  {
    m_facets.clear();
    const_cast<Mesher&>(mesher) (m_points.begin(), m_points.end(), std::back_inserter (m_facets));
  }

  /// \cond SKIP_IN_MANUAL
  void reconstruct_surface ()
  {
    CGAL_assertion (m_internal_squared_radius != 0.);
    reconstruct_surface (Alpha_shape_mesher(m_internal_squared_radius));
  }
  /// \endcond

  /// @}

  /// @{
  /// \name Output

  /// gives the number of points of the surface.
  std::size_t number_of_points() const { return m_points.size(); }
  
  /// gives an iterator to the first point at the current scale.
  /** \warning Changes to the scale-space do not cause an automatic update to
   *  the surface.
   */
  Point_iterator points_begin() { return m_points.begin(); }

  /// gives a past-the-end iterator of the points at the current scale.
  /** \warning Changes to the scale-space do not cause an automatic update to
   *  the surface.
   */
  Point_iterator points_end() { return m_points.end(); }
  
  /// gives an iterator to the first point at the current scale.
  Point_const_iterator points_begin() const { return m_points.begin(); }
  
  /// gives a past-the-end iterator of the points at the current scale.
  Point_const_iterator points_end() const { return m_points.end(); }

  /// gives the number of facets of the surface.
  std::size_t number_of_facets() const { return m_facets.size(); }

  /// gives an iterator to the first triple in the surface.
  /** \warning Changes to the surface may change its topology.
   */
  Facet_iterator facets_begin() { return m_facets.begin(); }

  /// gives a past-the-end iterator of the triples in the surface.
  /** \warning Changes to the surface may change its topology.
   */
  Facet_iterator facets_end() { return m_facets.end(); }

  /// gives an iterator to the first triple in the surface.
  Facet_const_iterator facets_begin() const { return m_facets.begin(); }

  /// gives a past-the-end iterator of the triples in the surface.
  Facet_const_iterator facets_end() const { return m_facets.end(); }

  /// @}
};


template <typename Geom_traits>
std::ostream& operator<< (std::ostream& os, const CGAL::Scale_space_surface_reconstruction_3<Geom_traits>& ss)
{
  typedef CGAL::Scale_space_surface_reconstruction_3<Geom_traits> Reconstruction;
  
  os << "OFF" << std::endl
     << ss.number_of_points() << " " << ss.number_of_facets() << " 0" << std::endl;
  for (typename Reconstruction::Point_const_iterator it = ss.points_begin();
       it != ss.points_end(); ++ it)
    os << *it << std::endl;
  for (typename Reconstruction::Facet_const_iterator it = ss.facets_begin();
       it != ss.facets_end(); ++ it)
    os << "3 " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << std::endl;
  return os;
}

} // namespace CGAL

template< typename T >
std::ostream&
operator<<( std::ostream& os, const CGAL::cpp11::array< T, 3 >& t ) {
    return os << t[0] << " " << t[1] << " " << t[2];
}

template< typename T >
std::istream&
operator>>( std::istream& is, CGAL::cpp11::array< T, 3 >& t ) {
    return is >> get<0>(t) >> get<1>(t) >> get<2>(t);
}

#endif // CGAL_SCALE_SPACE_SURFACE_RECONSTRUCTION_3_H

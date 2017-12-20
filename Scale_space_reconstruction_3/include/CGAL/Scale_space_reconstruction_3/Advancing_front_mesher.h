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
// Author(s):      Simon Giraudot


#ifndef CGAL_SCALE_SPACE_RECONSTRUCTION_3_ADVANCING_FRONT_MESHER_H
#define CGAL_SCALE_SPACE_RECONSTRUCTION_3_ADVANCING_FRONT_MESHER_H

#include <CGAL/license/Scale_space_reconstruction_3.h>

#include <CGAL/Advancing_front_surface_reconstruction.h>

#include <CGAL/Union_find.h>

namespace CGAL
{

namespace Scale_space_reconstruction_3
{
  
/** \ingroup PkgScaleSpaceReconstruction3Classes
 *
 *  Surface mesher for scale space reconstruction based on
 *  `CGAL::Advancing_front_surface_reconstruction`.
 *
 *  This class applies the advancing front reconstruction algorithm
 *  with the possibility of using an upper bound on the length of the
 *  produced facets.
 * 
 *  \cgalModels CGAL::Scale_space_reconstruction_3::Mesher
 *
 *  \tparam Geom_traits geometric traits class. It must be a
 *  model of `DelaunayTriangulationTraits_3`. It must have a
 *  `RealEmbeddable` field number type. Generally,
 *  `Exact_predicates_inexact_constructions_kernel` is preferred.
 */
template <typename Geom_traits>
class Advancing_front_mesher
{
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3                        Point;          ///< defines the point type.
  
  typedef CGAL::cpp11::array< unsigned int, 3 >       Facet;
private:

  class Priority
  {
    double bound;
  public:
    Priority (double bound)
      : bound(bound)
    {}

    template <typename AdvancingFront, typename Cell_handle>
    double operator() (const AdvancingFront& adv, Cell_handle& c,
                       const int& index) const
    {
      if(bound == 0){
        return adv.smallest_radius_delaunay_sphere (c, index);
      }

      double d = 0.;
      d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                c->vertex((index+2)%4)->point()));
      if(d>bound) return adv.infinity();
      d = sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                                c->vertex((index+3)%4)->point()));
      if(d>bound) return adv.infinity();
      d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                                c->vertex((index+3)%4)->point()));
      if(d>bound) return adv.infinity();

      return adv.smallest_radius_delaunay_sphere (c, index);
    }
  };

  Priority m_priority;
  FT m_radius_ratio_bound;
  FT m_beta;
  
public:

  /**
   * Constructs and advancing front mesher.
   *
   * \param maximum_facet_length upper bound on the length of the facets.
   * \param radius_ratio_bound candidates incident to surface
         triangles which are not in the beta-wedge are discarded, if
         the ratio of their radius and the radius of the surface
         triangle is larger than `radius_ratio_bound`.  Described in
         Section \ref AFSR_Boundaries
   * \param beta half the angle of the wedge in which only the radius
         of triangles counts for the plausibility of candidates.
         Described in Section \ref AFSR_Selection
   */
  Advancing_front_mesher (FT maximum_facet_length = 0.,
                          FT radius_ratio_bound = 5,
                          FT beta = 0.52)
    : m_priority (maximum_facet_length), m_radius_ratio_bound (radius_ratio_bound), m_beta (beta)
  {

  }

  /// \cond SKIP_IN_MANUAL
  template <typename InputIterator, typename OutputIterator>
  void operator() (InputIterator begin, InputIterator end, OutputIterator output)
  {
    CGAL::advancing_front_surface_reconstruction (begin, end, output,
                                                  m_priority,
                                                  m_radius_ratio_bound,
                                                  m_beta);
  }
  /// \endcond

};

  
} // namespace Scale_space_reconstruction_3

} // namespace CGAL

#endif // CGAL_SCALE_SPACE_RECONSTRUCTION_3_ADVANCING_FRONT_MESHER_H

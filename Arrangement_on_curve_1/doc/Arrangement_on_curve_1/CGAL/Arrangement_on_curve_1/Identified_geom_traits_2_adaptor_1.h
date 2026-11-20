// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TraitsClasses
 *
 * \brief `Identified_geom_traits_2_adaptor_1` adapts a 2D geometry traits
 * class of the 2D Arrangements package, whose parameter space has identified
 * left and right sides, into a 1D geometry traits class whose parameter space
 * is identified.
 *
 * The points of the 1D arrangement are points of the 2D traits that lie on a
 * single `x`-monotone master curve whose two ends lie on the 2D
 * identification curve: its minimal end on the left boundary and its maximal
 * end on the right boundary of the 2D parameter space.  Both ends are the same
 * point of the surface, namely the identification point of the 1D
 * arrangement.  A point is the identification point if and only if it lies on
 * the 2D identification curve.  The remaining points are ordered
 * lexicographically, by `x` and then by `y`, which is the order along the
 * master curve.  Pair this traits with `Identified_topology_traits`.
 *
 * An example of a suitable 2D traits is
 * `Arr_geodesic_arc_on_sphere_traits_2`, where the 2D identification curve is
 * the meridian that separates the left and right sides of the parameter space.
 *
 * \tparam GeomTraits_2 a model of both `AosBasicTraits_2` and
 *   `AosIdentifiedVerticalTraits_2`, whose `Left_side_category` and
 *   `Right_side_category` are both `Arr_identified_side_tag`.
 *
 * \cgalModels{AocIdentifiedBoundaryTraits_1}
 *
 * \sa `Open_geom_traits_2_adaptor_1<GeomTraits_2>`
 * \sa `Left_closed_geom_traits_2_adaptor_1<GeomTraits_2>`
 * \sa `Right_closed_geom_traits_2_adaptor_1<GeomTraits_2>`
 * \sa `Closed_geom_traits_2_adaptor_1<GeomTraits_2>`
 */
template <typename GeomTraits_2>
class Identified_geom_traits_2_adaptor_1 {
public:
  //! \name Types
  //! @{

  //! the adapted 2D geometry traits type.
  typedef GeomTraits_2 Geom_traits_2;

  //! the point type; a 2D point that lies on the master curve.
  typedef typename Geom_traits_2::Point_2 Point_1;

  //! a shared-ownership pointer to a const 2D geometry traits object.
  typedef std::shared_ptr<const Geom_traits_2> Shared_geom_traits_2;

  //! @}

  //! \name Creation
  //! @{

  /*! constructs from a default-constructed 2D geometry traits object.
   */
  Identified_geom_traits_2_adaptor_1();

  /*! constructs from a shared 2D geometry traits object.
   * \param traits_2 a shared pointer to the 2D geometry traits.
   */
  explicit Identified_geom_traits_2_adaptor_1(Shared_geom_traits_2 traits_2);

  /*! constructs from a copy of a 2D geometry traits object.
   * \param traits_2 the 2D geometry traits to copy.
   */
  explicit Identified_geom_traits_2_adaptor_1(const Geom_traits_2& traits_2);

  //! @}

  //! \name Accessors
  //! @{

  /*! obtains the adapted 2D geometry traits object.
   */
  const Geom_traits_2& geom_traits_2() const;

  /*! obtains the shared-ownership pointer to the adapted 2D geometry traits object.
   */
  const Shared_geom_traits_2& shared_geom_traits_2() const;

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

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
 * \brief `Right_closed_geom_traits_2_adaptor_1` adapts a 2D geometry traits
 * class of the 2D Arrangements package into a 1D geometry traits class whose
 * parameter space is open on the left side and closed on the right side.
 *
 * The points of the 1D arrangement are points of the 2D traits that lie on a
 * single master curve, which is either an `x`-monotone curve or a vertical
 * segment.  The right boundary point, provided at construction, is the
 * lexicographically largest point of the master curve; the left endpoint of
 * the master curve is not a point of the 1D arrangement.  Points are ordered
 * lexicographically, by `x` and then by `y`.  Pair this traits with
 * `Right_closed_topology_traits`, which `Default_topology_traits` selects for
 * it.
 *
 * \tparam GeomTraits_2 a model of `AosBasicTraits_2`.
 *
 * \cgalModels{AocClosedBoundaryTraits_1}
 *
 * \sa `Open_geom_traits_2_adaptor_1<GeomTraits_2>`
 * \sa `Left_closed_geom_traits_2_adaptor_1<GeomTraits_2>`
 * \sa `Closed_geom_traits_2_adaptor_1<GeomTraits_2>`
 * \sa `Identified_geom_traits_2_adaptor_1<GeomTraits_2>`
 */
template <typename GeomTraits_2>
class Right_closed_geom_traits_2_adaptor_1 {
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

  /*! constructs from a shared 2D geometry traits object and the right boundary point.
   * \param traits_2 a shared pointer to the 2D geometry traits.
   * \param right_bound the right boundary point.
   */
  Right_closed_geom_traits_2_adaptor_1(Shared_geom_traits_2 traits_2, const Point_1& right_bound);

  /*! constructs from a copy of a 2D geometry traits object and the right boundary point.
   * \param traits_2 the 2D geometry traits to copy.
   * \param right_bound the right boundary point.
   */
  Right_closed_geom_traits_2_adaptor_1(const Geom_traits_2& traits_2, const Point_1& right_bound);

  //! @}

  //! \name Accessors
  //! @{

  /*! obtains the adapted 2D geometry traits object.
   */
  const Geom_traits_2& geom_traits_2() const;

  /*! obtains the shared-ownership pointer to the adapted 2D geometry traits object.
   */
  const Shared_geom_traits_2& shared_geom_traits_2() const;

  /*! obtains the right boundary point.
   */
  const Point_1& right_bound() const;

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

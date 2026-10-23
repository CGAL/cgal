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
 * \brief `Closed_geom_traits_2_adaptor_1` adapts a 2D geometry traits class
 * of the 2D Arrangements package into a 1D geometry traits class whose
 * parameter space is closed on both sides.
 *
 * The points of the 1D arrangement are points of the 2D traits that lie on a
 * single master curve, which is either an `x`-monotone curve or a vertical
 * segment.  The left and right boundary points, provided at construction, are
 * the lexicographically smallest and largest points of the master curve,
 * respectively.  Points are ordered lexicographically, by `x` and then by
 * `y`.  Pair this traits with `Closed_topology_traits`.
 *
 * \tparam GeomTraits_2 a model of `AosBasicTraits_2`.
 *
 * \cgalModels{AocClosedBoundaryTraits_1}
 *
 * \sa `Open_geom_traits_2_adaptor_1<GeomTraits_2>`
 * \sa `Left_closed_geom_traits_2_adaptor_1<GeomTraits_2>`
 * \sa `Right_closed_geom_traits_2_adaptor_1<GeomTraits_2>`
 * \sa `Identified_geom_traits_2_adaptor_1<GeomTraits_2>`
 */
template <typename GeomTraits_2>
class Closed_geom_traits_2_adaptor_1 {
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

  /*! constructs from a shared 2D geometry traits object and the two boundary points.
   * \param traits_2 a shared pointer to the 2D geometry traits.
   * \param left_bound the left boundary point.
   * \param right_bound the right boundary point.
   * \pre `left_bound` is lexicographically smaller than `right_bound`.
   */
  Closed_geom_traits_2_adaptor_1(Shared_geom_traits_2 traits_2, const Point_1& left_bound,
                                 const Point_1& right_bound);

  /*! constructs from a copy of a 2D geometry traits object and the two boundary points.
   * \param traits_2 the 2D geometry traits to copy.
   * \param left_bound the left boundary point.
   * \param right_bound the right boundary point.
   * \pre `left_bound` is lexicographically smaller than `right_bound`.
   */
  Closed_geom_traits_2_adaptor_1(const Geom_traits_2& traits_2, const Point_1& left_bound,
                                 const Point_1& right_bound);

  //! @}

  //! \name Accessors
  //! @{

  /*! obtains the adapted 2D geometry traits object.
   */
  const Geom_traits_2& geom_traits_2() const;

  /*! obtains the shared-ownership pointer to the adapted 2D geometry traits object.
   */
  const Shared_geom_traits_2& shared_geom_traits_2() const;

  /*! obtains the left boundary point.
   */
  const Point_1& left_bound() const;

  /*! obtains the right boundary point.
   */
  const Point_1& right_bound() const;

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

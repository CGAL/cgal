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
 * \brief `Direction_geom_traits_2_adaptor_1` adapts a 2D geometry traits
 * class that provides directions in the plane into a 1D geometry traits class
 * of the circle of directions, whose parameter space is identified.
 *
 * The points of the 1D arrangement are directions in the plane. The
 * identification point is a given reference direction \f$s\f$, and the
 * directions are ordered counterclockwise, starting right after \f$s\f$:
 * \f$s < d_1 < d_2 < \ldots < d_n < s\f$ (one full turn). Pair this traits
 * with `Identified_topology_traits`.
 *
 * The order is determined by the predicate
 * `Kernel::CounterclockwiseInBetween_2` alone; unlike `Circle_traits_1`, it
 * refers to no axis, and no direction is constructed. A direction \f$a\f$
 * precedes a direction \f$b\f$ if and only if \f$a\f$ is counterclockwise
 * strictly in between \f$s\f$ and \f$b\f$; \f$a\f$ and \f$b\f$ are equal if
 * and only if neither precedes the other. A direction \f$a\f$ is the
 * identification point if and only if it is equal to \f$s\f$, that is, if
 * and only if \f$a\f$ is not counterclockwise strictly in between \f$s\f$ and
 * \f$s\f$ (the open range from \f$s\f$ to \f$s\f$ is the entire circle except
 * for \f$s\f$).
 *
 * \tparam GeomTraits_2 a geometry traits class that provides the type
 *   `Direction_2`, a model of `Kernel::Direction_2`, and the functor type
 *   `Counterclockwise_in_between_2`, a model of
 *   `Kernel::CounterclockwiseInBetween_2`, along with the member function
 *   `counterclockwise_in_between_2_object()`. Every model of `Kernel` is
 *   such a class, and so are the geometry traits classes derived from a
 *   kernel, such as `Arr_segment_traits_2<Kernel>`.
 *
 * \cgalModels{AocIdentifiedBoundaryTraits_1}
 *
 * \sa `Circle_traits_1<Kernel>`
 * \sa `Identified_geom_traits_2_adaptor_1<GeomTraits_2>`
 */
template <typename GeomTraits_2>
class Direction_geom_traits_2_adaptor_1 {
public:
  //! \name Types
  //! @{

  //! the adapted 2D geometry traits type.
  typedef GeomTraits_2 Geom_traits_2;

  //! the point type; a direction in the plane.
  typedef typename Geom_traits_2::Direction_2 Point_1;

  //! a shared-ownership pointer to a const 2D geometry traits object.
  typedef std::shared_ptr<const Geom_traits_2> Shared_geom_traits_2;

  //! @}

  //! \name Functor Types
  //! @{

  /*! a functor that determines whether two directions are equal. It
   * provides the operator `bool operator()(const Point_1& d1, const Point_1& d2) const`.
   */
  typedef unspecified_type Equal_1;

  //! @}

  //! \name Creation
  //! @{

  /*! constructs from a shared 2D geometry traits object and a reference
   * direction.
   * \param traits_2 a shared pointer to the 2D geometry traits.
   * \param reference the reference direction, which is the identification
   *        point.
   * \pre `traits_2` is not null.
   */
  Direction_geom_traits_2_adaptor_1(Shared_geom_traits_2 traits_2, const Point_1& reference);

  /*! constructs from a copy of a 2D geometry traits object and a reference
   * direction.
   * \param traits_2 the 2D geometry traits to copy.
   * \param reference the reference direction, which is the identification
   *        point.
   */
  Direction_geom_traits_2_adaptor_1(const Geom_traits_2& traits_2, const Point_1& reference);

  //! @}

  //! \name Accessors
  //! @{

  /*! obtains the adapted 2D geometry traits object.
   */
  const Geom_traits_2& geom_traits_2() const;

  /*! obtains the shared-ownership pointer to the adapted 2D geometry traits object.
   */
  const Shared_geom_traits_2& shared_geom_traits_2() const;

  /*! obtains the reference direction, which is the identification point.
   */
  const Point_1& reference_direction() const;

  /*! obtains an `Equal_1` functor.
   */
  Equal_1 equal_1_object() const;

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

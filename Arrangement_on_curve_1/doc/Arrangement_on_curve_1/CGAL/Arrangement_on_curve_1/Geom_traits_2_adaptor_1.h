// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/*! \ingroup PkgArrangementOnCurve1TraitsClasses
 *
 * \brief `Geom_traits_2_adaptor_1` wraps any 2D geometry traits class that
 * models the `ArrangementTraits_2` concept (or a refinement thereof) into a 1D
 * geometry traits class suitable for use with `Arrangement_on_curve_1`.
 *
 * All points of the 1D arrangement are `Point_2` objects from the wrapped 2D
 * traits, all assumed to lie on a single master arc which is either an
 * x-monotone curve or a single vertical segment.
 *
 * The adaptor models an open (unbounded) parameter space — both ends of the
 * master arc extend to infinity — and should therefore be paired with
 * `Unbounded_topology_traits`.  Subclasses may override the side categories to
 * handle bounded arcs.
 *
 * \tparam GeomTraits_2 a model of `ArrangementTraits_2` (or a refinement).
 *
 * \cgalModels{AocTraits_1}
 */
template <typename GeomTraits_2>
class Geom_traits_2_adaptor_1 {
public:
  //! \name Types
  //! @{

  //! The wrapped 2D geometry traits type.
  typedef GeomTraits_2 Geom_traits_2;

  //! The 1D arrangement point type — a 2D point lying on the master arc.
  typedef typename Geom_traits_2::Point_2 Point_1;

  //! A 2D curve type used to define the master arc.
  typedef typename Geom_traits_2::Curve_2 Curve_2;

  //! An x-monotone 2D curve type used to define the master arc.
  typedef typename Geom_traits_2::X_monotone_curve_2 X_monotone_curve_2;

  //! Shared-ownership pointer to the 2D traits (allows sharing a single instance).
  typedef std::shared_ptr<const Geom_traits_2> Shared_geom_traits_2;

  //! Both sides of the arc are open (unbounded).
  typedef Open_side_tag Left_side_category;

  //! Both sides of the arc are open (unbounded).
  typedef Open_side_tag Right_side_category;

  //! @}

  //! \name Creation
  //! @{

  /*! default constructor; allocates a default-constructed `Geom_traits_2`.
   */
  Geom_traits_2_adaptor_1();

  /*! constructs from a shared pointer to the 2D geometry traits.
   * \param shared_traits_2 a shared pointer to the underlying 2D traits.
   */
  Geom_traits_2_adaptor_1(Shared_geom_traits_2 shared_traits_2);

  /*! constructs by copying the given 2D traits object.
   * \param traits_2 the underlying 2D traits to copy.
   */
  explicit Geom_traits_2_adaptor_1(const Geom_traits_2& traits_2);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

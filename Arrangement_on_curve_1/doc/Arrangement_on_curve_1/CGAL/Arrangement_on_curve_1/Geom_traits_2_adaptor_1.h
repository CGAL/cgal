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
 * adapts any 2D geometry traits class of the 2D Arrangements package (i.e., any
 * model of the `ArrangementTraits_2` concept, or a refinement thereof) into a 1D
 * geometry traits class, suitable for use with `Arrangement_on_curve_1`. The
 * points of the 1D arrangement are the `Point_2` objects of the 2D traits, all
 * assumed to lie on a single curve, the "master curve", which is either
 * \f$x\f$-monotone or a single vertical segment.
 *
 * \cgalModels{AocTraits_1}
 */
template <typename GeomTraits_2>
class Geom_traits_2_adaptor_1 {
public:
  /// \name Types
  /// @{
  typedef GeomTraits_2 Geom_traits_2;                                    ///< The 2D geometry traits type.
  typedef typename Geom_traits_2::Point_2 Point_1;                       ///< The 1D arrangement point type.
  typedef typename Geom_traits_2::X_monotone_curve_2 X_monotone_curve_2; ///< The 2D master curve type.
  typedef std::shared_ptr<const Geom_traits_2> Shared_geom_traits_2;     ///< A smart pointer to the 2D traits type.
  /// @}

  /// \name Creation
  /// @{

  /*! constructs a traits object that adapts a given 2D geometry traits.
   * \param shared_traits_2 a smart pointer to the 2D geometry traits
   */
  Geom_traits_2_adaptor_1(Shared_geom_traits_2 shared_traits_2);

  /// @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/** \ingroup PkgArrangementOnCurve1TraitsClasses
 *
 * \brief `Ft_traits_1` is a minimal, scalar model of the `AocTraits_1` concept.
 *
 * It parameterizes a 1D arrangement directly using a coordinate field type
 * (\c Kernel::FT) representing a continuous horizontal number line \f$(-\infty, +\infty)\f$.
 *
 * \tparam Kernel must be a standard \cgal `Kernel` model (e.g., `Exact_predicates_exact_constructions_kernel`).
 *
 * \cgalModels{AocTraits_1}
 */
template <typename Kernel>
class Ft_traits_1 {
public:
  /// \name Types
  /// @{

  /**
   * \brief The 1D geometric point type representation, mapping directly to the underlying field type scalar.
   */
  typedef typename Kernel::FT Point_1;
  /// @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

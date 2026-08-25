// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/** \ingroup PkgArrangementOnCurve1TraitsClasses
 *
 * \brief `Line_d_traits_1` is a geometric traits model of `AocTraits_1` designed for
 * subdivisions along an infinite supporting line embedded in a d-dimensional coordinate space.
 *
 * Points are represented using multi-dimensional dynamic coordinates (\c Kernel::Point_d).
 * The linear ordering is resolved by computing inner dot products across coordinate dimensions via
 * Cartesian iterator ranges.
 *
 * \tparam Kernel must be a standard multi-dimensional CGAL Kernel model (e.g., `Cartesian_d`).
 *
 * \cgalModels{AocTraits_1}
 */
template <typename Kernel>
class Line_d_traits_1 {
public:
  /// \name Types
  /// @{
  typedef typename Kernel::Point_d Point_1;   ///< The 1D arrangement point represented as a d-dimensional coordinate.
  typedef typename Kernel::Line_d Line_d;     ///< The d-dimensional master line supporting the topology.
  typedef typename Kernel::Vector_d Vector_d; ///< A d-dimensional direction vector.
  typedef typename Kernel::FT FT;             ///< The core algebraic field type.
  /// @}

  /// \name Creation
  /// @{

  /**
   * \brief Constructor from an embedding d-dimensional line space.
   * \param line The supporting geometric master line trajectory.
   */
  Line_d_traits_1(const Line_d& line);
  /// @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/**
 * \ingroup PkgArrangementOnCurve1Ref
 *
 * \brief `Line_2_traits_1` is a geometric traits model of `AocTraits_1` designed for
 * subdivisions along an infinite supporting line embedded in 2D Euclidean space.
 *
 * Points are represented using standard 2D Cartesian coordinates, and their ordering
 * is determined by projecting them onto the direction vector of the master line.
 *
 * \tparam Kernel must be a standard CGAL Kernel model.
 *
 * \cgalModels{AocTraits_1}
 */
template <typename Kernel>
class Line_2_traits_1 {
public:
  /// \name Types
  /// @{
  typedef typename Kernel::Point_2 Point_1;   ///< The 1D arrangement point represented as a 2D coordinate.
  typedef typename Kernel::Line_2 Line_2;     ///< The 2D master line supporting the topology.
  typedef typename Kernel::Vector_2 Vector_2; ///< A 2D structural direction vector.
  /// @}

  /// \name Creation
  /// @{

  /**
   * \brief Constructor from an embedding 2D line.
   * \param line The supporting geometric master line trajectory.
   */
  Line_2_traits_1(const Line_2& line);
  /// @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

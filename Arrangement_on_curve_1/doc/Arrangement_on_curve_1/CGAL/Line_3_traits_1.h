// Copyright (c) 2026 Tel-Aviv University (Israel). All rights reserved.
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace CGAL {
namespace Arrangement_on_curve_1 {

/**
 * \ingroup PkgArrangementOnCurve1Ref
 *
 * \brief `Line_3_traits_1` is a geometric traits model of `AocTraits_1` designed for
 * subdivisions along an infinite supporting line embedded in 3D Euclidean space.
 *
 * Points are represented using standard 3D Cartesian coordinates, and their relative sequence
 * is computed by mapping positions across the directional trajectory vector of the 3D line.
 *
 * \tparam Kernel must be a standard CGAL Kernel model.
 *
 * \cgalModels{AocTraits_1}
 */
template <typename Kernel>
class Line_3_traits_1 {
public:
  /// \name Types
  /// @{
  typedef typename Kernel::Point_3 Point_1;   ///< The 1D arrangement point represented as a 3D coordinate.
  typedef typename Kernel::Line_3 Line_3;     ///< The 3D master line supporting the topology.
  typedef typename Kernel::Vector_3 Vector_3; ///< A 3D structural direction vector.
  /// @}

  /// \name Creation
  /// @{

  /**
   * \brief Constructor from an embedding 3D line.
   * \param line The supporting geometric master line trajectory.
   */
  Line_3_traits_1(const Line_3& line);
  /// @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

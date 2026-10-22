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
 * \brief `Line_3_traits_1` is a geometric traits model of `AocTraits_1` designed
 * for subdivisions along an infinite supporting line embedded in 3D Euclidean
 * space.
 *
 * Points are represented using standard 3D Cartesian coordinates, and their
 * relative sequence is computed by projecting positions onto the directional
 * trajectory vector of the 3D line.  Both sides of the parameter space are
 * open, so `Line_3_traits_1` should be paired with `Open_topology_traits`.
 *
 * \tparam Kernel must be a standard \cgal `Kernel` model.
 *
 * \cgalModels{AocTraits_1}
 */
template <typename Kernel>
class Line_3_traits_1 {
public:
  //! \name Types
  //! @{

  //! The 1D arrangement point type is a 3D coordinate on the supporting line.
  typedef typename Kernel::Point_3 Point_1;

  //! The 3D master line supporting the topology.
  typedef typename Kernel::Line_3 Line_3;

  //! A 3D direction vector along the master line.
  typedef typename Kernel::Vector_3 Vector_3;

  //! The left side of the parameter space is open.
  typedef Open_side_tag Left_side_category;

  //! The right side of the parameter space is open.
  typedef Open_side_tag Right_side_category;

  //! @}

  //! \name Creation
  //! @{

  /*! constructs from an embedding 3D line.
   * \param line the supporting geometric master line trajectory.
   */
  Line_3_traits_1(const Line_3& line);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

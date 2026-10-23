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
 * \brief `Line_d_traits_1` is a geometric traits model of `AocTraits_1` designed
 * for subdivisions along an infinite supporting line embedded in a
 * \f$d\f$-dimensional coordinate space.
 *
 * Points are represented using multi-dimensional dynamic coordinates
 * (`Kernel::Point_d`).  The linear ordering is resolved by computing inner dot
 * products across coordinate dimensions via Cartesian iterator ranges.  Both
 * sides of the parameter space are open, so `Line_d_traits_1` should be paired
 * with `Open_topology_traits`.
 *
 * \tparam Kernel must be a standard multi-dimensional \cgal `Kernel` model
 *   (e.g., `Cartesian_d`).
 *
 * \cgalModels{AocTraits_1}
 */
template <typename Kernel>
class Line_d_traits_1 {
public:
  //! \name Types
  //! @{

  //! The 1D arrangement point type is a \f$d\f$-dimensional coordinate.
  typedef typename Kernel::Point_d Point_1;

  //! The \f$d\f$-dimensional master line supporting the topology.
  typedef typename Kernel::Line_d Line_d;

  //! A \f$d\f$-dimensional direction vector along the master line.
  typedef typename Kernel::Vector_d Vector_d;

  //! The core algebraic field type.
  typedef typename Kernel::FT FT;

  //! The left side of the parameter space is open.
  typedef Open_side_tag Left_side_category;

  //! The right side of the parameter space is open.
  typedef Open_side_tag Right_side_category;

  //! @}

  //! \name Creation
  //! @{

  /*! constructs from an embedding \f$d\f$-dimensional line.
   * \param line the supporting geometric master line trajectory.
   */
  Line_d_traits_1(const Line_d& line);

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

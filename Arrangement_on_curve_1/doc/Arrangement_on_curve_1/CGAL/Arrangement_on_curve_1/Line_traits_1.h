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
 * \brief `Line_traits_1` is a minimal, scalar model of the `AocTraits_1` concept.
 *
 * It parameterizes a 1D arrangement using a coordinate field type (\c Kernel::FT)
 * representing a continuous open interval \f$(-\infty, +\infty)\f$.
 * Both sides of the parameter space are open
 *
 * \tparam Kernel must be a standard \cgal `Kernel` model (e.g., `Exact_predicates_exact_constructions_kernel`).
 *
 * \cgalRefines{DefaultConstructible}

 * \cgalModels{AocTraits_1}
 */
template <typename Kernel>
class Line_traits_1 {
public:
  //! \name Types
  //! @{

  //! The 1D arrangement point type.
  using Point_1 = typename Kernel::FT;

  //! The left side of the parameter space is open.
  typedef Open_side_tag Left_side_category;

  //! The right side of the parameter space is open.
  typedef Open_side_tag Right_side_category;

  //! @}
};

} // namespace Arrangement_on_curve_1
} // namespace CGAL

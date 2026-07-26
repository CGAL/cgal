// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

#ifndef CGAL_AOC_TRAITS_1_H
#define CGAL_AOC_TRAITS_1_H

#include <CGAL/enum.h>

namespace CGAL {

/*! \ingroup PkgArrangementOnCurve1Concepts
 * \cgalConcept
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Ft_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Line_2_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Line_3_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Line_d_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Geom_traits_2_adaptor_1<GeomTraits_2>}
 * \cgalHasModelsEnd
 */
class AocTraits_1 {
public:
  /// A point guaranteed to lie on the master curve
  typedef unspecified_type Point_1;

  /// Functor evaluating relative positions along the curve parameterization
  class Compare_x_1 {
  public:
    Comparison_result operator()(const Point_1& p1, const Point_1& p2) const;
  };

  Compare_position_1 compare_x_1_object() const;
};

} // namespace CGAL

#endif

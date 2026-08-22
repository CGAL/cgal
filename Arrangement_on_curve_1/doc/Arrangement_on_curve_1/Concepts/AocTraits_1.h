// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

/*! \ingroup PkgArrangementOnCurve1ConceptsTraits
 * \cgalConcept
 *
 * The concept `AocTraits_1` defines the minimal set of geometric
 * predicates needed for the construction and maintenance of objects of the
 * class `Arrangement_on_curve_1`, as well as performing simple queries (such as
 * point-location queries) on such arrangements.
 *
 * \cgalRefines{CopyConstructible,Assignable,DefaultConstructible}
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
  /// \name Types
  /// @{

  /// models the concept AocTraits::Point_1.
  typedef unspecified_type Point_1;

  /// @}

  /// \name Functor Types
  /// @{

  /// models the concept AocTraits::CompareX_1.
  typedef unspecified_type Compare_x_1;

  /// @}

  /// \name Accessing Functor Objects
  /// @{

  Compare_x_1 compare_x_1_object() const;

  /// @}
};

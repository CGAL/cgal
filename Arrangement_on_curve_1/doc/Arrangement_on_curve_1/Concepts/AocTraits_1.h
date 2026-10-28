// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

/*! \ingroup PkgArrangementOnCurve1ConceptsTraits
 * \cgalConcept
 *
 * The concept `AocTraits_1` defines the minimal set of geometric predicates
 * needed for the construction and maintenance of objects of the class
 * `CGAL::Arrangement_on_curve_1::Arrangement_on_curve_1`, as well as performing
 * simple queries (such as point-location queries) on such arrangements.
 *
 * \cgalRefines{CopyConstructible, Assignable}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Line_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Line_2_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Line_3_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Line_d_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Circle_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Open_geom_traits_2_adaptor_1<GeomTraits_2>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Closed_geom_traits_2_adaptor_1<GeomTraits_2>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Left_closed_geom_traits_2_adaptor_1<GeomTraits_2>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Right_closed_geom_traits_2_adaptor_1<GeomTraits_2>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Identified_geom_traits_2_adaptor_1<GeomTraits_2>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Direction_geom_traits_2_adaptor_1<GeomTraits_2>}
 * \cgalHasModelsEnd
 */
class AocTraits_1 {
public:
  //! \name Types
  //! @{

  //! models the concept `AocTraits::Point_1`.
  typedef unspecified_type Point_1;

  //! Describes the topological behaviour of the left boundary of the parameter
  //! space.  Must be one of `CGAL::Arrangement_on_curve_1::Open_side_tag`,
  //! `CGAL::Arrangement_on_curve_1::Closed_side_tag`, or
  //! `CGAL::Arrangement_on_curve_1::Identified_side_tag`.
  typedef unspecified_type Left_side_category;

  //! Describes the topological behaviour of the right boundary of the parameter
  //! space.  Must be one of `CGAL::Arrangement_on_curve_1::Open_side_tag`,
  //! `CGAL::Arrangement_on_curve_1::Closed_side_tag`, or
  //! `CGAL::Arrangement_on_curve_1::Identified_side_tag`.
  typedef unspecified_type Right_side_category;

  //! @}

  //! \name Functor Types
  //! @{

  //! models the concept `AocTraits::CompareX_1`.
  typedef unspecified_type Compare_x_1;

  //! @}

  //! \name Accessing Functor Objects
  //! @{

  //! returns the `Compare_x_1` functor.
  Compare_x_1 compare_x_1_object() const;

  //! @}
};

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
 * Models of the concept `AocIdentifiedBoundaryTraits_1` handle curves on a circle
 * or a curve that is topological equivalent to a circle. The circle is oriented
 * in such a way that the boundary of the parameter space, the circle is the
 * mapping of which, is identified on the left and right sides.
 *
 * \cgalRefines{AocTraits_1}
 *
 * When both categories are `CGAL::Arrangement_on_curve_1::Identified_side_tag`, the traits must also provide
 * the functor type `Is_on_x_identification_1` (and its accessor), which tests
 * whether a point lies on the identification point.
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Circle_traits_1<Kernel>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Identified_geom_traits_2_adaptor_1<GeomTraits_2>}
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Direction_geom_traits_2_adaptor_1<GeomTraits_2>}
 * \cgalHasModelsEnd
 */
class AocIdentifiedBoundaryTraits_1 {
public:
  //! \name Functor Types
  //! @{

  /*! Functor that tests whether a point lies exactly on the identification
   * point (which identifies the left and right boundaries of the parameter
   * space.
   *
   * \pre `Left_side_category` and `Right_side_category` are both
   *      `CGAL::Arrangement_on_curve_1::Identified_side_tag`.
   */
  typedef unspecified_type Is_on_x_identification_1;

  //! @}

  //! \name Accessing Functor Objects
  //! @{

  /*! returns the `Is_on_x_identification_1` functor.
   */
  Is_on_x_identification_1 is_on_x_identification_1_object() const;

  //! @}
};

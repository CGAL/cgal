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
 * Models of the concept `AosClosedBoundaryTraits_1` handle curves theat are
 * closed at least on one side of the two sides of the parameter space.
 *
 * \cgalRefines{AosTraits_2}
 *
 * When at least one category is `Closed_side_tag`, the traits must also provide
 * the functor type `Parameter_space_in_x_1` (and its accessor), which
 * determines the position of a point in the parameter space. If the parameter
 * space is open on both sides, then the position of a point is trivial: it is
 * always in the interior of the parameter space.
 *
 * \cgalRefines{CopyConstructible,Assignable,DefaultConstructible}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{CGAL::Arrangement_on_curve_1::Circle_traits_1<Kernel>}
 * \cgalHasModelsEnd
 */
class AosClosedBoundaryTraits_1 {
public:
  //! \name Functor Types
  //! @{

  /*! Functor that determines the position of a point in the parameter space.;
   *! either `LEFT_BOUNDARY`, `INTERIOR`, or `RIGHT_BOUNDARY`.
   *! \pre either `Left_side_category` or `Right_side_category` are `Closed_side_tag`.
   */
  typedef unspecified_type Parameter_space_in_x_1;

  //! @}

  //! \name Accessing Functor Objects
  //! @{

  /*! returns the `Parameter_space_in_x_1` functor.
   *! \pre either `Left_side_category` or `Right_side_category` are `Closed_side_tag`.
   */
  Parameter_space_in_x_1 parameter_space_in_x_1_object() const;

  //! @}
};

// Copyright (c) 2026 Tel-Aviv University (Israel).
// All rights reserved.
// This file is part of CGAL (www.cgal.org).
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// Author(s): Efi Fogel         <efifogel@gmail.com>

namespace AocTraits {

/*! \ingroup PkgArrangementOnCurve1ConceptsFunctionObjects
 * \cgalConcept
 *
 * \cgalRefines{CopyConstructible}
 *
 * \cgalHasModelsBegin
 * \cgalHasModels{AocClosedBoundaryTraits_1::Parameter_space_in_x_1}
 * \cgalHasModelsEnd
 */
class ParameterSpaceInX_1 {
public:
  //! \name Operations
  //! A model of this concept must provide:
  //! @{

  /*! determines the location of the pre-image \f$t\f$ of the point `p` along
   * the interval parameter space.
   *
   * If \f$t\f$ lies on the left boundary, returns `CGAL::Arrangement_on_curve_1::LEFT_BOUNDARY`;
   * otherwise, if \f$t\f$ lies on the right boundary, returns `CGAL::Arrangement_on_curve_1::RIGHT_BOUNDARY`;
   * otherwise, returns `CGAL::Arrangement_on_curve_1::INTERIOR`.
   */
  Parameter_space operator()(const AocTraits::Point_1& p) const;

  //! @}
};

}

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
 * \cgalHasModels{AocTraits_1::Compare_x_1}
 * \cgalHasModelsEnd
 */
class CompareX_1 {
public:
  //! \name Operations
  //! A model of this concept must provide:
  //! @{

  /*! returns `CGAL::SMALLER`, `CGAL::EQUAL`, or `CGAL::LARGER` according to
   * the ordering of points `p1` and `p2` along the 1D carrier curve.
   *
   * \pre When the parameter space is identified, each of the points `p1` and
   *      `p2` must not coincide with the identification point.
   */
  Comparison_result operator()(const AocTraits::Point_1& p1, const AocTraits::Point_1& p2) const;

  //! @}
};

}

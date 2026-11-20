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
 * \cgalHasModels{AocIdentifiedBoundaryTraits_1::Is_on_x_identification_1}
 * \cgalHasModelsEnd
 */
class IsOnXIdentification_1 {
public:
  //! \name Operations
  //! A model of this concept must provide:
  //! @{

  /*! determines whether the point `p` coincides with the identification point.
   */
  bool operator()(AocTraits::Point_1& point) const;

  /*! determines whether the point `p` coincides with the identification curvepoint.
   */
  bool operator()(const AocTraits::Point_1& p) const;

  //! @}
};

}

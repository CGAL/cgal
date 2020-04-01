// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================


#ifndef CGAL_ACK_ENUMS_H
#define CGAL_ACK_ENUMS_H 1

namespace CGAL {

namespace internal {

  enum Three_valued {

    ROOT_OF_FIRST_SET = 1,
    ROOT_OF_BOTH_SETS = 0,
    ROOT_OF_SECOND_SET=-1

  };

} // namespace internal

/*!
 * \brief Represents different strategies how to handle
 * degenerate cases during the analysis
 *
 * Currently, there are two possible strategies implemented. See the
 * constructor of \c Curve_analysis_2 for more details.
 */
enum Degeneracy_strategy {

    SHEAR_STRATEGY = 0,
    EXCEPTION_STRATEGY = 1,
    SHEAR_ONLY_AT_IRRATIONAL_STRATEGY = 2
};

} //namespace CGAL

#endif

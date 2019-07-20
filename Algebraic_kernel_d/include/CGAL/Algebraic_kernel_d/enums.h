// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
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

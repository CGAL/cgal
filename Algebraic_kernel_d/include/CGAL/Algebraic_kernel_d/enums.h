// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://hemmer@scm.gforge.inria.fr/svn/cgal/trunk/Polynomial/include/CGAL/Polynomial.h $
// $Id: Polynomial.h 47254 2008-12-06 21:18:27Z afabri $
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

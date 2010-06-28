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
// Author(s)     : Michael Kerber    <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_KERNEL_D_2_H
#define CGAL_ALGEBRAIC_KERNEL_D_2_H

#include <CGAL/basic.h>

#include <CGAL/Algebraic_kernel_d_1.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_curve_kernel_2.h>

namespace CGAL {

template<typename Coefficient> class Algebraic_kernel_d_2
  : public CGAL::Algebraic_curve_kernel_2
               < CGAL::Algebraic_kernel_d_1 <Coefficient> >
{};
 
} //namespace CGAL



#endif // CGAL_ALGEBRAIC_KERNEL_D_1_H

// Copyright (c) 2006 Inria Lorraine (France). All rights reserved.
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
// $URL: $
// $Id: $
// 
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

#ifndef CGAL_GBRS_SOLVE_2_H
#define CGAL_GBRS_SOLVE_2_H

#include <CGAL/Gbrs_polynomial_2.h>

CGAL_BEGIN_NAMESPACE

// solve given the precision, returns de number of roots
int solve_2(const Rational_polynomial_2&,const Rational_polynomial_2&);

CGAL_END_NAMESPACE

#endif	// CGAL_GBRS_SOLVE_2_H

// Max-Planck-Institute Saarbruecken (Germany).
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
// $URL$
// $Id$
//
// Author(s)     : Michael Hemmer


#include <CGAL/Modular_arithmetic/Residue_type.h>

namespace CGAL{

int Residue::prime_int = 67111067;
double Residue::prime =67111067.0;
double Residue::prime_inv =1/67111067.0;
const double Residue::CST_CUT = std::ldexp( 3., 51 );

const CGALi::Modular_arithmetic_needs_ieee_double_precision 
Residue::modular_arithmetic_needs_ieee_double_precision 
= CGALi::Modular_arithmetic_needs_ieee_double_precision();
}

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

/*! \file CGAL/Modular.h
    \brief Defines the class CGAL::Modular and CGAL::Modular_traits.
 
    Provides the \c CGAL::Modular_traits specialization for the build in number 
    types. 
*/

#include <CGAL/Modular_arithmetic/Modular_type.h>

namespace CGAL{
    int Modular::prime_int = 67111067;
    double Modular::prime =67111067.0;
    double Modular::prime_inv =1/67111067.0;
    
    const double Modular::CST_CUT = std::ldexp( 3., 51 );

}

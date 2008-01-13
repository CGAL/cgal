// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
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
// $URL:$
// $Id:$
//
// Author(s)     : Michael Hemmer, Dominik Hülse
//
// ============================================================================

#ifndef NiX_EXTENDED_EUCLIDEAN_ALGORITHM_H
#define NiX_EXTENDED_EUCLIDEAN_ALGORITHM_H 1

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template< class NT > 
NT extended_euclidean_algorithm(const NT& a_, const NT& b_, NT& u, NT& v){

    typedef Algebraic_structure_traits<NT> AST;
    typename AST::Div_mod div_mod;
    typename AST::Unit_part unit_part;
    typename AST::Integral_division idiv;
    
    NT unit_part_a(unit_part(a_));
    NT unit_part_b(unit_part(b_));
    
    NT a(idiv(a_,unit_part_a));
    NT b(idiv(b_,unit_part_b));
    
    NT x(0),y(1),last_x(1),last_y(0);
    NT temp, quotient; 
  
    //TODO: unroll to avoid swapping 
    while (b != 0){
        temp = b;
        div_mod(a,b,quotient,b);
        a = temp;
     
        temp = x;
        x = last_x-quotient*x;
        last_x = temp;
      
        temp = y;
        y = last_y-quotient*y;
        last_y = temp;
    }
    u = last_x * unit_part_a;
    v = last_y * unit_part_b;

    CGAL_precondition(unit_part(a) == NT(1));
    CGAL_precondition(a == a_*u + b_*v);
    return a; 
}

CGAL_END_NAMESPACE

#endif // NiX_EXTENDED_EUCLIDEAN_ALGORITHM_H //

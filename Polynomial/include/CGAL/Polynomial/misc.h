// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
//
//
// Author(s)     : Michael Hemmer 

#ifndef CGAL_POLYNOMIAL_MISC_H
#define CGAL_POLYNOMIAL_MISC_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial/fwd.h>

namespace CGAL{
namespace internal{

// template meta function Innermost_coefficient_type
// returns the tpye of the innermost coefficient 
template <class T> struct Innermost_coefficient_type{ typedef T Type; };
template <class Coefficient_type> 
struct Innermost_coefficient_type<Polynomial<Coefficient_type> >{
    typedef typename Innermost_coefficient_type<Coefficient_type>::Type Type; 
};

// template meta function Dimension
// returns the number of variables 
template <class T> struct Dimension{ static const int value = 0;};
template <class Coefficient_type> 
struct Dimension<Polynomial<Coefficient_type> > {
    static const int value = Dimension<Coefficient_type>::value + 1 ; 
};

} // namespace internal
} // namespace CGAL

#endif // CGAL_POLYNOMIAL_MISC_H

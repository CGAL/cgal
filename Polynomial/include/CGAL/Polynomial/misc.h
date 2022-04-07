// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
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

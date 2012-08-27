// Copyright (c) 2005  Stanford University (USA).
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
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_INTERNAL_CONSTRUCT_FUNCTION_H
#define CGAL_POLYNOMIAL_INTERNAL_CONSTRUCT_FUNCTION_H
#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class Fn>
struct Construct_function
{

    typedef Fn result_type;
    typedef typename result_type::NT argument_type;

/*template <class It>
result_type operator()(It b,
           It e) const {
  return result_type(b,e);
  }*/

//! construct high degree polynomials
    result_type operator()(const argument_type &a0,
        const argument_type &a1=0) const
    {
        argument_type v[] = {a0, a1};
        return result_type(v, v+2);
    }

//! construct high degree polynomials
    result_type operator()(const argument_type &a0,
        const argument_type &a1,
        const argument_type &a2,
        const argument_type &a3=0) const
    {
        argument_type v[] = {a0, a1, a2, a3};
	int size=4;
	while (v[size-1]==argument_type(0)) --size;
        return result_type(v, v+size);
    }

//! construct high degree polynomials
    result_type operator()(const argument_type &a0,
        const argument_type &a1,
        const argument_type &a2,
        const argument_type &a3,
        const argument_type &a4,
        const argument_type &a5=0,
        const argument_type &a6=0,
        const argument_type &a7=0) const
    {
        argument_type v[] = {a0, a1, a2, a3, a4, a5, a6, a7};
	int size=8;
	while (v[size-1]==argument_type(0)) --size;
        return result_type(v, v+size);
    }

//! construct high degree polynomials
    result_type operator()(const argument_type &a0,
        const argument_type &a1,
        const argument_type &a2,
        const argument_type &a3,
        const argument_type &a4,
        const argument_type &a5,
        const argument_type &a6,
        const argument_type &a7,
        const argument_type &a8,
        const argument_type &a9=0,
        const argument_type &a10=0,
        const argument_type &a11=0,
        const argument_type &a12=0,
        const argument_type &a13=0,
        const argument_type &a14=0,
        const argument_type &a15=0,
        const argument_type &a16=0,
        const argument_type &a17=0,
        const argument_type &a18=0,
        const argument_type &a19=0) const
    {
        argument_type v[] = {
            a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10,
            a11, a12, a13, a14, a15, a16, a17, a18, a19
        };
	int size=20;
	while (v[size-1]==argument_type(0)) --size;
        return result_type(v, v+size);
    }
};

} } } //namespace CGAL::POLYNOMIAL::internal
#endif

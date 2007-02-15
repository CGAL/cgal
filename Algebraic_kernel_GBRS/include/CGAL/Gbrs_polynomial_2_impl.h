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
// $URL$
// $Id$
// 
//
// Author(s)     : Luis Peñaranda <penarand@loria.fr>

#ifndef CGAL_GBRS_POLYNOMIAL_2_IMPL_H
#define CGAL_GBRS_POLYNOMIAL_2_IMPL_H

CGAL_BEGIN_NAMESPACE

template<class T>
inline Rational_polynomial_2& Rational_polynomial_2::operator*=(const T &s){
	return(*this*=Gmpz(s));
};

template<class T>
Rational_polynomial_2 Rational_polynomial_2::operator*(const T &n)const{
	Rational_polynomial_2 r(*this);
	return (r*=n);
};

template<class T>
inline Rational_polynomial_2 operator*
(const T &n,const Rational_polynomial_2 &p){
	return (p*n);
};

CGAL_END_NAMESPACE

#endif	// CGAL_GBRS_POLYNOMIAL_2_IMPL_H

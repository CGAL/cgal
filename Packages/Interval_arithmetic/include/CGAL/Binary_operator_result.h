// Copyright (c) 2005  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_BINARY_OPERATOR_RESULT_H
#define CGAL_BINARY_OPERATOR_RESULT_H

#include <CGAL/basic.h>

// This class helps finding out the result type of mixed operators +-*/.
// For example it answers what the type of double+int is.

// This class is meant to be specialized for some number types pairs,
// when a mixed operator is defined.

CGAL_BEGIN_NAMESPACE

template < typename T1, typename T2 >
struct Binary_operator_result;

// T1 == T2
template < typename T >
struct Binary_operator_result <T, T>
{ typedef T type; };

// T1 == int
template < typename T >
struct Binary_operator_result <T, int>
{ typedef T type; };

// T2 == int
template < typename T >
struct Binary_operator_result <int, T>
{ typedef T type; };

CGAL_END_NAMESPACE

#endif // CGAL_BINARY_OPERATOR_RESULT_H

// Copyright (c) 2010  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_INTERNAL_GCC_CPP0X_H
#define CGAL_INTERNAL_GCC_CPP0X_H

// Enable C++0x features with GCC -std=c++0x (even when not specified at build time)
// See http://gcc.gnu.org/projects/cxx0x.html .

#if defined __GNUC__ && defined __GXX_EXPERIMENTAL_CXX0X__

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3) // GCC >= 4.3

// (tested with Fedora 10's g++ 4.3.2)

#undef CGAL_CFG_NO_CPP0X_ARRAY
#undef CGAL_CFG_NO_CPP0X_DECLTYPE
#undef CGAL_CFG_NO_CPP0X_DEFAULT_TEMPLATE_ARGUMENTS_FOR_FUNCTION_TEMPLATES
#undef CGAL_CFG_NO_CPP0X_ISFINITE
#undef CGAL_CFG_NO_CPP0X_LONG_LONG
#undef CGAL_CFG_NO_CPP0X_RVALUE_REFERENCE
#undef CGAL_CFG_NO_CPP0X_TUPLE
#undef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES

#endif // GCC >= 4.3

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4) // GCC >= 4.4

#undef CGAL_CFG_NO_CPP0X_AUTO
#undef CGAL_CFG_NO_CPP0X_DELETED_AND_DEFAULT_FUNCTIONS
#undef CGAL_CFG_NO_CPP0X_INITIALIZER_LISTS

#endif // GCC >= 4.4

#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 5) // GCC >= 4.5

#undef CGAL_CFG_NO_CPP0X_LAMBDAS

#endif // GCC >= 4.5

// Still not available in 4.5 :
// CGAL_CFG_NO_CPP0X_DELEGATING_CONSTRUCTORS

#endif // __GNUC__ && __GXX_EXPERIMENTAL_CXX0X__

#endif // CGAL_INTERNAL_GCC_CPP0X_H

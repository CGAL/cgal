// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Sylvain Pion

// CGAL_CFG_BOOL_IN_TEMPLATE_BUG.cpp
// ---------------------------------------------------------------------
// This program is used by cgal_configure.
// The following documentation will be pasted in the generated configfile.
// ---------------------------------------------------------------------

//| This flag is set, if the compiler ICEs on some combination of boolean
//| and (&&) within template boolean arguments. (e.g. g++ 3.3).

template < typename >
struct type_to_bool { static const bool value = false; };

template < bool >
struct bool_to_type { typedef void type; };

template < typename T >
struct P {
  template < typename U >
  P(U, typename bool_to_type< type_to_bool<U>::value && true >::type* = 0)
  {}
};

int main()
{
  P<int> p(0);
  return 0;
}

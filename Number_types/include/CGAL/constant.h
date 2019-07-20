// Copyright (c) 2006  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_CONSTANT_H
#define CGAL_CONSTANT_H

#include <CGAL/config.h>
#include <CGAL/tss.h>

namespace CGAL {

// The function constant<T, int i>() returns a const reference to T(i).
// TODO : is it worth documenting ?

template < typename T, int i >
inline
const T&
constant()
{
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(T, t,i);
  return t;
}

} //namespace CGAL

#endif // CGAL_CONSTANT_H

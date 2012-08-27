// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Nico Kruithof

//#define CGAL_TRIANGULATION_DONT_USE_SHORT_NAMES

#include <CGAL/Triangulation_3.h>

bool del = false;

#include <CGAL/_test_types.h>
#include <CGAL/_test_cls_triangulation_simplex_3.h>

// Explicit instantiation of the whole class :
template class CGAL::Triangulation_3<K>;

int main()
{
  typedef CGAL::Triangulation_3<K>                               Cls3;

  _test_cls_triangulation_simplex_3( Cls3() );

  return 0;
}

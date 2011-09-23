// Copyright (c) 2005, 2006  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
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
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_KERNEL_TYPE_MAPPER_H
#define CGAL_KERNEL_TYPE_MAPPER_H

#include <CGAL/basic.h>

namespace CGAL {

// This is a tool to obtain the K2::Point_2 from K1 and K1::Point_2.
// Similarly for other kernel types.

// TODO : add more specializations ?  Use a different mechanism ?

template < typename T, typename K1, typename K2 >
struct Type_mapper
{
  typedef T type; // By default, assume same type (e.g. Object).
};

// Then we specialize for all kernel objects.

#define CGAL_Kernel_obj(X) \
  template < typename K1, typename K2 > \
  struct Type_mapper < typename K1::X, K1, K2 > \
  { typedef typename K2::X type; };

#include <CGAL/Kernel/interface_macros.h>

} //namespace CGAL

#endif // CGAL_KERNEL_TYPE_MAPPER_H

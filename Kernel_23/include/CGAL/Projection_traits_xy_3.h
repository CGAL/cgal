// Copyright (c) 1997-2010  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mariette Yvinec

#ifndef CGAL_PROJECTION_TRAITS_XY_3_H
#define CGAL_PROJECTION_TRAITS_XY_3_H

#include <CGAL/internal/Projection_traits_3.h>

namespace CGAL { 

template < class R >
class Projection_traits_xy_3
  : public internal::Projection_traits_3<R,2>
{};
  
} //namespace CGAL 

#endif // CGAL_PROJECTION_TRAITS_XY_3_H

// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Sylvain Pion
//                 Andreas Fabri   <Andreas.Fabri@geometryfactory.com>

#ifndef CGAL_TRIANGULATION_UTILS_2_H
#define CGAL_TRIANGULATION_UTILS_2_H

#include <CGAL/license/TDS_2.h>


#include <CGAL/triangulation_assertions.h>

namespace CGAL { 
template < class T = void >
struct Triangulation_cw_ccw_static_2 {

static const int ccw_map[3];
static const int cw_map[3];
};
template < class T >
const int Triangulation_cw_ccw_static_2<T>::ccw_map[3] = {1, 2, 0};

template < class T >
const int Triangulation_cw_ccw_static_2<T>::cw_map[3] = {2, 0, 1};

class Triangulation_cw_ccw_2 
  : public  Triangulation_cw_ccw_static_2<>
{
public:
  static int ccw(const int i) 
    {
      CGAL_triangulation_precondition( i >= 0 && i < 3);
      return ccw_map[i];
    }

  static int cw(const int i)
    {
      CGAL_triangulation_precondition( i >= 0 && i < 3);
      return cw_map[i];
    }
};

} //namespace CGAL

#endif //CGAL_TRIANGULATION_UTILS_2_H

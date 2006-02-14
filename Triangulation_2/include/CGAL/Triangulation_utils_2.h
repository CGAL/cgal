// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//                 Sylvain Pion    <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri   <Andreas.Fabri@geometryfactory.com>
#ifndef CGAL_TRIANGULATION_UTILS_2_H
#define CGAL_TRIANGULATION_UTILS_2_H

#include <CGAL/triangulation_assertions.h>

CGAL_BEGIN_NAMESPACE 
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

CGAL_END_NAMESPACE

#endif //CGAL_TRIANGULATION_UTILS_2_H

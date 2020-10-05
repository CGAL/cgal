// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

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
// Author(s)     : Carine Bonetto
//                 Mariette Yvinec

#ifndef CGAL_TRIANGULATION_PS_STREAM_H
#define CGAL_TRIANGULATION_PS_STREAM_H

#include <CGAL/license/Triangulation_2.h>



#ifdef CGAL_TRIANGULATION_2_H
namespace CGAL {
template <class Gt,class Tds>
PS_Stream& operator << (PS_Stream& ps, const Triangulation_2<Gt,Tds> &t)

{
 t.draw_triangulation(ps);
  return ps;
}
} //namespace CGAL
#endif // CGAL_TRIANGULATION_2_H


#ifdef CGAL_DELAUNAY_TRIANGULATION_2_H
namespace CGAL {
template < class Gt, class Tds >
PS_Stream& operator << (PS_Stream& ps,
                        const Delaunay_triangulation_2<Gt,Tds> &t)
{
 t.draw_triangulation(ps);
 return ps;
}
} //namespace CGAL
#endif // CGAL_DELAUNAY_TRIANGULATION_2_H

#ifdef CGAL_CONSTRAINED_TRIANGULATION_2_H
namespace CGAL {
template < class Gt, class Tds>
PS_Stream& operator<<(PS_Stream& ps,
                      const Constrained_triangulation_2<Gt,Tds> &t)
{

 t.draw_triangulation(ps);
 return ps;
}
} //namespace CGAL
#endif // CGAL_CONSTRAINED_TRIANGULATION_2_H


#ifdef CGAL_REGULAR_TRIANGULATION_2_H
namespace CGAL {
template < class Gt, class Tds >
PS_Stream& operator << (PS_Stream& ps,
                        Regular_triangulation_2<Gt,Tds> &t)
{
  t.draw_triangulation(ps);
  return ps;
}
} //namespace CGAL
#endif // CGAL_REGULAR_TRIANGULATION_2_H

#endif //CGAL_TRIANGULATION_PS_STREAM





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
// $URL$
// $Id$
// 
//
// Author(s)     : Carine Bonetto
//                 Mariette Yvinec

#ifndef CGAL_TRIANGULATION_PS_STREAM_H
#define CGAL_TRIANGULATION_PS_STREAM_H


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



	

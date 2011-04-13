// ============================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/IO/Triangulation_ps_stream.h
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Carine Bonetto
//                 Mariette Yvinec
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_TRIANGULATION_PS_STREAM_H
#define CGAL_TRIANGULATION_PS_STREAM_H


#ifdef CGAL_TRIANGULATION_2_H
CGAL_BEGIN_NAMESPACE
template <class Gt,class Tds>
PS_Stream& operator << (PS_Stream& ps, const Triangulation_2<Gt,Tds> &t)

{
 t.draw_triangulation(ps);
  return ps;
}
CGAL_END_NAMESPACE
#endif // CGAL_TRIANGULATION_2_H


#ifdef CGAL_DELAUNAY_TRIANGULATION_2_H
CGAL_BEGIN_NAMESPACE
template < class Gt, class Tds >
PS_Stream& operator << (PS_Stream& ps, 
			const Delaunay_triangulation_2<Gt,Tds> &t)
{
 t.draw_triangulation(ps);
 return ps; 
}
CGAL_END_NAMESPACE
#endif // CGAL_DELAUNAY_TRIANGULATION_2_H

#ifdef CGAL_CONSTRAINED_TRIANGULATION_2_H
CGAL_BEGIN_NAMESPACE
template < class Gt, class Tds>
PS_Stream& operator<<(PS_Stream& ps,
		      const Constrained_triangulation_2<Gt,Tds> &t)
{

 t.draw_triangulation(ps);
 return ps;
}
CGAL_END_NAMESPACE
#endif // CGAL_CONSTRAINED_TRIANGULATION_2_H


#ifdef CGAL_REGULAR_TRIANGULATION_2_H
CGAL_BEGIN_NAMESPACE
template < class Gt, class Tds >
PS_Stream& operator << (PS_Stream& ps, 
			Regular_triangulation_2<Gt,Tds> &t)
{
  t.draw_triangulation(ps);
  return ps;
}
CGAL_END_NAMESPACE
#endif // CGAL_REGULAR_TRIANGULATION_2_H

#endif //CGAL_TRIANGULATION_PS_STREAM



	





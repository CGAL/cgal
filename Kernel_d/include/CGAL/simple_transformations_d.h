// Copyright (c) 2000,2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
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
// Author(s)     : Michael Seel

#ifndef CGAL_SIMPLE_TRANSFORMATIONS_D_H
#define CGAL_SIMPLE_TRANSFORMATIONS_D_H

#include <CGAL/Homogeneous_d.h>

namespace CGAL {
/*{\Moptions outfile=Simple_transformations.man}*/
/*{\Manpage{}{}{Simple Affine Transformations}{}}*/

template <class R>
Point_d<R> reflect(const Point_d<R>& p, const Point_d<R>& q)
/*{\Mfunc  returns |p| reflected at point |q|. }*/
{ return Point_d<R>(p + 2*(q-p)); }

template <class R>
Point_d<R> reflect(const Point_d<R>& p, const Line_d<R>& l)
/*{\Mfunc  returns |p| reflected across line |l|.}*/
{ // reflect point across line through r and s
  Point_d<R> c = l.projection(p);
  return Point_d<R>(p + 2*(c-p));
}

template <class R>
Sphere_d<R> reflect(const Sphere_d<R>& s, const Point_d<R>& p)
/*{\Mfunc returns |s| reflected across point |p|.}*/
{ std::vector< Point_d<R> > B(s.points_begin(),s.points_end());
  typename std::vector< Point_d<R> >::iterator it;
  for (it = B.begin(); it != B.end(); ++it) *it = reflect(*it,p);
  return Sphere_d<R>(p.dimension(),B.begin(),B.end());
}

template <class R>
Sphere_d<R> reflect(const Sphere_d<R>& s, const Line_d<R>& l)
/*{\Mfunc returns |s| reflected across line |l|.}*/
{ std::vector< Point_d<R> > B(s.points_begin(),s.points_end());
  typename std::vector< Point_d<R> >::iterator it;
  for (it = B.begin(); it != B.end(); ++it) *it = reflect(*it,l);
  return Sphere_d<R>(l.dimension(),B.begin(),B.end());
}

template <class R>
Segment_d<R> reflect(const Segment_d<R>& s, const Point_d<R>& p)
/*{\Mfunc returns |s| reflected across point |p|.}*/
{ return Segment_d<R>(reflect(s.point(0),p),reflect(s.point(1),p)); }

template <class R>
Segment_d<R> reflect(const Segment_d<R>& s, const Line_d<R>& l)
/*{\Mfunc returns |s| reflected across line |l|.}*/
{ return Segment_d<R>(reflect(s.point(0),l),reflect(s.point(1),l)); }

template <class R>
Ray_d<R> reflect(const Ray_d<R>& r, const Point_d<R>& p)
/*{\Mfunc returns |s| reflected across point |p|.}*/
{ return Ray_d<R>(reflect(r.point(0),p),reflect(r.point(1),p)); }

template <class R>
Ray_d<R> reflect(const Ray_d<R>& r, const Line_d<R>& l)
/*{\Mfunc returns |s| reflected across line |l|.}*/
{ return Ray_d<R>(reflect(r.point(0),l),reflect(r.point(1),l)); }

template <class R>
Line_d<R> reflect(const Line_d<R>& l, const Point_d<R>& p)
/*{\Mfunc returns |s| reflected across point |p|.}*/
{ return Line_d<R>(reflect(l.point(0),p),reflect(l.point(1),p)); }

template <class R>
Line_d<R> reflect(const Line_d<R>& s, const Line_d<R>& l)
/*{\Mfunc returns |s| reflected across line |l|.}*/
{ return Line_d<R>(reflect(s.point(0),l),reflect(s.point(1),l)); }

} //namespace CGAL
#endif // CGAL_SIMPLE_TRANSFORMATIONS_D_H

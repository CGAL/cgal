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
// release       : $CGAL_Revision: CGAL-0.9-I-06 $
// release_date  : $CGAL_Date: 1998/03/11 $
//
// file          : include/CGAL/Polygon_2.C
// source        :
// revision      : 1.8a
// revision_date : 13 Mar 1998
// author(s)     : Wieger Wesselink <wieger@cs.ruu.nl>
//
// coordinator   : Utrecht University
//
// ============================================================================

#ifndef CGAL_POLYGON_2_H
#include <CGAL/Polygon_2.h>
#endif // CGAL_POLYGON_2_H

//-----------------------------------------------------------------------//
//                          operator==
//-----------------------------------------------------------------------//

CGAL_BEGIN_NAMESPACE

template <class _Traits, class _Container1, class _Container2>
bool operator==( const Polygon_2<_Traits,_Container1> &x,
                 const Polygon_2<_Traits,_Container2> &y )
{
  CGAL_polygon_precondition( (x.size() != 0) || (y.size() != 0));

  if (x.size() != y.size()) return false;

  typename Polygon_2<_Traits,_Container1>::Vertex_const_iterator x_iter =
    x.vertices_begin();

  typename Polygon_2<_Traits,_Container2>::Vertex_const_iterator y_iter =
    std::find(y.vertices_begin(), y.vertices_end(), *x.vertices_begin());

  // if y doesn't contain the first point of x ...
  if (y_iter == y.vertices_end()) return false;

  ++x_iter;
  ++y_iter;

  while (y_iter != y.vertices_end()) {
    if (*x_iter != *y_iter) return false;
    ++x_iter;
    ++y_iter;
  }

  y_iter = y.vertices_begin();
  while (x_iter != x.vertices_end()) {
    if (*x_iter != *y_iter) return false;
    ++x_iter;
    ++y_iter;
  }

  return true;
}

//-----------------------------------------------------------------------//
//                          operator>>
//-----------------------------------------------------------------------//

template <class _Traits, class _Container>
istream &operator>>(istream &is, Polygon_2<_Traits,_Container>& p)
{
  int n; // number of vertices
  is >> n;

  typename _Traits::Point_2 point;

  for (int i=0; i<n; i++) {
    is >> point;
    p.push_back(point);
  }

  return is;
}

//-----------------------------------------------------------------------//
//                          operator<<
//-----------------------------------------------------------------------//

template <class _Traits, class _Container>
ostream &operator<<(ostream &os, const Polygon_2<_Traits,_Container>& p)
{
  typename Polygon_2<_Traits,_Container>::Vertex_const_iterator i;

  switch(os.iword(IO::mode)) {
    case IO::ASCII :
      os << p.size() << ' ';
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << *i << ' ';
      }
      return os;

    case IO::BINARY :
      os << p.size();
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << *i;
      }
      return os;

    default:
      os << "Polygon_2(" << endl;
      for (i = p.vertices_begin(); i != p.vertices_end(); ++i) {
        os << "  " << *i << endl;
      }
      os << ")" << endl;
      return os;
  }
}

CGAL_END_NAMESPACE

//-----------------------------------------------------------------------//
//                          transform
//-----------------------------------------------------------------------//

#ifdef CGAL_REP_CLASS_DEFINED
#ifndef CGAL_POLYGON_TRAITS_2_H
#include <CGAL/Polygon_traits_2.h>
#endif // CGAL_POLYGON_TRAITS_2_H

CGAL_BEGIN_NAMESPACE

template <class Transformation, class _Traits, class _Container>
Polygon_2<_Traits,_Container>
transform(const Transformation& t, const Polygon_2<_Traits,_Container>& p)
{
  typedef typename Polygon_2<_Traits,_Container>::Vertex_const_iterator VI;
  Polygon_2<_Traits,_Container> result;
  for (VI i = p.vertices_begin(); i != p.vertices_end(); ++i)
    result.push_back(t(*i));
  return result;
}

CGAL_END_NAMESPACE

#endif // CGAL_REP_CLASS_DEFINED


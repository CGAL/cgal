// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 01
//
// file          : include/CGAL/bops_dcel_out.h
// package       : bops (2.2)
// source        : include/CGAL/bops_dcel_out.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef DCEL_OUT_H
#define DCEL_OUT_H

#include <iostream>

CGAL_BEGIN_NAMESPACE

#ifdef CGAL__DCEL_DEFS_H
std::ostream& operator<<(std::ostream& o, const _Dcel_Color& c) {
  o << int(c);
  return o;
}
#endif

#ifdef CGAL__DCEL_VERTEX_H
template <class I>
std::ostream& operator<<(std::ostream& o, const _Dcel_vertex_type<I>* e) {
  if( e != NULL )
    o << 'v' << (*e).index();
  else
    o << "nil";
  return o;
}
#endif

#ifdef CGAL__DCEL_FACE_H
template <class I>
std::ostream& operator<<(std::ostream& o, const _Dcel_face_type<I>* e) {
  if( e != NULL )
    o << 'f' << (*e).index();
  else
    o << "nil";
  return o;
}
#endif

#ifdef CGAL__DCEL_EDGE_H
template <class I>
std::ostream& operator<<(std::ostream& o, const _Dcel_edge_type<I>* e) {
  if( e != NULL )
    o << '[' << (*e).index() << ']';
  else
    o << "nil";
  return o;
}
#endif


#ifdef CGAL__DCEL_VERTEX_H
template <class I>
std::ostream& operator<<(std::ostream& o, const _Dcel_vertex_type<I>& vertex) {
  o << "V("    << vertex.index()
    << ", d"   << vertex.degree()
    << ", c"   << vertex.color()
    << ", HV=" << vertex.header()
    << ", pt=" << vertex.point() << ')';
  return o;
}
#endif

#ifdef CGAL__DCEL_FACE_H
template <class I>
std::ostream& operator<<(std::ostream& o, const _Dcel_face_type<I>& face) {
  o << "F("    << face.index()
    << ", c"   << face.color()
    << ", HF=" << face.header() << ')';
  return o;
}
#endif

#ifdef CGAL__DCEL_EDGE_H
template <class I>
std::ostream& operator<<(std::ostream& o, const _Dcel_edge_type<I>& edge) {
#ifdef CGAL__DCEL_DEBUG_ON_ALL
  o << "{ c=" << edge.color() << ", ";
  if( edge.V1() == NULL )
    o << "nil";
  else
    o << (_Dcel_vertex_type<I>)*(edge.V1());

  o << ", ";
  if( edge.V2() == NULL )
    o << "nil";
  else
    o << (_Dcel_vertex_type<I>)*(edge.V2());

  o << ", ";
  if( edge.F1() == NULL )
    o << "nil";
  else
    o << (_Dcel_face_type<I>)*(edge.F1());
  o << ", ";
  if( edge.F2() == NULL )
    o << "nil";
  else
    o << (_Dcel_face_type<I>)*(edge.F2());

  o << ", " << edge.P1() << ", " << edge.P2();

  o << '}';
#else
  o << "{ " << edge.index() << ": c" << edge.color() << ", "
    << edge.V1() << ", " << edge.V2() << ", "
    << edge.F1() << ", " << edge.F2() << ", "
    << edge.P1() << ", " << edge.P2() << '}';
#endif
  return o;
}
#endif


template<class T1, class T2>
std::ostream& operator<<(std::ostream& o, std::pair<T1,T2> p){
  o << p.first << ' ' << p.second;
  return o;
}


template<class Iterator>
std::ostream& print( std::ostream& o, char* text, Iterator begin, Iterator end) {
  Iterator it;
  int i= 0;
  for( it= begin; it != end; it++ )
    o << text << '[' << i++ << "]= " << *it << endl;
  return o;
}

CGAL_END_NAMESPACE

#endif // DCEL_OUT_H

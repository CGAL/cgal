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
// file          : include/CGAL/bops_dcel_vertex.h
// package       : bops (2.2)
// source        : include/CGAL/bops_dcel_vertex.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL__DCEL_VERTEX_H
#define CGAL__DCEL_VERTEX_H

#include <list>
#include <vector>

#include <CGAL/bops_dcel_element.h>

CGAL_BEGIN_NAMESPACE

/*
  VERTEX in the DCEL:
  -------------------
  vertex_type:      point, header, color;
  vertex:           typedef const _Dcel_vertex_type* _Dcel_vertex;
  container:        vector<_Dcel_vertex_type>
  vertex_iterator:  typedef vector<_Dcel_vertex_type>::const_iterator
                    vertex_iterator;
  conversion:       vertex and vertex_iterator are type-identical

  *---------------------------------------------------------------------*
  *A Point< R > is stored as a Object to avoid templatizing of*
  *class _Dcel_vertex_type!                                        *
  *---------------------------------------------------------------------*

  CGAL-point:       Point< R > pt;
                    if( assign( pt, vertex->point()) )
                      this is a CGAL-point;
                    else
                      this is not a CGAL-point;

                    ...

  CGAL-object:      Object obj= make_object(pt);
  
*/

template <class I>
class _Dcel_vertex_type : public _Dcel_element_type {
public:
  typedef typename I::Point                                 Point;

#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_4
  typedef _Dcel_point_smaller_x<Point>        Point_smaller;
  typedef set<Point, Point_smaller>                Points_container;
  typedef typename Points_container::const_iterator const_points_iterator;
  typedef typename Points_container::iterator       points_iterator;
  //typedef Point*                   const_points_iterator;
  typedef _Dcel_edge_type<I>* edge;
#else
  //typedef typename I::points_iterator       points_iterator;
  typedef typename I::const_points_iterator const_points_iterator;
  typedef typename I::const_edges_iterator  edge;
#endif

  _Dcel_vertex_type() : _Dcel_element_type(), _degree(0) {}

  _Dcel_vertex_type(
     //const Point& pt,
     const_points_iterator pt,
     int ind,
     _Dcel_Color col = _NO_COLOR
  )
     : _Dcel_element_type(ind, col),
     _point(pt), _header(NULL), _degree(0) {}

  _Dcel_vertex_type(
     //const Point& pt,
     const_points_iterator pt,
     _Dcel_Color col = _NO_COLOR)
        : _Dcel_element_type(col),
          _point(pt), _header(NULL), _degree(0) {}

  Point point() const { return *_point; }
  //Point& point() { return (Point)*(points_iterator)_point; }
  //const Point& point() const { return _point; }
  //Point& point() { return _point; }

  int  degree() const { return _degree; }
  edge header() const { return _header; }

protected:
  edge&  header() { return _header; }
  void   header(edge h)       { _header= h; }

  friend class _Dcel_base<I>;

private:
  //Point _point;
  const_points_iterator _point;
  edge _header;
  int  _degree; /* the degree of this node */
};

CGAL_END_NAMESPACE

#endif /* CGAL__DCEL_VERTEX_H */


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
// file          : include/CGAL/bops_dcel_defs.h
// package       : bops (2.2)
// source        : include/CGAL/bops_dcel_defs.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL__DCEL_DEFS_H
#define CGAL__DCEL_DEFS_H

#include <CGAL/basic.h>

#include <utility>
#include <list>
#include <vector>
#include <set>
#include <algorithm>

#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Object.h>


#ifdef CGAL__DCEL_DEBUG_ON

# define  _DOUT cout

#  define CGAL__BOPS_DCEL_DEBUG(_t_) _DOUT << _t_ << flush
#  define CGAL__BOPS_DCEL_DEBUG_LN(_t_) _DOUT << _t_ << endl << flush
#  define CGAL__BOPS_DCEL_DEBUG_VAR(_t_,_v_) \
          _DOUT << _t_ << _v_ << endl << flush
#  define CGAL__BOPS_DCEL_DEBUG_PAIR(_t_,_p1_,_p2_)\
          _DOUT << _t_ << '(' << _p1_ << ',' << _p2_ << ')' << flush
#  define CGAL__BOPS_DCEL_DEBUG_ITERATOR(_text_, _begin_, _end_)\
          print(_DOUT, _text_, _begin_, _end_) << flush

#else

#  define CGAL__BOPS_DCEL_DEBUG(_t_)
#  define CGAL__BOPS_DCEL_DEBUG_LN(_t_) 
#  define CGAL__BOPS_DCEL_DEBUG_VAR(_t_,_v_)
#  define CGAL__BOPS_DCEL_DEBUG_PAIR(_t_,_p1_,_p2_)
#  define CGAL__BOPS_DCEL_DEBUG_ITERATOR(_text_, _begin_, _end_) 

#endif /* CGAL__DCEL_DEBUG_ON */

CGAL_BEGIN_NAMESPACE

/* Colors for the DCEL */
class _Dcel_Color {
public:
  _Dcel_Color() : _c(0) {}
  _Dcel_Color( const _Dcel_Color& c ) : _c(c._c) {}
  _Dcel_Color( int c ) : _c(c) {}

  void operator=(int i) { _c= i; }
  void operator=(_Dcel_Color other) { _c= other._c; }
  operator int() const {return _c;}
  _Dcel_Color operator|( const _Dcel_Color& c) const {
    return _Dcel_Color( _c | c._c );
  }
private:
  int _c;
};

const _Dcel_Color   _NO_COLOR= 0;
const _Dcel_Color   _UNCOLORED= 0;
const _Dcel_Color   _RED= 1;
const _Dcel_Color   _BLACK= 2;
const _Dcel_Color   _RED_AND_BLACK= 3;
const _Dcel_Color   _TWICE_COLORED= 3;



template<class I> class _Dcel_base;
template<class I> class _Dcel_vertex_type;
template<class I> class _Dcel_face_type;
template<class I> class _Dcel_edge_type;


/* _Dcel_point_compare: 
 * ------------------------
 *   predicate structure for the comparison of points in a container
 *   see also in _Dcel<I>::insert_new_vertex()
 */
template <class point_type, class vertex_type>
struct _Dcel_point_compare {
    _Dcel_point_compare() {}
    _Dcel_point_compare(const point_type& p) : _pt(p) {}
    bool operator()( const vertex_type& v ) {
      return _pt == v.point();
    }
    point_type _pt;
};

template <class point_type>
struct _Dcel_point_smaller_x :
  CGAL_STD::binary_function<point_type, point_type, bool>
{
  bool operator()(const point_type& p1, const point_type& p2) const {
    return p1.x() <= p2.x();
  }
};


template <class _I>
struct _Dcel_defs : public _I {
  typedef _Dcel_defs<_I>                I;

#ifdef CGAL_CFG_INCOMPLETE_TYPE_BUG_5
  _Dcel_edge_type<I>   _dummy_edge;
  _Dcel_face_type<I>   _dummy_face;
  _Dcel_vertex_type<I> _dummy_vertex;
#endif // ...

  typedef std::vector<_Dcel_edge_type<I> >   Edges_container;
  typedef std::vector<_Dcel_face_type<I> >   Faces_container;
  typedef std::vector<_Dcel_vertex_type<I> > Vertices_container;

  typedef typename Edges_container::const_iterator    const_edges_iterator;
  typedef typename Faces_container::const_iterator    const_faces_iterator;
  typedef typename Vertices_container::const_iterator const_vertices_iterator;

  typedef typename Edges_container::iterator          edges_iterator;
  typedef typename Faces_container::iterator          faces_iterator;
  typedef typename Vertices_container::iterator       vertices_iterator;

#ifndef CGAL_CFG_INCOMPLETE_TYPE_BUG_4
  typedef typename _I::Point                          Point;
#endif

  typedef _Dcel_point_smaller_x<Point>  Point_smaller;
  typedef std::set<Point, Point_smaller>          Points_container;
  typedef typename Points_container::const_iterator   const_points_iterator;
  typedef typename Points_container::iterator         points_iterator;


/* old typedefs (they should not be used further!!!) */
  typedef typename Edges_container::const_iterator edge_iterator;
  typedef typename Faces_container::const_iterator face_iterator;
  typedef typename Vertices_container::const_iterator vertex_iterator;
};

CGAL_END_NAMESPACE

#endif /* CGAL__DCEL_DEFS_H */

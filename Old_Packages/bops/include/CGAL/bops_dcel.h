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
// file          : include/CGAL/bops_dcel.h
// package       : bops (2.2)
// source        : include/CGAL/bops_dcel.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL__DCEL_H
#define CGAL__DCEL_H

#include <CGAL/bops_dcel_defs.h>
#include <CGAL/bops_dcel_base.h>
#include <CGAL/min_sqr_distance.h>


#ifdef CGAL_CFG_NO_SCOPE_MEMBER_FUNCTION_PARAMETERS
#define CGAL_CFG_RETURN_TYPE_BUG_2  /* to avoid troubles with sun-pro */
#endif

CGAL_BEGIN_NAMESPACE

template<class NT>
inline bool _number_type_is_non_exact(const NT&) { return false; }

inline bool _number_type_is_non_exact(const double&) { return true; }
inline bool _number_type_is_non_exact(const float&) { return true; }


/*
  Implementation of the Double Connected Edge List (DCEL)
  -------------------------------------------------------

  template <class I>
  class _Dcel : public _Dcel_base<I>;

  template <class I>
  struct Bops_dcel : public _Dcel< _Dcel_defs<I> >;

*/



template <class I>
class _Dcel : public _Dcel_base<I>
{
public:
  typedef _Dcel_base<I>  dcel_base;
  typedef typename I::R R;
  typedef typename I::Point Point;
  typedef std::pair<int,int> epair;

  _Dcel() { use_epsilon= false; }
  _Dcel(const _Dcel<I>& dl) { *this= dl; }
  _Dcel(const std::list<epair>& eds, const std::list<Point>& pts ) {
    use_epsilon= false;
    insert(eds,pts);
  }
  _Dcel(const std::list<epair>& eds, const std::vector<Point>& pts ) {
    use_epsilon= false;
    insert(eds,pts);
  }

# ifdef  CGAL_CFG_RETURN_TYPE_BUG_2
  bool colorize(const std::list<Point>& pgon, const _Dcel_Color& col) {
    __point_list= &pgon;
    return colorize(col);
  }
# else
  bool colorize(const std::list<Point>& pgon, const _Dcel_Color& col);
# endif  



  Point point(const_vertices_iterator v) const {
    return (*v).point();
  }

  const_vertices_iterator find(const Point& pt) const {
     //_Dcel_point_compare<I> pred(pt);
     //_Dcel_point_compare<Point, _Dcel_vertex_type<I> > pred(pt);
     //return  find_if(_v_list.begin(), _v_list.end(), pred );
     for(const_vertices_iterator it= _v_list.begin(); it != _v_list.end(); it++)
       if( compare_points(pt, (*it).point()) ) return it;
     return _v_list.end();
  }

  
  const_vertices_iterator null_vertex(void) const {
     return _v_list.end();
  }

private:
  double pts_epsilon;
  bool   use_epsilon;

  bool compare_points(const Point& p1, const Point& p2) const {
    if( use_epsilon ) {
      return abs(to_double(p1.x()-p2.x())) < pts_epsilon &&
             abs(to_double(p1.y()-p2.y())) < pts_epsilon;
    }
    else
      return p1 == p2;
  }

  void insert_points(const std::list<Point>& points) {
    if( (use_epsilon= _number_type_is_non_exact(R::FT(0))) ) {
      /* calculates epsilon for non-exact number types */
      I traits;
      pts_epsilon= minimal_square_distance(
		    points.begin(), points.end(), traits);
      pts_epsilon= sqrt(pts_epsilon)/3.0;
    }

    //std::vector<vertex_iterator> c_it; /* help-array */
    /* insert vertices (also in help array) */
    _v_list.reserve(points.size());
    int n= 0;
    std::list<typename I::Point>::const_iterator it;
    for( it= points.begin(); it != points.end(); it++ )          
        c_it.push_back( insert_point( *it, n++ ) );
  
    CGAL__BOPS_DCEL_DEBUG_ITERATOR("c_it", c_it.begin(), c_it.end() );
#   ifdef CGAL__BOPS_DEBUG_ON
    print(cout, "c_it", c_it.begin(), c_it.end() );
    print(cout, "V", vertex_begin(), vertex_end() );
    print(cout, "P", _point_list.begin(), _point_list.end() ); // from dcel-base
    copy(_point_list.begin(), _point_list.end(),
         ostream_iterator<Point>(cout, "\n" ) );
#   endif
    return;
  }

  void insert_points(const std::vector<Point>& points)
  {
    typedef typename R::FT FT;
    if( (use_epsilon= _number_type_is_non_exact(FT(0))) ) {
      /* calculates epsilon for non-exact number types */
      I traits;
      pts_epsilon= minimal_square_distance(
		    points.begin(), points.end(), traits);
      pts_epsilon= std::sqrt(pts_epsilon)/3.0;
    }
    // std::vector<vertex_iterator> c_it; /* help-array */
    /* insert vertices (also in help array) */
    _v_list.reserve(points.size());
    int n= 0;
    std::vector< typename I::Point >::const_iterator it;
    for( it= points.begin(); it != points.end(); it++ ) 
        c_it.push_back( insert_point( *it, n++ ) );

    CGAL__BOPS_DCEL_DEBUG_ITERATOR("c_it", c_it.begin(), c_it.end() );
#   ifdef CGAL__BOPS_DEBUG_ON
    print(cout, "c_it", c_it.begin(), c_it.end() );
    print(cout, "V", vertex_begin(), vertex_end() );
    print(cout, "P", _point_list.begin(), _point_list.end() );
#   endif
    return;
  }

  std::vector<const_vertices_iterator> c_it; /* help-array */

public:
  void insert(const std::list< epair >& eds, const std::list< Point >& pts ) {
      insert_points(pts);
#     ifdef  CGAL_CFG_RETURN_TYPE_BUG_2
      // because of compiler (g++, Solaris)  problems,
      // we have to do this:
        __edges= &eds;
        insert_edges();
#     else
        insert_edges(eds);
#     endif  
      c_it= std::vector<const_vertices_iterator>();
  }

  void insert(const std::list< epair >& eds, const std::vector< Point >& pts ) {
      insert_points(pts);
#     ifdef  CGAL_CFG_RETURN_TYPE_BUG_2
      // because of compiler (g++, Solaris)  problems,
      // we have to do this
        __edges= &eds;
        insert_edges();
#     else
        insert_edges(eds);
#     endif  
      c_it= std::vector<const_vertices_iterator>();
  }

private:

  const_vertices_iterator insert_point( const Point& pt, int index ) {
    return insert_new_vertex(pt, index);
  }
  
  const_vertices_iterator insert_new_vertex(
    const Point& pt,
    int index,
    _Dcel_Color col = _NO_COLOR
  )
  {
     const_vertices_iterator it= find(pt);

     if( it == null_vertex() ) {
       /* vertex does not exist, hence append vertex */
       const_points_iterator pt_it= dcel_base::insert(pt);
       //it= dcel_base::insert( _Dcel_vertex_type<I>(pt, index, col) );
       it= dcel_base::insert( _Dcel_vertex_type<I>(pt_it, index, col) );
     }

     return it;
  }

private:
#ifdef  CGAL_CFG_RETURN_TYPE_BUG_2
  bool colorize(const _Dcel_Color& col);
  void insert_edges();
  const std::list<Point>* __point_list;
  const std::list<epair>* __edges;
#else
  void insert_edges(const std::list<epair>& eds);
#endif

};

template <class I>
struct Bops_dcel : public _Dcel< _Dcel_defs<I> > {
  typedef _Dcel< _Dcel_defs<I> >  dcel;
  typedef typename dcel::const_edges_iterator    edge_iterator;
  typedef typename dcel::const_faces_iterator    face_iterator;
  typedef typename dcel::const_vertices_iterator vertex_iterator;
};

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/bops_dcel.C>
#endif 

#endif /* CGAL__DCEL_H */


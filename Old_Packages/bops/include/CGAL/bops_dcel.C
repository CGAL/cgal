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
// file          : include/CGAL/bops_dcel.C
// package       : bops (2.2)
// source        : include/CGAL/bops_dcel.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL__DCEL_C
#define CGAL__DCEL_C

#ifndef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/bops_dcel.h>
#endif

#include <CGAL/bops_V2E_rep.h>
#include <CGAL/bops_dcel_out.h>
#include <CGAL/bops_V2E_rep_out.h>

//#define CGAL__DCEL_DEBUG_ON

CGAL_BEGIN_NAMESPACE

template<class I>
struct _Dcel_V2E_compare : public I {
   typedef const _Dcel_vertex_type<I>* vertex;
   typedef typename I::Point Point;
   typedef typename I::Direction Direction;
   Point p0;

   _Dcel_V2E_compare(vertex v0) : p0((*v0).point()) { }
   bool operator()(vertex v1, vertex v2 ) const { return compare(v1,v2); }
   bool compare(vertex v1, vertex v2 ) const {
     Point p1= (*v1).point(), p2= (*v2).point();
     Direction d1= I::direction_of(p0,p1);
     Direction d2= I::direction_of(p0,p2);
     return d1 < d2; 
   }
};




/*
  Implementation of the Double Connected Edge List (DCEL)
  -------------------------------------------------------

  template <class I> class DCEL__Dcel;
*/


#ifdef  CGAL_CFG_RETURN_TYPE_BUG_2 
template <class I> void _Dcel<I> :: insert_edges() {
  const std::list<epair>& edges= *__edges;
#else 
template <class I>
void _Dcel<I> :: insert_edges( const std::list<epair>& edges) {
#endif // CGAL_CFG_RETURN_TYPE_BUG_2

  typedef _V2E_rep_type<const_vertices_iterator,edge_iterator,
          _Dcel_V2E_compare<I> > V2E_rep_dcel;
  typedef _V2E_rep_base_type<const_vertices_iterator,edge_iterator>
          V2E_rep_base;


    /* insert edges */
    _e_list.reserve(edges.size());
    std::list< std::pair<int,int> >::const_iterator ee;
    for( ee= edges.begin(); ee != edges.end(); ee++ )
      _Dcel_base<I>::insert(
	_Dcel_edge_type<I>(c_it[(*ee).first], c_it[(*ee).second])
      );
    
    CGAL__BOPS_DCEL_DEBUG_LN("VERTICES inserted: ");
    CGAL__BOPS_DCEL_DEBUG_ITERATOR("dcel", this->begin(), this->end() );
    CGAL__BOPS_DCEL_DEBUG_ITERATOR("PT",
	  _point_list.begin(), _point_list.end() );

    /* template< class T_vertex, class T_edge, class T_compare > */
    V2E_rep_dcel v2e_rep( c_it.size(), edges.size() );

    /* insertion in vertex-to-edge rep. */
    /* const_vertices_iterator v0= _v_list.begin(); */
    const_vertices_iterator v1, v2;
    int iv1, iv2;
    for( edge_iterator e= _e_list.begin(); e != _e_list.end(); e++ ) {
      v1= (*e).V1();
      v2= (*e).V2();
      iv1= (*v1).index();
      iv2= (*v2).index();
      CGAL__BOPS_DCEL_DEBUG_VAR("insert: ", e);
      CGAL__BOPS_DCEL_DEBUG_PAIR("", iv1, iv2 );
      CGAL__BOPS_DCEL_DEBUG_LN("");
      v2e_rep.insert( iv1, (*e).V1(), iv2, (*e).V2(), e);
    }

    CGAL__BOPS_DCEL_DEBUG((V2E_rep_base)v2e_rep);

    v2e_rep.sort_vertices_CCW();
    CGAL__BOPS_DCEL_DEBUG((V2E_rep_base)v2e_rep);

    /* construction of vertex cycles */
    construct_vertex_cycles( v2e_rep );
    CGAL__BOPS_DCEL_DEBUG_LN("POINTERS inserted: ");
    CGAL__BOPS_DCEL_DEBUG_ITERATOR("dcel", this->begin(), this->end() );
    CGAL__BOPS_DCEL_DEBUG_ITERATOR("V", vertex_begin(), vertex_end() );
  
    /* construction of face cycles */
    construct_face_cycles();
    CGAL__BOPS_DCEL_DEBUG_LN("FACES inserted: ");
    CGAL__BOPS_DCEL_DEBUG_ITERATOR("dcel", this->begin(), this->end() );
    CGAL__BOPS_DCEL_DEBUG_ITERATOR("F", face_begin(), face_end() );

    return;
  }


template <class I>
#ifdef CGAL_CFG_RETURN_TYPE_BUG_2
bool _Dcel<I>::colorize(const _Dcel_Color& col) {
  typedef typename I::Point Point;
  const std::list<Point>& pgon= *__point_list;
#else 
bool _Dcel<I>::colorize(
			     const std::list<typename I::Point>& pgon,
			     const _Dcel_Color& col)
{
  typedef typename I::Point Point;
  typedef typename I::vertices_iterator vertices_iterator;
  typedef typename I::edges_iterator  edges_iterator;
  typedef typename I::faces_iterator  faces_iterator;
  typedef typename I::edge_iterator   edge_iterator;

#endif // CGAL_CFG_RETURN_TYPE_BUG_2
  
/* This routine colorizes the interior of a simple polygon 'pgon'
   given by a list of its vertices in counterclockwise order(!).
                                      ^^^^^^^^^^^^^^^^^^^^^^^^^
   The color is given by the input variable 'col' and will be set using
   the or operation, i.e. vertex::set_color( vertex::color() | col );

   return: true  iff everything ok
           false if e.g. polygon not valid
*/


  if( pgon.size() < 3 ) return false;

  std::list<Point>::const_iterator it_p;

  /* get last point */
  std::list<Point>::const_reverse_iterator rit_p= pgon.rbegin();
  Point pt= (*rit_p);
  Point v1_pt;

  const_vertices_iterator v0, v1;
  v0= find( pt );       // find first vertex
  vertices_iterator v= (vertices_iterator)v0;
  (*v).set_color( (*v0).color() | col);
  edge_iterator a, e0;
  edges_iterator e;
  faces_iterator f;

  std::list<const_vertices_iterator> vlist_deg2;
  // all vertices with degree greater than 2
  
  for(it_p= pgon.begin(); it_p != pgon.end(); it_p++) { /* for all vertices */
    v0= v;
    if( (*v0).degree() > 2 ) 
      vlist_deg2.push_back(v0);

    /* find second vertex 'v1' and edge [v0,v1] */
    a= e0= begin(v0);
    v1= (*a).opposite_vertex(v0);
    //if( (*v1).point() != *it_p )
    if( !compare_points((*v1).point(), *it_p ) )
    {
      for( a= next(e0,v0); a != e0; a= next(a,v0) ) {
        v1= (*a).opposite_vertex(v0);
        //if( (*v1).point() == *it_p ) break; /* found !! */
        if( compare_points((*v1).point(), *it_p) ) break; /* found !! */
      }
      if( a==e0) return false; /* ABORT: --> polygon not valid */
    }

    e= (edges_iterator)a;
    (*e).set_color( (*e).color()|col);
    v= (vertices_iterator)v1;
    (*v).set_color( (*v).color()|col);
    face_iterator face0= (*e).left_face(v0); /* v0 -> v1 ** pgon is CCW */
    face_iterator face1= (*e).opposite_face(face0);
    f= (faces_iterator)face0;
    (*f).set_color( (*f).color()|col);
    if( face_inside_in_face(face1, face0) ) {
      f= (faces_iterator)face1;
      (*f).set_color( (*f).color()|col);
    }
  }

  
  std::list<const_vertices_iterator>::const_iterator vit;
  face_iterator f1, f2;
  for( vit= vlist_deg2.begin(); vit != vlist_deg2.end(); vit++) {
    v0= *vit;
    a= e0= begin(v0);
    do {
      if( !((*a).color() & col ) ) {
        f1= (*a).F1();
        f2= (*a).F2();
        if( (*f1).color()&col && (*f2).color()&col ) {
          e= (edges_iterator)a;
          (*e).set_color( (*a).color()|col );
          v1= (*a).opposite_vertex(v0);
          v= (vertices_iterator)v1;
          (*v).set_color( (*v).color()|col);
          vlist_deg2.push_back(v1);
        }
      }
      a= next(a,v0);
    }
    while( a != e0 );
  }

  CGAL__BOPS_DCEL_DEBUG_LN("colorized: ");
  CGAL__BOPS_DCEL_DEBUG_ITERATOR("E", this->begin(), this->end() );
  CGAL__BOPS_DCEL_DEBUG_ITERATOR("V", vertex_begin(), vertex_end() );
  CGAL__BOPS_DCEL_DEBUG_ITERATOR("F", face_begin(), face_end() );

  return true;
}

CGAL_END_NAMESPACE

#endif /* CGAL__DCEL_C */

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
// file          : include/CGAL/bops_dcel_base.C
// package       : bops (2.2)
// source        : include/CGAL/bops_dcel_base.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Wolfgang Freiseisen <Wolfgang.Freiseisen@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ======================================================================

#ifndef CGAL__DCEL_BASE_C
#define CGAL__DCEL_BASE_C

#ifndef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/bops_dcel_base.h>
#endif

#include <CGAL/bops_V2E_rep.h>

CGAL_BEGIN_NAMESPACE

//#defineCGAL_CFG_RETURN_TYPE_BUG_2
#ifdef CGAL_CFG_RETURN_TYPE_BUG_2
template <class I>
void _Dcel_base<I>::construct_vertex_cycles() {
  //typedef typename vector<_Dcel_vertex_type<I> >::const_iterator _vertex;
  //typedef typename vector<_Dcel_edge_type<I> >  ::const_iterator _edge;
  typedef const_vertices_iterator _vertex;
  typedef const_edges_iterator    _edge;
  typedef _V2E_rep_base_type<_vertex,_edge> V2E_rep_base_dcel;
  V2E_rep_base_dcel& v2e= (*_v2e);
#else
template <class I>
void _Dcel_base<I>::construct_vertex_cycles(
  typename _Dcel_base<I>::V2E_rep_base_dcel& v2e
)
{
#endif

    /* construction of vertex cycles */
    /* preconditions:
                vertex list _v_list is initialized: 
                edge   list _e_list is initialized: V1, V2
                vertex-to-edge representation is initialized and sorted
    */


  typedef typename I::const_edges_iterator    edge_iterator;
  typedef typename I::const_vertices_iterator const_vertices_iterator;

  typedef typename I::faces_iterator       faces_iterator;
  typedef typename I::vertices_iterator    vertices_iterator;
  typedef typename I::edges_iterator       edges_iterator;



    typename V2E_rep_base_dcel::header_iterator v1;
    typename V2E_rep_base_dcel::vertex_iterator v2, v2_begin;

    vertices_iterator vtex;
    edges_iterator e;
    _Dcel_edge_id  u;


    /* for all vertices v1 */
    for( v1= v2e.header_begin(); v1 != v2e.header_end(); v1++) {
      vtex= (vertices_iterator)(*(*v1).vertex()).vertex();
      (*vtex)._degree= (*v1).size();

      /* for all vertices v2 in "v1-list" */
      v2_begin= v2e.vertex_begin(v1);
      (*vtex)._header= (edges_iterator)(*v2_begin).edge();
      for( v2= v2_begin; v2 != v2e.vertex_end(v2); v2= v2e.vertex_next(v2) )
      { /* edge [v1,v2] */
        e= (edges_iterator)(*v2).edge();
        u= (*e).V1()==vtex ? EDGE_1 : EDGE_2;
        (*e).P(u)= (*v2e.vertex_next(v2)).edge(); 
      }
      (*e).P(u)= (*v2_begin).edge(); /* closing the circle */
      
    }
    return;
  }


template <class I>
void _Dcel_base<I>::construct_face_cycles() {
    /* construction of face cycles */
    /* preconditions:
          face list _f_list is not initialized, it will be cleared and filled
          edge list _e_list is initialized: V1, V2, P1, P2
    */
  typedef typename I::const_vertices_iterator const_vertices_iterator;
  typedef typename I::faces_iterator       faces_iterator;
  typedef typename I::vertices_iterator    vertices_iterator;
  typedef typename I::edges_iterator       edges_iterator;

    int m= number_of_edges();
    edges_iterator a= _e_list.begin();

    /* create face-list */
    int f= 2 + m - number_of_vertices();
    _f_list.erase(_f_list.begin(), _f_list.end());
    _f_list.reserve(f+1);
    for(int i= 0; i < f; i++)
      _f_list.push_back( _Dcel_face_type<I>(i) ); 
    
    faces_iterator s= _f_list.begin();; 
    int k= 1; /* edge-counter */
  
    _Dcel_edge_id   u;

    while ( k<= 2*m ) {
      if( (*a).F1() != NULL && (*a).F2() != NULL )
        /* if F1(a) and F2(a) are set, then take next edge a */
        a++;
      else { /* at least one of F1(a) or F2(a) is not initialized */
        /* which face is not initialized ? -> u */
        u= (*a).F1() == NULL ? FACE_1 : FACE_2;

        /* initialize face and its header */
        k++;
        (*a).F(u)= s;
        (*s).header()= a;

        /* get the corresponding vertex */
        const_vertices_iterator c= (*a).V(u);
      
        /* initialize all unset faces by face s and vertex c */
        edges_iterator e= a;
        for( e= (edges_iterator)(*e).P(u);
	     e != a;
	     e= (edges_iterator)(*e).P(u) )
	{
          u= (*e).V1() == c  ? VERTEX_2 : VERTEX_1;
          /* initialize face Fu(e) */
          k++;
          (*e).F(u)= s;
          c= (*e).V(u); /* set (new) vertex */
        }
        s++;
      }
    } /* while */

    return;
  }

CGAL_END_NAMESPACE
  
#endif /* CGAL__DCEL_C */

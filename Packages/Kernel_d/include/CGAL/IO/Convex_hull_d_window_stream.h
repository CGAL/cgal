// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Convex_hull_d_window_stream.h
// package       : Kernel_d
// chapter       : Basic
//
// source        : ddgeo/Convex_hull_d.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Susan Hert <hert@mpi-sb.mpg.de>
//
// implementation: Higher dimensional geometry
// ============================================================================
#ifndef CGAL_CONVEX_HULL_D_WINDOW_STREAM_H
#define CGAL_CONVEX_HULL_D_WINDOW_STREAM_H

#include <CGAL/basic.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/IO/Regular_complex_d_window_stream.h>


CGAL_BEGIN_NAMESPACE

/*{\Mtext \headerline{Low Dimensional Output Routines}
include |<CGAL/IO/Convex_hull_d_window_stream.h>|
\setopdims{2cm}{3cm}}*/

template <class R>
void d2_show(const Convex_hull_d<R>& C, CGAL::Window_stream& W);
/*{\Mfunc draws the convex hull |C| in window |W|.\\
\precond |dim == 2|. }*/

template <class R>
void d3_surface_map(const Convex_hull_d<R>& C, 
                    GRAPH< typename Convex_hull_d<R>::Point_d ,int>& G);
/*{\Mfunc constructs the representation of the surface of |\Mvar| as a 
bidirected LEDA graph |G|.\\ \precond |dim == 3|.}*/


template <class R>
void d2_show(const Convex_hull_d<R>& C, CGAL::Window_stream& W)
{ /* We first draw every simplex*/
  typedef typename Convex_hull_d<R>::Simplex_const_handle Simplex_handle;
  typedef typename Convex_hull_d<R>::Vertex_const_handle Vertex_handle;
  Simplex_handle S;
  forall_rc_simplices(S,C) {
    for (int v = ( C.is_unbounded_simplex(S)  ? 1 : 0); 
         v <= C.current_dimension(); v++) {
      // for each vertex except the anti - origin 
 
      for (int e = v + 1; e <= C.current_dimension(); e++) {
        // draw undrawn edges incident to vertex 
        if ( C.is_unbounded_simplex(S) ) W.set_line_width(3); 
        else W.set_line_width(1); 
        W.draw_segment(to_leda_point(C.point_of_simplex(S,v)),
                       to_leda_point(C.point_of_simplex(S,e))); 
      }
    }
  }
  /* Now we draw every point */
  typename Convex_hull_d<R>::Point_const_iterator pit;
  for (pit = C.points_begin(); pit != C.points_end(); ++pit) {
    W.draw_point(to_leda_point(*pit)); 
  }
}


template <class R> 
void d3_surface_map(const Convex_hull_d<R>& C, 
                    GRAPH< typename Convex_hull_d<R>::Point_d ,int>& G)
{ 
  typedef typename Convex_hull_d<R>::Vertex_const_handle  Vertex_handle;
  typedef typename Convex_hull_d<R>::Simplex_const_handle Simplex_handle;
  typedef typename Convex_hull_d<R>::Facet_const_handle   Facet_handle;
  typedef typename Convex_hull_d<R>::Point_d Point_d;
  typedef typename R::RT RT;

  G.clear();
  if (C.dimension() != 3) 
    CGAL_assertion_msg(0,"d3_surface_map: dim must be 3");
  if (C.current_dimension() < 3) {
    Hash_map<Vertex_handle,leda_node> node_for(nil);
    Vertex_handle v; Simplex_handle s;
    forall_rc_vertices(v,C) {
      node_for[v] = G.new_node(C.associated_point(v));
    }
    if (C.current_dimension() <= 0) { return; }
    if (C.current_dimension() == 1) {
      forall_rc_simplices(s,C) {
        if (C.is_bounded_simplex(s)) {
          leda_node v0 = node_for[C.vertex(s,0)];
          leda_node v1 = node_for[C.vertex(s,1)];
          leda_edge e01 = G.new_edge(v0,v1);
          leda_edge e10 = G.new_edge(v1,v0);
          G.set_reversal(e01,e10);
        }
      }
      return;
    }

    if (C.current_dimension() == 2) {
      leda_node_array<bool> untreated(G,true);

      typename R::Orthogonal_vector_d ortho_vector =
        C.kernel().orthogonal_vector_d_object();
      typename R::Point_d pc = C.center() + 
        ortho_vector(C.base_facet_plane(C.origin_simplex()));
      forall_rc_simplices(s,C) {
        if (C.is_bounded_simplex(s)) {
          for (int i = 0; i <= C.current_dimension(); i++) {
            leda_node vi = node_for[C.vertex(s,i)];
            if ( untreated[vi] ) {
              untreated[vi] = false;
              int j = (i + 1) % (C.current_dimension() + 1);  
              // a vertex different from i;
              int k = (i + 2) % (C.current_dimension() + 1); 
              leda_node vj = node_for[C.vertex(s,j)];
              leda_node vk = node_for[C.vertex(s,k)];
              std::vector< Point_d > V(4);
              V[0]=G[vi]; V[1]=G[vj]; V[2]=G[vk]; V[3]=pc;
              typename R::Orientation_d orientation =
                C.kernel().orientation_d_object();
              if ( orientation(V.begin(),V.end()) == POSITIVE ) {
                std::swap(vj,vk);
                std::swap(j,k); 
              }

              leda_edge efirst = G.new_edge(vi,vk); 
              // first edge incident to vi
              Simplex_handle scur = s; 
              int jcur = j, kcur = k, icur = i;

              while ( C.is_bounded_simplex(C.opposite_simplex(scur,jcur)) && 
                      C.opposite_simplex(scur,jcur) != s ) {
                // we have not reached the end nor closed the cycle
                kcur = C.index_of_opposite_vertex(scur,jcur);
                scur = C.opposite_simplex(scur,jcur);
                for (icur = 0; icur <= 2; icur++)
                  if (node_for[C.vertex(scur,icur)] == vi) break;
                jcur = 3 - icur - kcur;
                vk = node_for[C.vertex(scur,kcur)];
                G.new_edge(vi,vk);
              } 

              if (C.is_unbounded_simplex(C.opposite_simplex(scur,jcur))) {
                /* we also need to walk in the other direction */
 
                efirst = G.new_edge(efirst,vj,0,LEDA::before);  // 0 is etype
                scur = s; jcur = j; kcur = k; icur = i;  
                // restore initial situation
              
                while ( C.is_bounded_simplex(
                        C.opposite_simplex(scur,kcur)) ) {
                  // we have not reached the end 
                  jcur = C.index_of_opposite_vertex(scur,kcur);
                  scur = C.opposite_simplex(scur,kcur);
                  for (icur = 0; icur <= 2; icur++)
                    if (node_for[C.vertex(scur,icur)] == vi) break;
                  kcur = 3 - jcur -icur;
                  vj = node_for[C.vertex(scur,jcur)];
                  efirst = G.new_edge(efirst,vj,0,LEDA::before); //as above
                } //end while
              }// end if
            }// end if untreated
          }// end for i
        }// end if bounded
      }// end forall
    if (!G.make_map())
      error_handler(1,"chull::surface_graphrep: not bidirected"); 
    return;
    }
  }

  Facet_handle f; 
  Vertex_handle v;
  Hash_map<Vertex_handle,leda_node> node_for(nil);
  int facet_num = 0;

  std::list<Facet_handle> Surface = C.all_facets();
  typename std::list<Facet_handle>::iterator it;
  for(it = Surface.begin(); it != Surface.end(); ++it) {
    f = *it;
    ++facet_num;
    for (int i=0; i < C.current_dimension(); i++) {
      v = C.vertex_of_facet(f,i);
      if (!node_for[v])
        node_for[v] = G.new_node(C.associated_point(v));
    }
  }
  if ( 2*G.number_of_nodes() != facet_num + 4)
    error_handler(1,"d3_surface_map: node equation wrong.");
 
  leda_node_array<bool> untreated(G,true);
  for(it = Surface.begin(); it != Surface.end(); ++it) {
    f = *it;
    for (int i = 0; i < C.current_dimension(); i++) {
      leda_node vi = node_for[C.vertex_of_facet(f,i)];
      if ( untreated[vi] ) {
        untreated[vi] = false;
        int j = (i + 1) % C.current_dimension();
        // a vertex different from i;
        int k = (i + 2) % C.current_dimension();
        leda_node vj = node_for[C.vertex_of_facet(f,j)];
        leda_node vk = node_for[C.vertex_of_facet(f,k)];
        typename R::Orientation_d orientation_ =
          C.kernel().orientation_d_object();
        std::vector< Point_d > V(4);
        V[0]=G[vi]; V[1]=G[vj]; V[2]=G[vk]; V[3]=C.center();
        if ( orientation_(V.begin(),V.end()) == POSITIVE ) {
          std::swap(vk,vj); std::swap(k,j); 
        }

        G.new_edge(vi,vk);  // first edge incident to vi
        Facet_handle fcur = f; 
        int jcur = j; int kcur = k; int icur = i;

        while ( C.opposite_facet(fcur,jcur) != f ) {
          // we have not reached the end
          kcur = C.index_of_vertex_in_opposite_facet(fcur,jcur);
          fcur = C.opposite_facet(fcur,jcur);
          for (icur = 0; icur < 3; icur++)
             if ( node_for[C.vertex_of_facet(fcur,icur)] == vi ) break;
          jcur = 3 - icur - kcur;
          vk = node_for[C.vertex_of_facet(fcur,kcur)];
          G.new_edge(vi,vk);
        } 

      } // end if untreated
    } // end for i
  } // end forall
  if (G.number_of_edges() != (3*facet_num))
    error_handler(1,"d3_surface_map: wrong number of edges");
  if (!G.make_map())
    error_handler(1,"d3_surface_map: not bidirected"); 
}




CGAL_END_NAMESPACE
#endif //CGAL_CONVEX_HULL_D_WINDOW_STREAM_H



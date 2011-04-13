#ifndef CGAL_REGULAR_COMPLEX_D_WINDOW_STREAM_H
#define CGAL_REGULAR_COMPLEX_D_WINDOW_STREAM_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#ifndef _MSC_VER
#include <CGAL/Regular_complex_d.h>
#else
#include <CGAL/Regular_complex_d_MSC.h>
#endif
#include <CGAL/IO/Window_stream.h>
#include <LEDA/graph.h>
#include <LEDA/node_map.h>
#include <LEDA/node_map2.h>
#include <LEDA/d3_point.h>


CGAL_BEGIN_NAMESPACE

/*{\Mtext \headerline{Visualization and conversion in low dimensions}
The corresponding operations can be found in |Regular_complex_d_window_stream.h|}*/

template <class R>
void d2_show(const Regular_complex_d<R>& RC, Window_stream& W);
/*{\Mfunc displays the regular complex R in window W.
\precond |dim == 2|.}*/

template <class R>
void d2_map(const Regular_complex_d<R>& RC, GRAPH<leda_point,int>& G);
/*{\Mfunc constructs the representation of |R| as a bidirected graph |G|.
\precond |dim == 2|.}*/

template <class R>
void d3_graph(const Regular_complex_d<R>& RC, GRAPH<leda_d3_point,int>& G);
/*{\Mfunc constructs the representation of |R| as a bidirected graph |G|.
\precond |dim == 3|.}*/


template <class Pnt_2> 
leda_point to_leda_point(const Pnt_2& p) 
{ return leda_point(CGAL::to_double(p.x()), CGAL::to_double(p.y())); }

template <class R>
void d2_show(const Regular_complex_d<R>& RC, CGAL::Window_stream& W)
{ CGAL_assertion_msg(RC.dimension() == 2,"show_rc: dimension not 2.");
  W.set_line_width(1); 

  typename Regular_complex_d<R>::Simplex_const_iterator s;
  forall_rc_simplices(s,RC) {
    for (int v = 0; v <= RC.current_dimension(); v++) {
      W << to_leda_point(RC.associated_point(s,v));
      for (int e = v + 1; e <= RC.current_dimension(); e++) {
        leda_segment seg(to_leda_point(RC.associated_point(s,v)),
                         to_leda_point(RC.associated_point(s,e)));
        W << seg;
      }
    }
  } 
}


template <class R>
void d2_map(const Regular_complex_d<R>& RC, GRAPH<leda_point,int>& G)
{ 
  typedef Regular_complex_d<R>::Simplex_const_iterator Simplex_iterator;
  typedef Regular_complex_d<R>::Vertex_const_iterator Vertex_iterator;
  typedef Regular_complex_d<R>::Vertex_const_handle Vertex_handle;
  typedef Regular_complex_d<R>::Simplex_const_handle Simplex_handle;
  Vertex_iterator v;
  Simplex_iterator s;

  if (RC.dimension() != 2) CGAL_assertion_msg(0,"d2_map: dim must be 2.");
  G.clear();
  Unique_hash_map<Vertex_handle,leda_node> node_for(nil);
  
  forall_rc_vertices(v,RC) {
    node_for[v] = G.new_node(to_leda_point(RC.associated_point(v)));
  }
 
  if (RC.current_dimension() <= 0) return;
  if (RC.current_dimension() == 1) {
    forall_rc_simplices(s,RC) {
      leda_node v0 = node_for[RC.vertex(s,0)];
      leda_node v1 = node_for[RC.vertex(s,1)];
      leda_edge e01 = G.new_edge(v0,v1);
         // every dart is a clockwise boundary dart
      leda_edge e10 = G.new_edge(v1,v0);
      G.set_reversal(e01,e10);
    }
    return;
  }
  
  int T_num(0); // number of triangles 
  int B_num(0); // number of boundary edges

  forall_rc_simplices(s,RC) {
    ++T_num;
    for (int i = 0; i <= RC.current_dimension(); i++)
      if (RC.opposite_simplex(s,i) == Simplex_handle()) ++B_num;
  }
  int N = G.number_of_nodes();
  if (((B_num + T_num) % 2 != 0) || (N != 1 + (B_num +T_num)/2)) 
    CGAL_assertion_msg(0,"d2_map: wrong number of vertices");
     
  leda_node_array<bool> untreated(G,true);

  forall_rc_simplices(s,RC) {
    for (int i = 0; i <= RC.current_dimension(); i++) {
      leda_node vi = node_for[RC.vertex(s,i)];
      if ( untreated[vi] ) {
        untreated[vi] = false;
        int j = (i + 1) % (RC.current_dimension() + 1);  // a vertex different from i;
        int k = (i + 2) % (RC.current_dimension() + 1); 
        leda_node vj = node_for[RC.vertex(s,j)];
        leda_node vk = node_for[RC.vertex(s,k)];
        if ( orientation(G[vi],G[vj],G[vk])<0) {
          leda_swap(vk,vj); leda_swap(j,k); 
        }

        leda_edge efirst = G.new_edge(vi,vk);  // first edge incident to vi
        Simplex_handle scur = s; 
        int jcur = j; int kcur = k; int icur = i;

        while ( RC.opposite_simplex(scur,jcur) != Simplex_handle() && 
                RC.opposite_simplex(scur,jcur) != s ) {
          // we have not reached the end nor closed the cycle
          kcur = RC.index_of_opposite_vertex(scur,jcur);
          scur = RC.opposite_simplex(scur,jcur);
          for (icur = 0; icur <= 2; icur++)
             if (node_for[RC.vertex(scur,icur)] == vi) break;
          jcur = 3 - icur - kcur;
          vk = node_for[RC.vertex(scur,kcur)];
          G.new_edge(vi,vk);
        } 

        if (RC.opposite_simplex(scur,jcur) == Simplex_handle()) {
          /* we also need to walk in the other direction */
 
          efirst = G.new_edge(efirst,vj,0,LEDA::before);  // 0 is etype
          scur = s; 
          jcur = j; kcur = k; icur = i;  // restore initial situation

          while ( RC.opposite_simplex(scur,kcur) != Simplex_handle() ) {
            // we have not reached the end 
            jcur = RC.index_of_opposite_vertex(scur,kcur);
            scur = RC.opposite_simplex(scur,kcur);
            for (icur = 0; icur <= 2; icur++)
               if (node_for[RC.vertex(scur,icur)] == vi) break;
            kcur = 3 - jcur -icur;
            vj = node_for[RC.vertex(scur,jcur)];
            efirst = G.new_edge(efirst,vj,0,LEDA::before); //as above
          } //end while
        }// end if
      }// end if untreated
    }// end for i
  }// end forall
  if (G.number_of_edges() != (3*T_num + B_num))
    CGAL_assertion_msg(0,"Regular_complex_d::d2_map: wrong number of edges");
  if (!G.make_map())
    CGAL_assertion_msg(0,"Regular_complex_d::d2_map:not bidirected"); 
}


template <class R>
void d3_graph(const Regular_complex_d<R>& RC, 
              GRAPH< typename Regular_complex_d<R>::Point_d ,int>& G)
{ 
  typedef Regular_complex_d<R>::Simplex_const_iterator Simplex_iterator;
  typedef Regular_complex_d<R>::Vertex_const_iterator Vertex_iterator;
  typedef Regular_complex_d<R>::Vertex_const_handle Vertex_handle;

  Simplex_iterator s;
  Vertex_iterator v;

  CGAL_assertion_msg(RC.dimension() == 3,"d3_graph: dim must be 3.");
  G.clear();
  node_map2<bool> connected(G);
  Unique_hash_map<Vertex_handle,leda_node> node_for(nil);

  forall_rc_vertices(v,RC) {
    node_for[v] = G.new_node(RC.associated_point(v));
  }
 
  forall_rc_simplices(s,RC) {
    for (int i = 0; i <= RC.current_dimension(); i++)
      for (int j = 0; j <= RC.current_dimension(); j++) {
        Vertex_handle vert1 = RC.vertex(s,i);
        Vertex_handle vert2 = RC.vertex(s,j);
        if (vert1 != Vertex_handle() && vert2 != Vertex_handle()) {
          leda_node v1 = node_for[vert1];
          leda_node v2 = node_for[vert2];
          if (v1 == nil || v2 == nil)
            CGAL_assertion_msg(0,"why is this shitty node not initialized?");
          if (!connected(v1,v2)) {
            connected(v1,v2) = connected(v2,v1) = true;
            leda_edge e1 = G.new_edge(v1,v2);
            leda_edge e2 = G.new_edge(v2,v1);
            G.set_reversal(e1,e2);
          }
        }
      }
  }
}



CGAL_END_NAMESPACE
#endif //CGAL_REGULAR_COMPLEX_D_WINDOW_STREAM_H


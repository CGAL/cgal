// ======================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 21
//
// file          : 
// package       : Convex_hull_3 (2.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// source        : LEP dd_geo_kernel
// revision      : 2.1
// revision_date : 
// author(s)     : Kurt Mehlhorn
//                 Michael Seel
//
// coordinator   : MPI, Saarbruecken <Stefan.Schirra@mpi-sb.mpg.de>
// ======================================================================
/*******************************************************************************
+
+  LEP dd_geokernel 2.1
+
+  This file is part of the research version of a LEDA extension package,
+  that can be used free of charge in academic research and teaching. 
+  Any commercial use of this software requires a commercial license,
+  which is distributed by the Algorithmic Solutions GmbH, 
+  Postfach 151101, 66041 Saarbruecken, FRG (fax +49 681 31104).
+
+  Copyright (c) 1997-1998  by  Max-Planck-Institut fuer Informatik
+  Im Stadtwald, 66123 Saarbruecken, Germany     
+  All rights reserved.
+ 
*******************************************************************************/


template <class TRAITS,class POINT>
void regl_complex<TRAITS,POINT>::
check_topology() const
{ 
  rc_simplex s,t; 
  rc_vertex v;
  int i,j,k; 

  if (dcur == -1) {
    if (!all_verts.empty() || !all_simps.empty() ) 
      error_handler(1,
      "check_topology: dcur is -1 but there are vertices or simplices");
  }

  forall(v,all_verts) {
    if ( v != vertex(simplex(v),index(v)) )
      error_handler(1,
      "check_topology: vertex-simplex relationship corrupted");
  }

  forall(s,all_simps) {
    for(i = 0; i <= dcur; i++) {
      for (j = i + 1; j <= dcur; j++) {
        if (vertex(s,i) == vertex(s,j))
          error_handler(1,
          "check_topology: a simplex with two equal vertices"); 
      }
    }
  }

  forall(s,all_simps) {
    for(i = 0; i <= dcur; i++) {
      if ((t = opposite_simplex(s,i)) != nil) { 
        int l = index_of_opposite_vertex(s,i); 
        if (s != opposite_simplex(t,l) || 
            i != index_of_opposite_vertex(t,l))
          error_handler(1,
          "check_topology: neighbor relation is not symmetric"); 

        for (j = 0; j <= dcur; j++) {
          if (j != i) {
            // j must also occur as a vertex of t
            for (k = 0; k <= dcur && 
                   ( vertex(s,j) != vertex(t,k) || k == l); k++); 
            // forloop has no body
            if (k > dcur) 
              error_handler(1,
              "check_topology: too few shared vertices."); 
          }
        }
      }
    }
  }
}

template <class TRAITS,class POINT>
void regl_complex<TRAITS,POINT>::
check_topology_and_geometry() const
{ 
  check_topology();
  rc_vertex v;
  forall(v,all_verts) {
    if ( v == nil || identical(associated_point(v),regl_complex::nil_point) )
      error_handler(1,"check_topology_and_geometry: \
      vertex with nil_point or no associated point.");
  }

  rc_simplex s;
  forall(s,all_simps) {
    array<POINT> A(dcur + 1);
    for (int i = 0; i <= dcur; i++) 
      A[i] = associated_point(s,i);
    if ( !TRAITS::affinely_independent(A) )
      error_handler(1,"check_topology_and_geometry: \
      corners of some simplex are not affinely independent");
  }
}


template <class TRAITS,class POINT>
void d2_show(const regl_complex<TRAITS,POINT>& R, window& W)
{ /* We first draw every simplex*/
  rc_Simplex<TRAITS,POINT>* s;

  if (R.dim() != 2) 
    error_handler(1,"show_rc: dimension not 2.");
  W.set_line_width(1); 

#ifdef DDGEO_STL_ITERATORS
  regl_complex<TRAITS,POINT>::rc_simplex_iterator sit;
  for(sit =  R.simplices_begin(); 
      sit != R.simplices_end();
      ++sit) {
    s = *sit;
#else
  forall(s,R.all_simplices()) {
#endif
    for (int v = 0; v <= R.dcurrent(); v++) {
      W.draw_point(TRAITS::to_d2_point((R.associated_point(s,v))));
      for (int e = v + 1; e <= R.dcurrent(); e++) {
        W.draw_segment(TRAITS::to_d2_point((R.associated_point(s,v))),
                       TRAITS::to_d2_point((R.associated_point(s,e)))); 
      }
    }
  } 
}


template <class TRAITS,class POINT>
void d2_map(const regl_complex<TRAITS,POINT>& R, GRAPH<POINT,int>& G)
{ 
  typedef rc_Simplex<TRAITS,POINT>* rc_simplex;
  typedef rc_Vertex<TRAITS,POINT>*  rc_vertex;
  rc_vertex v;
  rc_simplex s;

  if (R.dim() != 2) 
    error_handler(1,"d2_map: dim must be 2.");
  G.clear();
  map<rc_vertex,node> node_for(nil);
  
#ifdef DDGEO_STL_ITERATORS
  regl_complex<TRAITS,POINT>::rc_vertex_iterator  vit;
  regl_complex<TRAITS,POINT>::rc_simplex_iterator sit;
  for(vit = R.vertices_begin(); 
      vit != R.vertices_end();
      ++vit) {
    v = *vit;
#else
  forall(v,R.all_vertices()) {
#endif
    node_for[v] = G.new_node(R.associated_point(v));
  }
 
  if (R.dcurrent() <= 0) return;
  
  if (R.dcurrent() == 1) {
#ifdef DDGEO_STL_ITERATORS
    for(sit =  R.simplices_begin(); 
        sit != R.simplices_end();
        ++sit) {
      s = *sit;
#else
    forall(s,R.all_simplices()) {
#endif
      node v0 = node_for[R.vertex(s,0)];
      node v1 = node_for[R.vertex(s,1)];
      edge e01 = G.new_edge(v0,v1);
         // every dart is a clockwise boundary dart
      edge e10 = G.new_edge(v1,v0);
      G.set_reversal(e01,e10);
    }
    return;
  }
  
  int T(0); // number of triangles 
  int B(0); // number of boundary edges
#ifdef DDGEO_STL_ITERATORS
  for(sit =  R.simplices_begin(); 
      sit != R.simplices_end();
      ++sit) {
    s = *sit;
#else
  forall(s,R.all_simplices()) {
#endif
    ++T;
    for (int i = 0; i <= R.dcurrent(); i++)
      if (R.opposite_simplex(s,i) == nil) ++B;
  }
  int N = G.number_of_nodes();
  if (((B + T) % 2 != 0) || (N != 1 + (B +T)/2)) 
    error_handler(1,"d2_map: wrong number of vertices");
     
  node_array<bool> untreated(G,true);
#ifdef DDGEO_STL_ITERATORS
  for(sit =  R.simplices_begin(); 
      sit != R.simplices_end();
      ++sit) {
    s = *sit;
#else
  forall(s,R.all_simplices()) {
#endif
    for (int i = 0; i <= R.dcurrent(); i++) {
      node vi = node_for[R.vertex(s,i)];
      if ( untreated[vi] ) {
        untreated[vi] = false;
        int j = (i + 1) % (R.dcurrent() + 1);  // a vertex different from i;
        int k = (i + 2) % (R.dcurrent() + 1); 
        node vj = node_for[R.vertex(s,j)];
        node vk = node_for[R.vertex(s,k)];
        if (TRAITS::orientation(G[vi],G[vj],G[vk])<0) {
          leda_swap(vk,vj);
          leda_swap(j,k); 
        }

        edge efirst = G.new_edge(vi,vk);  // first edge incident to vi
        rc_simplex scur = s; 
        int jcur = j; int kcur = k; int icur = i;

        while ( R.opposite_simplex(scur,jcur) && 
                R.opposite_simplex(scur,jcur) != s ) {
          // we have not reached the end nor closed the cycle
          kcur = R.index_of_opposite_vertex(scur,jcur);
          scur = R.opposite_simplex(scur,jcur);
          for (icur = 0; icur <= 2; icur++)
             if (node_for[R.vertex(scur,icur)] == vi) break;
          jcur = 3 - icur - kcur;
          vk = node_for[R.vertex(scur,kcur)];
          G.new_edge(vi,vk);
        } 

        if (R.opposite_simplex(scur,jcur) == nil) {
          /* we also need to walk in the other direction */
 
          efirst = G.new_edge(efirst,vj,0,LEDA_PREFIXLI before);  // 0 is etype
          scur = s; 
          jcur = j; kcur = k; icur = i;  // restore initial situation

          while ( R.opposite_simplex(scur,kcur) ) {
            // we have not reached the end 
            jcur = R.index_of_opposite_vertex(scur,kcur);
            scur = R.opposite_simplex(scur,kcur);
            for (icur = 0; icur <= 2; icur++)
               if (node_for[R.vertex(scur,icur)] == vi) break;
            kcur = 3 - jcur -icur;
            vj = node_for[R.vertex(scur,jcur)];
            efirst = G.new_edge(efirst,vj,0,LEDA_PREFIXLI before); //as above
          } //end while
        }// end if
      }// end if untreated
    }// end for i
  }// end forall
  if (G.number_of_edges() != (3*T + B))
    error_handler(1,"regl_complex::d2_map: wrong number of edges");
  if (!G.make_map())
    error_handler(1,"regl_complex::d2_map:not bidirected"); 
}

template <class TRAITS,class POINT>
void d3_graph(const regl_complex<TRAITS,POINT>& R, GRAPH<POINT,int>& G)
{ 
  typedef rc_Simplex<TRAITS,POINT>* rc_simplex;
  typedef rc_Vertex<TRAITS,POINT>*  rc_vertex;
  rc_simplex s;
  rc_vertex v;

  if (R.dim() != 3) 
    error_handler(1,"d3_graph: dim must be 3.");
  G.clear();
  node_map2<bool> connected(G);
  map<rc_vertex,node> node_for(nil);

#ifdef DDGEO_STL_ITERATORS
  regl_complex<TRAITS,POINT>::rc_vertex_iterator  vit;
  for(vit = R.vertices_begin(); 
      vit != R.vertices_end();
      ++vit) {
    v = *vit;
#else
  forall(v,R.all_vertices()) {
#endif
    node_for[v] = G.new_node(R.associated_point(v));
  }
 
#ifdef DDGEO_STL_ITERATORS
  regl_complex<TRAITS,POINT>::rc_simplex_iterator sit;
  for(sit =  R.simplices_begin(); 
      sit != R.simplices_end();
      ++sit) {
    s = *sit;
#else
  forall(s,R.all_simplices()) {
#endif
    for (int i = 0; i <= R.dcurrent(); i++)
      for (int j = 0; j <= R.dcurrent(); j++) {
        rc_vertex vert1 = R.vertex(s,i);
        rc_vertex vert2 = R.vertex(s,j);
        if (vert1 != nil && vert2 != nil) {
          node v1 = node_for[vert1];
          node v2 = node_for[vert2];
          if (v1 == nil || v2 == nil)
            error_handler(1,"why is this shitty node not initialized?");
          if (!connected(v1,v2)) {
            connected(v1,v2) = connected(v2,v1) = true;
            edge e1 = G.new_edge(v1,v2);
            edge e2 = G.new_edge(v2,v1);
            G.set_reversal(e1,e2);
          }
        }
      }
  }
}


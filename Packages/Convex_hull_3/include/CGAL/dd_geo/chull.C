// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 2000, August 21
// 
// source        : chull.fw
// file          : include/CGAL/dd_geo/chull.C
// package       : Convex_hull_3 (2.6)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 2.6
// revision_date : 21 Aug 2000 
// author(s)     : Kurt Mehlhorn
//                 Michael Seel
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 



template <class CHTRAITS,class POINT,class PLANE>
chull<CHTRAITS,POINT,PLANE>::
chull(int d) :
  regl_complex<CHTRAITS,POINT>(d)
{
  origin_simplex = nil;
  start_facet = nil;
  anti_origin = nil;
  number_of_points = number_of_vertices = 0;
  number_of_unbounded_simplices = number_of_bounded_simplices = 0;
  number_of_visibility_tests = 0;
  quasi_center = CHTRAITS::zero_vector(d);
}


template <class CHTRAITS,class POINT,class PLANE>
bool chull<CHTRAITS,POINT,PLANE>::
contains_in_base_facet(CHSIMPLEX s, const POINT& x) const
{
  array<POINT> A(dcur);
  for (int i = 1; i <= dcur; i++)
    A[i - 1] = point_of_simplex(s,i);
  return CHTRAITS::contained_in_simplex(A,x);
}

template <class CHTRAITS,class POINT,class PLANE>
void chull<CHTRAITS,POINT,PLANE>::
compute_equation_of_base_facet(CHSIMPLEX S)
{

  array<POINT> P(dcur);
  for (int i = 0; i < dcur; i++)
    P[i] = point_of_simplex(S,i + 1);
  S->base_facet = CHTRAITS::hyperplane_through_points(P,center(), - 1);

  #ifdef SELFCHECK
    { /* Let us check */
      for ( int i = 1; i <= dcur; i++)
        if (CHTRAITS::side(S->base_facet,point_of_simplex(S,i)) != 0)
          error_handler(1," does not support base ");
      if (CHTRAITS::side(S->base_facet,center()) != - 1)
        error_handler(1," quasi center on wrong side  ");
    }
  #endif


}


template <class CHTRAITS,class POINT,class PLANE>
rc_Vertex<CHTRAITS,POINT>* chull<CHTRAITS,POINT,PLANE>::
insert(const POINT& x)
{
  ch_vertex z = nil;
  all_pnts.append(x);
  number_of_points++;
  if (dcur == -1) { // |x| is the first point to be inserted

    ch_simplex outer_simplex; // a pointer to the outer simplex
    dcur = 0;  // we jump from dimension - 1 to dimension 0
    origin_simplex = new_simplex();  number_of_bounded_simplices ++;
    outer_simplex  = new_simplex();  number_of_unbounded_simplices ++;
    start_facet = origin_simplex;
    z = new_vertex(x);               number_of_vertices ++;
    assoc_vertex_with_simplex(origin_simplex,0,z);
      // z is the only point and the peak
    assoc_vertex_with_simplex(outer_simplex,0,anti_origin);
    set_neighbor(origin_simplex,0,outer_simplex,0);
    quasi_center = CHTRAITS::to_vector(x);


  }
  else if ( is_dimension_jump(x) ) {

    dcur++;
    z = new_vertex(x); number_of_vertices++;
    quasi_center = quasi_center + CHTRAITS::to_vector(x);
    dimension_jump(origin_simplex, z);
    clear_visited_marks(origin_simplex);
    list_item it;
    forall_items(it,all_simps) {
      ch_simplex S = (ch_simplex) all_simps[it];
      compute_equation_of_base_facet(S);
    }
    number_of_unbounded_simplices += number_of_bounded_simplices;
    if (dcur > 1) {
      start_facet = opposite_simplex(origin_simplex,dcur);
      ASSERT(!vertex_of_simplex(start_facet,0),start_facet_update_dj);
    }


  }
  else {

    if ( dcur == 0 ) {
      z = vertex_of_simplex(origin_simplex,0);
      assoc_point_with_vertex(z,x);
      return z;
    }
    list<ch_simplex> visible_simplices;
    int location = -1;
    ch_facet f = nil;


    int number_of_visited_simplices = 0;

    visibility_search(origin_simplex, x, visible_simplices,
                      number_of_visited_simplices, location, f);

    number_of_visibility_tests += number_of_visited_simplices;

    #ifdef COUNTS
      cout << "\nthe number of visited simplices in this iteration is ";
      cout << number_of_visited_simplices << endl;
    #endif

    clear_visited_marks(origin_simplex);


    #ifdef COUNTS
      cout << "\nthe number of bounded simplices constructed ";
      cout << " in this iteration is  " << visible_simplices.length() << endl;
    #endif

    number_of_bounded_simplices += visible_simplices.length();

    switch (location) {
      case -1:
        return nil;
      case 0:
        { for (int i = 0; i < dcur; i++) {
            if ( x == point_of_facet(f,i) ) {
              z = vertex_of_facet(f,i);
              assoc_point_with_vertex(z,x);
              return z;
            }
          }
          return nil;
        }
      case 1:
        { number_of_vertices++;
          z = new_vertex(x);
          list<ch_simplex> NewSimplices; // list of new simplices
          ch_simplex S;
          forall(S, visible_simplices) {
            assoc_vertex_with_simplex(S,0,z);
            for (int k = 1; k <= dcur; k++) {
              if (!is_base_facet_visible(opposite_simplex(S,k),x)) {

                ch_simplex T = new_simplex();
                NewSimplices.append(T);

                /* set the vertices of T as described above */
                for (int i = 1; i < dcur; i++) {
                  if ( i != k )
                    assoc_vertex_with_simplex(T,i,vertex_of_simplex(S,i)); }
                  if (k != dcur)
                    assoc_vertex_with_simplex(T,k,vertex_of_simplex(S,dcur));
                  assoc_vertex_with_simplex(T,dcur,z);
                  assoc_vertex_with_simplex(T,0,anti_origin);
                /* in the above, it is tempting to drop the tests ( i != k ) and ( k
                   != dcur ) since the subsequent lines after will correct the
                   erroneous assignment.  This reasoning is fallacious as the
                   procedure assoc_vertex_with_simplex also the internal data of the
                   third argument. */

                /* compute the equation of its base facet */
                  compute_equation_of_base_facet(T);

                /* record adjacency information for the two known neighbors */
                  set_neighbor(T,dcur,opposite_simplex(S,k),
                               index_of_vertex_in_opposite_simplex(S,k));
                  set_neighbor(T,0,S,k);


              }
            }
          }

          number_of_unbounded_simplices += NewSimplices.length();
          if ( vertex_of_simplex(start_facet,0) )
            start_facet = NewSimplices.back();
          ASSERT(!vertex_of_simplex(start_facet,0), start_facet_update_ndj);


          ch_simplex Af;
          forall(Af, NewSimplices) {
            for (int k = 1; k < dcur ; k++) {
              // neighbors 0 and dcur are already known
              if (opposite_simplex(Af,k) == nil) {
                // we have not performed the walk in the opposite direction yet
                ch_simplex T = opposite_simplex(Af,0);
                int y1 = 0;
                while ( vertex_of_simplex(T,y1) != vertex_of_simplex(Af,k) )
                  y1++;
                // exercise: show that we can also start with y1 = 1
                int y2 = index_of_vertex_in_opposite_simplex(Af,0);

                while ( vertex_of_simplex(T,0) == z ) {
                  // while T has peak x do find new y_1 */
                  int new_y1 = 0;
                  while (vertex_of_simplex(opposite_simplex(T,y1),new_y1) !=
                         vertex_of_simplex(T,y2))
                    new_y1++;
                    // exercise: show that we can also start with new_y1 = 1
                  y2 = index_of_vertex_in_opposite_simplex(T,y1);
                  T =  opposite_simplex(T,y1);
                  y1 = new_y1;
                }
                set_neighbor(Af,k,T,y1); // update adjacency information
              }
            }
          }


        }
    }

  }
#ifdef SELFCHECK
  check();
#endif
  return z;
}


template <class CHTRAITS,class POINT,class PLANE>
void chull<CHTRAITS,POINT,PLANE>::
visibility_search(CHSIMPLEX S, const POINT& x,
                  list< CHSIMPLEX >& visible_simplices,
                  int& number_of_visited_simplices,
                  int& location,
                  CHSIMPLEX& f) const
{
  number_of_visited_simplices ++;
  S->visited = true; // we have visited S and never come back ...
  for(int i = 0; i <= dcur; i++) {
    ch_simplex T = opposite_simplex(S,i); // for all neighbors T of S
    if (!T->visited ) {
      int side = CHTRAITS::side(T->base_facet,x);
      if ( is_unbounded_simplex(T) ) {
        if ( side > 0 ) {
          // T is an unbounded simplex with x-visible base facet
          visible_simplices.push(T);
          location = 1;
        }
        if ( side == 0 && location == -1 && contains_in_base_facet(T,x) ) {
          location = 0;
          f = T;
          return;
        }
      }
      if ( side > 0 || (side == 0 && location == -1) ) {
        visibility_search(T,x,visible_simplices,
                          number_of_visited_simplices,location,f);
        // do the recursive search
      }
    } // end visited
    else {
    }
  } // end for
}


template <class CHTRAITS,class POINT,class PLANE>
void chull<CHTRAITS,POINT,PLANE>::
clear_visited_marks(ch_Simplex<CHTRAITS,POINT,PLANE>* S) const
{
  S->visited = false; // clear the visit - bit
  for(int i = 0; i <= dcur; i++) // for all neighbors of S
    if (opposite_simplex(S,i)->visited)
      // if the i - th neighbor has been visited
      clear_visited_marks(opposite_simplex(S,i));
      // clear its bit recursively
}


template <class CHTRAITS,class POINT,class PLANE>
list< ch_Simplex<CHTRAITS,POINT,PLANE>* > chull<CHTRAITS,POINT,PLANE>::
facets_visible_from(const POINT& x)
/* returns the list of all facets that are visible from |x|. */
{
  list<ch_simplex> visible_simplices;
  int location = -1;                       // intialization is important
  int number_of_visited_simplices = 0;     // irrelevant
  ch_facet f;                              // irrelevant

  visibility_search(origin_simplex, x, visible_simplices,
                   number_of_visited_simplices, location, f);
  clear_visited_marks(origin_simplex);
  return visible_simplices;
}

template <class CHTRAITS,class POINT,class PLANE>
int chull<CHTRAITS,POINT,PLANE>::
is_where(const POINT& x)
{
  if ( is_dimension_jump(x) ) return 1;

  list<ch_simplex> visible_simplices;
  int location = -1;                       // intialization is important
  int number_of_visited_simplices = 0;     // irrelevant
  ch_facet f;

  visibility_search(origin_simplex, x, visible_simplices,
                   number_of_visited_simplices, location, f);
  clear_visited_marks(origin_simplex);
  return location;
}


template <class CHTRAITS,class POINT,class PLANE>
void chull<CHTRAITS,POINT,PLANE>::
dimension_jump(CHSIMPLEX S, rc_Vertex<CHTRAITS,POINT>* x)
{
  ch_simplex S_new;

  S->visited = true;
  assoc_vertex_with_simplex(S,dcur,x);
  if (! is_unbounded_simplex(S) ) {
    // S is bounded
    S_new = new_simplex();
    set_neighbor(S,dcur,S_new,0);
    assoc_vertex_with_simplex(S_new,0,anti_origin);
    for (int k = 1; k <= dcur; k++) {
      assoc_vertex_with_simplex(S_new,k,vertex_of_simplex(S,k-1));
    }

  }
  /* visit unvisited neighbors 0 to dcur - 1 */
  for (int k = 0; k <= dcur - 1; k++) {
    if (! opposite_simplex(S,k) -> visited) {
      dimension_jump(opposite_simplex(S,k), x);
    }
  }
  /* the recursive calls ensure that all neighbors exist */
  if ( is_unbounded_simplex(S) ) {
    set_neighbor(S,dcur,opposite_simplex(opposite_simplex(S,0),dcur),
                 index_of_vertex_in_opposite_simplex(S,0) + 1);

  } else {
    for (int k = 0; k < dcur; k++) {
      if ( is_unbounded_simplex(opposite_simplex(S,k)) ) {
        // if F' is unbounded
        set_neighbor(S_new,k + 1,opposite_simplex(S,k),dcur);
        // the neighbor of S_new opposite to v is S' and x stands in position dcur
      } else { // F' is bounded
        set_neighbor(S_new,k + 1,opposite_simplex(opposite_simplex(S,k),dcur),
                     index_of_vertex_in_opposite_simplex(S,k) + 1);
        // neighbor of S_new opposite to v is S_new'
        // the vertex opposite to v remains the same but ...
        // remember the shifting of the vertices one step to the right
      }
    }

  }
}



template <class CHTRAITS,class POINT,class PLANE>
void chull<CHTRAITS,POINT,PLANE>::
check() const
{
  check_topology();
  if (dcur < 1) return;

  /* Recall that center() gives us the center-point of the origin
     simplex. We check whether it is locally inside with respect to
     all hull facets.  */

  IPOINT centerpoint = center();
  list_item it;

  forall_items(it, all_simps) {
    ch_simplex S = (ch_simplex) all_simps[it];
    if ( is_unbounded_simplex(S) &&
         CHTRAITS::side(S->base_facet,centerpoint) >= 0) {
      error_handler(1,"check: center on wrong side of a hull facet");
    }
  }

  /* next we check convexity at every ridge. Let |S| be any hull
     simplex and let |v| be any vertex of its base facet. The vertex
     opposite to |v| must not be on the positive side of the base
     facet.*/

  forall_items(it,all_simps) {
    ch_simplex S = (ch_simplex) all_simps[it];
    if ( is_unbounded_simplex(S) ) {
      for (int i = 1; i <= dcur; i++) {
        int k = index_of_vertex_in_opposite_simplex(S,i);
        if (CHTRAITS::side(S->base_facet,
                         point_of_simplex(opposite_simplex(S,i),k)) > 0) {
          error_handler(1,"check: detected local non-convexity.");
        }
      }
    }
  }

  /* next we select one hull facet */

  ch_simplex selected_hull_simplex;
  forall_items(it,all_simps) {
    ch_simplex S = (ch_simplex) all_simps[it];
    if ( is_unbounded_simplex(S) )
     { selected_hull_simplex = S; break; }
  }

  /* we compute the center of gravity of the base facet of the hull simplex */

  VECTOR center_of_hull_facet = CHTRAITS::zero_vector(dmax);
  for (int i = 1; i <= dcur; i++) {
    center_of_hull_facet +=
    CHTRAITS::to_vector(point_of_simplex(selected_hull_simplex,i));
  }
  IPOINT center_of_hull_facet_point =
    CHTRAITS::to_ipoint(center_of_hull_facet/RT(dcur));

  /* we set up the ray from the center to the center of the hull facet */

  IRAY l(centerpoint, center_of_hull_facet_point);

  /* and check whether it intersects the interior of any hull facet */

  forall_items(it,all_simps) {
    ch_simplex S = (ch_simplex) all_simps[it];
    if ( is_unbounded_simplex(S) && S != selected_hull_simplex ) {
      IPOINT p;
      PLANE  h = S->base_facet;
      if ( CHTRAITS::intersection(h,l,p) ) {
        array<POINT> A(dcur);
        for (int i = 0; i < dcur; i++) {
          A[i] = point_of_simplex(S,i + 1);
        }
        if (CHTRAITS::contained_in_simplex(A,p))
          error_handler(1,"check: current hull has double coverage.");
      }
    }
  }
}


template <class CHTRAITS,class POINT,class PLANE>
list< ch_Simplex<CHTRAITS,POINT,PLANE>* >  chull<CHTRAITS,POINT,PLANE>::
all_facets() const
{
  list<ch_facet> result;
  if (dcur > 1) {
    stack<ch_facet> candidates;
    candidates.push(start_facet);
    start_facet -> visited = true;
    ch_facet current;
    while ( !candidates.empty() ) {
      current = candidates.pop();
      ASSERT(!vertex_of_simplex(current,0), all_facets);
      result.append(current);
      for(int i = 1; i <= dcur; ++i) {
        ch_facet f = opposite_simplex(current,i);
        if (!f -> visited) {
          candidates.push(f);
          f -> visited = true;
        }
      }
    }
    clear_visited_marks(start_facet);
  }
  else if ( dcur == 1 ) {
    result.append(start_facet);
  }
  return result;
}


template <class CHTRAITS,class POINT,class PLANE> void
d2_show(const chull<CHTRAITS,POINT,PLANE>& C,window& W)
{ /* We first draw every simplex*/
  typedef chull<CHTRAITS,POINT,PLANE>::ch_simplex ch_simplex;
  typedef chull<CHTRAITS,POINT,PLANE>::ch_vertex ch_vertex;
  ch_simplex S;

#ifdef DDGEO_STL_ITERATORS
  chull<CHTRAITS,POINT,PLANE>::ch_simplex_iterator sit;
  for(sit =  C.simplices_begin();
      sit != C.simplices_end();
      ++sit) {
    S = *sit;
#else
  list<ch_simplex> Simpl = C.all_simplices();
  forall(S,Simpl) {
#endif
    for (int v = ( C.is_unbounded_simplex(S)  ? 1 : 0);
         v <= C.dcurrent(); v++) {
      // for each vertex except the anti - origin

      for (int e = v + 1; e <= C.dcurrent(); e++) {
        // draw undrawn edges incident to vertex
        if ( C.is_unbounded_simplex(S) )
          W.set_line_width(3); // thick lines for hull edges
        else
          W.set_line_width(1);
        W.draw_segment(CHTRAITS::to_d2_point(C.point_of_simplex(S,v)),
                       CHTRAITS::to_d2_point(C.point_of_simplex(S,e)));
      }
    }
  }
  /* Now we draw every point */
  POINT x;
#ifdef DDGEO_STL_ITERATORS
  chull<CHTRAITS,POINT,PLANE>::ch_point_iterator pit;
  for(pit =  C.points_begin(); pit != C.points_end(); ++pit) {
    x=*pit;
#else
  forall(x,C.all_points()) {
#endif
   W.draw_point(CHTRAITS::to_d2_point(x));
  }
}



template <class CHTRAITS,class POINT,class PLANE>
void
d3_surface_map(const chull<CHTRAITS,POINT,PLANE>& C, GRAPH<POINT,int>& G)
{
  typedef chull<CHTRAITS,POINT,PLANE>::ch_vertex  ch_vertex;
  typedef chull<CHTRAITS,POINT,PLANE>::ch_simplex ch_simplex;
  typedef chull<CHTRAITS,POINT,PLANE>::ch_facet   ch_facet;
  typedef typename CHTRAITS::RT RT;
  typedef typename CHTRAITS::IPOINT CHIPOINT;
  typedef typename CHTRAITS::PLANE  CHPLANE;
#ifdef DDGEO_STL_ITERATORS
  chull<CHTRAITS,POINT,PLANE>::ch_vertex_iterator vit;
  chull<CHTRAITS,POINT,PLANE>::ch_simplex_iterator sit;
#else
  list<ch_simplex> Simps = C.all_simplices();
  list<ch_vertex>  Verts = C.all_vertices();
#endif
  G.clear();
  if (C.dim() != 3)
    error_handler(1,"d3_surface_map: dim must be 3");
  if (C.dcurrent() < 3) {
    map<ch_vertex,node> node_for(nil);
    ch_vertex v;
    ch_simplex s;
#ifdef DDGEO_STL_ITERATORS
    for(vit  = C.vertices_begin();
        vit != C.vertices_end();
        ++vit) {
      v = *vit;
#else
    forall(v, Verts) {
#endif
      node_for[v] = G.new_node(C.associated_point(v));
    }
    if (C.dcurrent() <= 0)
    { return; }
    if (C.dcurrent() == 1) {
#ifdef DDGEO_STL_ITERATORS
      for(sit =  C.simplices_begin();
          sit != C.simplices_end();
          ++sit) {
        s = *sit;
#else
      forall(s,Simps) {
#endif
        if (C.is_bounded_simplex(s)) {
          node v0 = node_for[C.vertex(s,0)];
          node v1 = node_for[C.vertex(s,1)];
          edge e01 = G.new_edge(v0,v1);
          edge e10 = G.new_edge(v1,v0);
          G.set_reversal(e01,e10);
        }
      }
      return;
    }

    
    if (C.dcurrent() == 2)
    {
      node_array<bool> untreated(G,true);
    /*
      typename CHTRAITS::IPOINT pc =
        CHTRAITS::to_ipoint(C.qcenter()/RT(3) +
        CHTRAITS::normal(C.base_facet_plane(C.osimplex())));
    */
      CHIPOINT pc =
        CHTRAITS::to_ipoint(( C.qcenter()/RT(3) )
      + CHTRAITS::normal(CHPLANE(
                           C.associated_point( C.vertex( C.osimplex(),0)),
                           C.associated_point( C.vertex( C.osimplex(),1)),
                           C.associated_point( C.vertex( C.osimplex(),2))))
                           );
    
      forall(s,Simps)
      {
        if (C.is_bounded_simplex(s))
        {
          for (int i = 0; i <= 2; i++)
          {
            node vi = node_for[C.vertex(s,i)];
            if ( untreated[vi] )
            {
              int j = (i + 1) % 3;
              // a vertex different from i;
              int k = (i + 2) % 3;
              node vj = node_for[C.vertex(s,j)];
              node vk = node_for[C.vertex(s,k)];
              if (CHTRAITS::orientation(G[vi],G[vj],G[vk],pc) > 0)
              {
                leda_swap(vj,vk);
                leda_swap(j,k);
              }
    
              edge efirst = G.new_edge(vi,vk);  // first edge incident to vi
              ch_simplex scur = s;
              int jcur = j;
              int kcur = k;
              int icur = i;
    
              while (    C.is_bounded_simplex(C.opposite_simplex(scur,jcur))
                      && C.opposite_simplex(scur,jcur) != s )
              {
                // we have not reached the end nor closed the cycle
                kcur = C.index_of_opposite_vertex(scur,jcur);
                scur = C.opposite_simplex(scur,jcur);
                for (icur = 0; icur <= 2; icur++)
                { if (node_for[C.vertex(scur,icur)] == vi) break; }
                jcur = 3 - icur - kcur;
                vk = node_for[C.vertex(scur,kcur)];
                G.new_edge(vi,vk);
              }
    
              if (C.is_unbounded_simplex(C.opposite_simplex(scur,jcur)))
              {
                /* we also need to walk in the other direction */
    
                efirst = G.new_edge(efirst,vj,0,LEDA_PREFIXLI before);  // 0 is etype
                scur = s;
                jcur = j;
                kcur = k;
                icur = i;
                // restore initial situation
    
                while ( C.is_bounded_simplex(C.opposite_simplex(scur,kcur)))
                {
                  // we have not reached the end
                  jcur = C.index_of_opposite_vertex(scur,kcur);
                  scur = C.opposite_simplex(scur,kcur);
                  for (icur = 0; icur <= 2; icur++)
                  { if (node_for[C.vertex(scur,icur)] == vi) break; }
                  kcur = 3 - jcur -icur;
                  vj = node_for[C.vertex(scur,jcur)];
                  efirst = G.new_edge(efirst,vj,0,LEDA_PREFIXLI before); //as above
                } //end while
              } // end if
              untreated[vi] = false;
            } // end if untreated
          } // end for i
        } // end if bounded
      } // end forall
    if (!G.make_map())
      error_handler(1,"chull::surface_graphrep: not bidirected");
    return;
    }
    
  }

  ch_facet f;
  ch_vertex v;
  map<ch_vertex,node> node_for(nil);
  int facet_num = 0;
#ifdef DDGEO_STL_ITERATORS
  chull<CHTRAITS,POINT,PLANE>::ch_facet_iterator fit;
  for(fit  = C.facets_begin();
      fit != C.facets_end();
      ++fit) {
    f = *fit;
#else
  list<ch_facet> Surface = C.all_facets();
  forall(f,Surface) {
#endif
    ++facet_num;
    for (int i=0; i < C.dcurrent(); i++) {
      v = C.vertex_of_facet(f,i);
      if (!node_for[v])
        node_for[v] = G.new_node(C.associated_point(v));
    }
  }
  if ( 2*G.number_of_nodes() != facet_num + 4)
    error_handler(1,"d3_surface_map: node equation wrong.");

  node_array<bool> untreated(G,true);
#ifdef DDGEO_STL_ITERATORS
  for(fit  = C.facets_begin();
      fit != C.facets_end();
      ++fit) {
    f = *fit;
#else
  forall(f,Surface) {
#endif
    for (int i = 0; i < C.dcurrent(); i++) {
      node vi = node_for[C.vertex_of_facet(f,i)];
      if (untreated[vi]) {
        untreated[vi] = false;
        int j = (i + 1) % C.dcurrent();
        // a vertex different from i;
        int k = (i + 2) % C.dcurrent();
        node vj = node_for[C.vertex_of_facet(f,j)];
        node vk = node_for[C.vertex_of_facet(f,k)];
        if (CHTRAITS::orientation(G[vi],G[vj],G[vk],C.center()) > 0) {
          leda_swap(vk,vj);
          leda_swap(j,k);
        }

        edge efirst = G.new_edge(vi,vk);  // first edge incident to vi
        ch_facet fcur = f;
        int jcur = j; int kcur = k; int icur = i;

        while ( C.opposite_facet(fcur,jcur) != f ) {
          // we have not reached the end
          kcur = C.index_of_vertex_in_opposite_facet(fcur,jcur);
          fcur = C.opposite_facet(fcur,jcur);
          for (icur = 0; icur < 3; icur++)
             if (node_for[C.vertex_of_facet(fcur,icur)] == vi) break;
          jcur = 3 - icur - kcur;
          vk = node_for[C.vertex_of_facet(fcur,kcur)];
          G.new_edge(vi,vk);
        }

      } // end if untreated
    } // end for i
  } // end forall
  if (G.number_of_edges() != (3*facet_num))
    error_handler(1,"chull::surface_graphrep: wrong number of edges");
  if (!G.make_map())
    error_handler(1,"chull::surface_graphrep: not bidirected");
}


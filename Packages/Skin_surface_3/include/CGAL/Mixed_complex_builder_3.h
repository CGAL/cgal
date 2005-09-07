#ifndef MIXED_COMPLEX_BUILDER_H
#define MIXED_COMPLEX_BUILDER_H

#include <CGAL/Unique_hash_map.h>
#include <CGAL/Compute_anchor_3.h>
#include <CGAL/Simplex_3.h>

#include <map>

CGAL_BEGIN_NAMESPACE

template < class Skin_Traits_ >
class Mixed_complex_builder_3 {
public:
  typedef Mixed_complex_builder_3<Skin_Traits_> Self;

  typedef Skin_Traits_                          Skin_Traits;
  typedef typename Skin_Traits::Regular         Regular;
  typedef typename Skin_Traits::Simplicial      Simplicial;

private:
  typedef typename Regular::Vertex_handle            Rt_Vertex_handle;
  typedef typename Regular::Edge                     Rt_Edge;
  typedef typename Regular::Facet                    Rt_Facet;
  typedef typename Regular::Cell_handle              Rt_Cell_handle;
	
  typedef typename Regular::Finite_vertices_iterator Rt_Finite_vertices_iterator;
  typedef typename Regular::Finite_edges_iterator    Rt_Finite_edges_iterator;
  typedef typename Regular::Finite_facets_iterator   Rt_Finite_facets_iterator;
  typedef typename Regular::All_cells_iterator       Rt_All_cells_iterator;
  typedef typename Regular::Finite_cells_iterator    Rt_Finite_cells_iterator;
	
  typedef typename Regular::Cell_circulator          Rt_Cell_circulator;

  typedef Simplex_3<Regular>                    Rt_Simplex;
  typedef typename Skin_Traits::Regular_K       Rt_Geom_traits;
  typedef typename Regular::Bare_point          Rt_Point;
  typedef typename Regular::Geom_traits::RT     Rt_RT;
  typedef typename Regular::Weighted_point      Rt_Weighted_point;

  typedef typename Simplicial::Vertex_handle    Sc_Vertex_handle;
  typedef typename Simplicial::Edge             Sc_Edge;
  typedef typename Simplicial::Facet            Sc_Facet;
  typedef typename Simplicial::Cell_handle      Sc_Cell_handle;

  typedef typename Simplicial::Finite_vertices_iterator Sc_Finite_vertices_iterator;
  typedef typename Simplicial::Finite_edges_iterator    Sc_Finite_edges_iterator;
  typedef typename Simplicial::Finite_facets_iterator   Sc_Finite_facets_iterator;
  typedef typename Simplicial::All_cells_iterator       Sc_All_cells_iterator;
  typedef typename Simplicial::Finite_cells_iterator    Sc_Finite_cells_iterator;
  typedef typename Simplicial::Cell_circulator          Sc_Cell_circulator;
	
  typedef typename Skin_Traits::Simplicial_K            Sc_Geom_traits;
  typedef typename Sc_Geom_traits::Point_3              Sc_Point;
  typedef typename Sc_Geom_traits::RT                   Sc_RT;
  //   typedef Weighted_point<Sc_Point>                Sc_Weighted_point;

  typedef Triangulation_incremental_builder_3<Simplicial>
                                             Triangulation_incremental_builder;

  //   typedef typename Skin_Traits::Mesh_K               M_traits;
  //   typedef Skin_surface_quadratic_surface_3< M_traits > Quadratic_surface;
  //   typedef Skin_surface_sphere_3< M_traits >            Sphere_surface;
  //   typedef Skin_surface_hyperboloid_3< M_traits >       Hyperboloid_surface;
  //   typedef typename M_traits::RT                        M_RT;

private:
  Sc_Vertex_handle add_vertex(Rt_Simplex &sDel, Rt_Simplex &sVor); 
  Sc_Cell_handle add_cell(Sc_Vertex_handle vh[], int orient, Rt_Simplex s);
	
  Sc_Vertex_handle get_vertex(Rt_Simplex &sDel, Rt_Simplex &sVor);


  void construct_vertices();
  Sc_Point get_anchor(Rt_Simplex &sDel, Rt_Simplex &sVor);

  void construct_0_cells();
  void construct_1_cells();
  void construct_2_cells();
  void construct_3_cells();
  void compute_value_on_vertices();

  typedef std::map<Rt_Simplex, typename Sc_Geom_traits::Point_3> Focus_map;
  Focus_map foc_map;

public:
  Skin_Traits traits;
	
  Mixed_complex_builder_3(Regular &T, Simplicial &sc, double shrink)
    : T(T), sc(sc), triangulation_incr_builder(sc), shrink(shrink), compute_anchor_obj(T) {
    edge_index[0][0] = -1; edge_index[0][1] =  0;
    edge_index[0][2] =  1; edge_index[0][3] =  2;
    edge_index[1][0] =  0; edge_index[1][1] = -1;
    edge_index[1][2] =  3; edge_index[1][3] =  4;
    edge_index[2][0] =  1; edge_index[2][1] =  3;
    edge_index[2][2] = -1; edge_index[2][3] =  5;
    edge_index[3][0] =  2; edge_index[3][1] =  4;
    edge_index[3][2] =  5; edge_index[3][3] = -1;
	  
    assert((0<shrink) && (shrink<1));
    //     shrink_m = shrink;

    build();
  }
private:
  void build() {

    triangulation_incr_builder.begin_triangulation(3);

    construct_vertices();

    construct_0_cells(); // mixed cells corresponding to regular vertices
    construct_1_cells(); // mixed cells corresponding to regular edges
    construct_2_cells(); // mixed cells corresponding to regular facets
    construct_3_cells(); // mixed cells corresponding to regular cells

    triangulation_incr_builder.end_triangulation();

    compute_value_on_vertices();
  }
  // Reference back to the regular triangulation
  std::map<Sc_Cell_handle, Rt_Simplex> backRef;
	
private:
  Regular &T;
  Simplicial &sc;
  Triangulation_incremental_builder triangulation_incr_builder;
  Sc_RT shrink;
  // M_RT shrink_m;

    Construct_weighted_circumcenter_3<
      Regular_triangulation_euclidean_traits_3<
      Sc_Geom_traits> >                                orthocenter_obj;
    Compute_squared_radius_smallest_orthogonal_sphere_3<
      Regular_triangulation_euclidean_traits_3<
      typename Simplicial::Geom_traits> >                                orthoweight_obj;
  Compute_anchor_3<Regular> compute_anchor_obj;

  int edge_index[4][4];
    struct Index_c4 { Sc_Vertex_handle V[4]; };
    struct Index_c6 { Sc_Vertex_handle V[6]; };
    struct Index_c44 { Sc_Vertex_handle V[4][4]; };
    struct Index_v {
      Unique_hash_map < Rt_Vertex_handle, Sc_Vertex_handle > V;
    };
    // Facets on the border of the simplicial complex:
    // name is given by (dim del,dim vor)

    // index to vertex
    Unique_hash_map < Rt_Cell_handle, Index_c4 > index_03;
    // index to edge
    Unique_hash_map < Rt_Cell_handle, Index_c6 > index_13;
    // index to facet
    Unique_hash_map < Rt_Cell_handle, Index_c4 > index_23;
    // first index points to facet, second to vertex
    Unique_hash_map < Rt_Cell_handle, Index_c44 > index_02;
    // From the vertex to the adjacent vertex
    Unique_hash_map < Rt_Vertex_handle, Index_v > index_01;
    // One for every cell
    Unique_hash_map < Rt_Cell_handle, Sc_Vertex_handle > index_33;
    // One for every vertex
    Unique_hash_map < Rt_Vertex_handle, Sc_Vertex_handle > index_00;
    // First index points to the facet,
    // the second to the vertex that is removed to obtain the edge
    Unique_hash_map < Rt_Cell_handle, Index_c44 > index_12;
    // index to the facet
    Unique_hash_map < Rt_Cell_handle, Index_c4 > index_22;
    // index to the vertex
    Unique_hash_map < Rt_Vertex_handle, Index_v > index_11;

    std::list<Sc_Vertex_handle> degenerate_vertices;
};

// Constructs the vertices of the simplicial complex
template <class Skin_Traits_>
void
Mixed_complex_builder_3<Skin_Traits_>::construct_vertices() {
  Rt_All_cells_iterator acit;
  Rt_Finite_cells_iterator cit;
  Rt_Finite_facets_iterator fit;
  Rt_Finite_edges_iterator eit;
  Rt_Finite_vertices_iterator vit;
  Rt_Cell_circulator ccir, cstart;
  Rt_Vertex_handle v1, v2, v3;
  Rt_Edge e;
  Rt_Cell_handle c1, c2;
  Rt_Simplex sDel, sVor;
  Sc_Vertex_handle vh;

  // anchor dimDel=0, dimVor=3
  for (cit=T.finite_cells_begin(); cit!=T.finite_cells_end(); cit++) {
    sVor = compute_anchor_obj.anchor_vor(cit);
    for (int i=0; i<4; i++) {
      sDel = compute_anchor_obj.anchor_del(cit->vertex(i));
      assert(sDel.dimension() == 0);
      assert(sVor.dimension() == 3);
      vh = add_vertex(sDel, sVor);
      index_03[cit].V[i] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }

  // anchor dimDel=1, dimVor=3
  for (cit=T.finite_cells_begin(); cit!=T.finite_cells_end(); cit++) {
    sVor = compute_anchor_obj.anchor_vor(cit);
    assert(sVor == Rt_Simplex(cit));
    for (int i=0; i<3; i++) {
      for (int j=i+1; j<4; j++) {
	sDel = compute_anchor_obj.anchor_del(Rt_Edge(cit,i,j));
	if (sDel.dimension() == 1) {
	  vh = add_vertex(sDel, sVor);
	} else {
	  vh = get_vertex(sDel, sVor);
	}
	index_13[cit].V[edge_index[i][j]] = vh;
	assert(vh == get_vertex(sDel, sVor));
      }
    }
  }

  // anchor dimDel=2, dimVor=3 and dimDel=0, dimVor=2
  for (fit=T.finite_facets_begin(); fit!=T.finite_facets_end(); fit++) {
    // anchor dimDel=2, dimVor=3
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);
    int n_index = c2->index(c1);

    sDel = compute_anchor_obj.anchor_del(*fit);
    assert(sDel.dimension() == compute_anchor_obj.anchor_del(Rt_Facet(c2,n_index)).dimension());
    assert(sDel == compute_anchor_obj.anchor_del(Rt_Facet(c2,n_index)));
    if (!T.is_infinite(c1)) {
      sVor = compute_anchor_obj.anchor_vor(c1);
      if (sDel.dimension() == 2) {
        vh = add_vertex(sDel, sVor);
      } else {
	vh = get_vertex(sDel, sVor);
      }
      index_23[c1].V[fit->second] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
    if (!T.is_infinite(c2)) {
      sVor = compute_anchor_obj.anchor_vor(c2);
      if (sDel.dimension() == 2) {
        vh = add_vertex(sDel, sVor);
      } else {
	vh = get_vertex(sDel, sVor);
      }
      index_23[c2].V[n_index] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
    // anchor dimDel=0, dimVor=2
    sVor = compute_anchor_obj.anchor_vor(*fit);
    for (int i=1; i<4; i++) {
      v1 = c1->vertex((fit->second+i)&3);
      sDel = compute_anchor_obj.anchor_del(v1);
      if (sVor.dimension() == 2) {
        vh = add_vertex(sDel, sVor);
      } else {
        vh = get_vertex(sDel, sVor);
      }
      index_02[c1].V[fit->second][(fit->second+i)&3] = vh;
      index_02[c2].V[n_index][c2->index(v1)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }
	
  // anchor dimDel=0, dimVor=1
  for (eit=T.finite_edges_begin(); eit!=T.finite_edges_end(); eit++) {
    sVor = compute_anchor_obj.anchor_vor(*eit);

    v1 = eit->first->vertex(eit->second);
    v2 = eit->first->vertex(eit->third);
    sDel = compute_anchor_obj.anchor_del(v1);
    if (sVor.dimension() == 1) {
      vh = add_vertex(sDel, sVor);
    } else {
      vh = get_vertex(sDel, sVor);
    }
    index_01[v1].V[v2] = vh;
    assert(vh == get_vertex(sDel, sVor));
			
    sDel = compute_anchor_obj.anchor_del(v2);
    if (sVor.dimension() == 1) {
      vh = add_vertex(sDel, sVor);
    } else {
      vh = get_vertex(sDel, sVor);
    }
    index_01[v2].V[v1] = vh;
    assert(vh == get_vertex(sDel, sVor));
  }
	
  // anchor dimDel=3, dimVor=3
  for (cit=T.finite_cells_begin(); cit!=T.finite_cells_end(); cit++) {
    sDel = compute_anchor_obj.anchor_del(cit);
    sVor = compute_anchor_obj.anchor_vor(cit);
    if (sDel.dimension() == 3) {
      vh = add_vertex(sDel, sVor);
    } else {
      vh = get_vertex(sDel, sVor);
    }
    index_33[cit] = vh;
    assert(vh == get_vertex(sDel, sVor));
  }

  // anchor dimDel=0, dimVor=0
  for (vit=T.finite_vertices_begin(); vit!=T.finite_vertices_end(); vit++) {
    sDel = compute_anchor_obj.anchor_del(vit);
    sVor = compute_anchor_obj.anchor_vor(vit);
    if (sVor.dimension() == 0) {
      vh = add_vertex(sDel, sVor);
    } else {
      vh = get_vertex(sDel, sVor);
    }
    index_00[vit] = vh;
    assert(vh == get_vertex(sDel, sVor));
  }
	
  // anchor dimDel=1, dimVor=2
  for (fit=T.finite_facets_begin(); fit!=T.finite_facets_end(); fit++) {
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);
    int n_index = c2->index(c1);

    sVor = compute_anchor_obj.anchor_vor(*fit);
    for (int i=1; i<3; i++) {
      for (int j=i+1; j<4; j++) {
        e.first = c1;
        e.second = (fit->second+i)&3;
        e.third = (fit->second+j)&3;
        sDel = compute_anchor_obj.anchor_del(e);
        if (sDel.dimension() == 1 && sVor.dimension() == 2) {
          vh = add_vertex(sDel, sVor);
        } else {
          vh = get_vertex(sDel, sVor);
        }
	int index = 6 - fit->second - e.second - e.third;
	index_12[c1].V[fit->second][index] = vh;
	index = 6-n_index-
	  c2->index(c1->vertex(e.second))-
	  c2->index(c1->vertex(e.third));
	index_12[c2].V[n_index][index] = vh;
	assert(vh == get_vertex(sDel, sVor));
      }
    }
  }
	
  // anchor dimDel=2, dimVor=2
  for (fit=T.finite_facets_begin(); fit!=T.finite_facets_end(); fit++) {
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);
    int n_index = c2->index(c1);

    sVor = compute_anchor_obj.anchor_vor(*fit);
    sDel = compute_anchor_obj.anchor_del(*fit);
    if (sDel.dimension() == 2 &&
        sVor.dimension() == 2) {
      vh = add_vertex(sDel, sVor);
    } else {
      vh = get_vertex(sDel, sVor);
    }
    index_22[c1].V[fit->second] = vh;
    index_22[c2].V[n_index] = vh;
    assert(vh == get_vertex(sDel, sVor));
  }
	
  // anchor dimDel=1, dimVor=1
  for (eit=T.finite_edges_begin(); eit!=T.finite_edges_end(); eit++) {
    v1 = eit->first->vertex(eit->second);
    v2 = eit->first->vertex(eit->third);

    sVor = compute_anchor_obj.anchor_vor(*eit);
    sDel = compute_anchor_obj.anchor_del(*eit);

    if (sDel.dimension() == 1 &&
        sVor.dimension() == 1) {
      vh = add_vertex(sDel, sVor);
    } else {
      vh = get_vertex(sDel, sVor);
    }
    
    index_11[v1].V[v2] = vh;
    index_11[v2].V[v1] = vh;
    assert(vh == get_vertex(sDel, sVor));
  }
}

// Constructs the cells of the mixed complex corresponding
// to Regular vertices
template <class Skin_Traits_>
void
Mixed_complex_builder_3<Skin_Traits_>::construct_0_cells() {
  Rt_Simplex sDel, sVor;
  Sc_Vertex_handle vh[4];
  
//   Quadratic_surface *surf;
  for (Rt_Finite_vertices_iterator vit=T.finite_vertices_begin();
       vit!=T.finite_vertices_end(); vit++) {
//     surf = new Sphere_surface(traits.r2m_converter(vit->point()), 
//       shrink_m, 1);
    
    Rt_Simplex simplex(vit);
    sDel = compute_anchor_obj.anchor_del(vit);
    sVor = compute_anchor_obj.anchor_vor(vit);
    vh[0] = get_vertex(sDel,sVor);
    
    std::list<Rt_Cell_handle> adj_cells;
    typename std::list<Rt_Cell_handle>::iterator adj_cell;
    T.incident_cells(vit, std::back_inserter(adj_cells));
    
    // Construct cells:
    for (adj_cell = adj_cells.begin();
	 adj_cell != adj_cells.end();
	 adj_cell ++) {
      if (!T.is_infinite(*adj_cell)) {
	sVor = compute_anchor_obj.anchor_vor(*adj_cell);
	vh[3] = get_vertex(sDel,sVor);
	int index = (*adj_cell)->index(vit);
	for (int i=1; i<4; i++) {
	  sVor = compute_anchor_obj.anchor_vor(
	    Rt_Facet(*adj_cell,(index+i)&3));
	  vh[2] = get_vertex(sDel,sVor);
	  
	  for (int j=1; j<4; j++) {
	    if (j!=i) {
	      sVor = compute_anchor_obj.anchor_vor(Rt_Edge(*adj_cell,index,(index+j)&3));
	      vh[1] = get_vertex(sDel,sVor);
	      if ((vh[0] != vh[1]) && (vh[1] != vh[2]) && (vh[2] != vh[3])) {
		
		Sc_Cell_handle ch;
		
		ch = add_cell(vh,(index + (j==(i%3+1)? 1:0))&1,simplex);
// 		ch->surf = surf;
	      }
	    }
	  }
	}
      }
    }
  }
}

// Constructs 1-cells of the mixed complex corresponding to edges
// of the regular triangulation
template <class Skin_Traits_>
void
Mixed_complex_builder_3<Skin_Traits_>::construct_1_cells() {
  Rt_Simplex sDel, sVor, sDel0;
  Sc_Vertex_handle vh[4];
  Rt_Vertex_handle v[2];
  Sc_Cell_handle ch;
  
//   Quadratic_surface *surf;
  for (Rt_Finite_edges_iterator eit=T.finite_edges_begin();
       eit!=T.finite_edges_end(); eit++) {
//     surf = new Hyperboloid_surface(
//       traits.s2m_converter(
// 	Sc_Weighted_point(
// 	  orthocenter_obj(
// 	    traits.r2s_converter(eit->first->vertex(eit->second)->point()),
// 	    traits.r2s_converter(eit->first->vertex(eit->third)->point())),
// 	  orthoweight_obj(
// 	    traits.r2s_converter(eit->first->vertex(eit->second)->point()),
// 	    traits.r2s_converter(eit->first->vertex(eit->third)->point())))),
//       traits.r2m_converter(eit->first->vertex(eit->second)->point())-
//       traits.r2m_converter(eit->first->vertex(eit->third)->point()),
//       shrink_m, 1);
    Rt_Simplex simplex(*eit);
    sDel = compute_anchor_obj.anchor_del(*eit);
    sVor = compute_anchor_obj.anchor_vor(*eit);
    Sc_Vertex_handle vh1 = get_vertex(sDel,sVor);
    
    v[0] = eit->first->vertex(eit->second);
    v[1] = eit->first->vertex(eit->third);
    
    // Construct cells on the side of v[0]:
    for (int vi=0; vi<2; vi++) {
      sDel0 = compute_anchor_obj.anchor_del(v[0]);
      if (!(sDel == sDel0)) {
	Rt_Cell_circulator ccir, cstart;
	ccir = cstart = T.incident_cells(*eit);
	do {
	  if (!T.is_infinite(ccir)) {
	    int index0 = ccir->index(v[0]);
	    int index1 = ccir->index(v[1]);

	    for (int fi=1; fi<4; fi++) {
	      if (((index0+fi)&3) != index1) {
		vh[0] = index_01[v[0]].V[v[1]];
		vh[1] = vh1;
		vh[2] = index_12[ccir].V[(index0+fi)&3]
		  [6-((index0+fi)&3)-index0-index1];
		vh[3] = index_13[ccir].V[edge_index[index0][index1]];
		if ((vh[0]!=vh[1]) && (vh[1]!=vh[2]) && (vh[2]!=vh[3])) {
		  int orient;
		  if (((4+index1-index0)&3) == 1) {
		    orient = (index1 + (fi==2))&1;
		  } else {
		    orient = (index1 + (fi==1))&1;
		  }
		  // vh: dimension are (01,11,12,13)
		  ch = add_cell(vh,orient,simplex);
// 		  ch->surf = surf;
									
		  vh[1] = index_02[ccir].V[(index0+fi)&3][index0];
		  // vh: dimension are (01,02,12,13)
		  ch = add_cell(vh,1-orient,simplex);
// 		  ch->surf = surf;
									
		  vh[2] = index_03[ccir].V[index0];
		  // vh: dimension are (01,02,03,13)
		  ch = add_cell(vh,orient,simplex);
// 		  ch->surf = surf;
		}
	      }
	    }
	  }
	  ccir ++;
	} while (ccir != cstart);
      }

      // Swap v[0] and v[1]:
      Rt_Vertex_handle tmp_v = v[0]; v[0]=v[1]; v[1]=tmp_v;
    }
		
  }
}


// Constructs 2-cells of the mixed complex corresponding to facets
// of the regular triangulation
template <class Skin_Traits_>
void
Mixed_complex_builder_3<Skin_Traits_>::construct_2_cells() {
  Rt_Simplex sDel, sVor;//, sDel0;
  Sc_Vertex_handle vh[4]; // Implicit function over vLabels is increasing ...
  // 	Rt_Vertex_handle v[2];
  Rt_Cell_handle rt_ch;
  Sc_Cell_handle ch;
//   Quadratic_surface *surf;
	
  for (Rt_Finite_facets_iterator fit = T.finite_facets_begin();
       fit != T.finite_facets_end();
       fit ++) {
//     surf = new Hyperboloid_surface(
//       traits.s2m_converter(
// 	Sc_Weighted_point(
// 	  orthocenter_obj(
// 	    traits.r2s_converter(fit->first->vertex((fit->second+1)&3)->point()),
// 	    traits.r2s_converter(fit->first->vertex((fit->second+2)&3)->point()),
// 	    traits.r2s_converter(fit->first->vertex((fit->second+3)&3)->point())),
// 	  orthoweight_obj(
// 	    traits.r2s_converter(fit->first->vertex((fit->second+1)&3)->point()),
// 	    traits.r2s_converter(fit->first->vertex((fit->second+2)&3)->point()),
// 	    traits.r2s_converter(fit->first->vertex((fit->second+3)&3)->point())))),
//       typename M_traits::Construct_orthogonal_vector_3()(
// 	traits.r2m_converter(fit->first->vertex((fit->second+1)&3)->point()),
// 	traits.r2m_converter(fit->first->vertex((fit->second+2)&3)->point()),
// 	traits.r2m_converter(fit->first->vertex((fit->second+3)&3)->point())),
//       1-shrink_m, -1);
    
    rt_ch = fit->first;
    int index = fit->second;
    Rt_Simplex simplex(*fit);
		
    for (int i=0; i<2; i++) { // Do this twice
      if (!T.is_infinite(rt_ch)) {
	vh[3] = index_23[rt_ch].V[index];
	Sc_Vertex_handle vh2 = index_22[rt_ch].V[index]; // Cache
	if (vh2 != vh[3]) {
	  // Facet and cell do not coincide ..
	  for (int vi=1; vi<4; vi++) {
	    vh[0] = index_02[rt_ch].V[index][(index+vi)&3];
	    for (int ei=1; ei<4; ei++) {
	      if (vi != ei) {
		vh[2] = vh2;
		assert(((index+6-vi-ei)&3) == (6-index-((index+vi)&3)-((index+ei)&3)));
		assert(((index+6-vi-ei)&3) == ((6-3*index-vi-ei)&3));
		vh[1] = index_12[rt_ch].V[index][(6+index-vi-ei)&3];
		if ((vh[0] != vh[1]) && (vh[1] != vh[2])) {
		  // index0: v0
		  // index1: v1
		  // index0+fi&3 == facet
		  int index0 = (index+vi)&3;
		  int index1 = (index+ei)&3;
		  int fi = (6+index-vi-ei)&3;//6-index-index0-index1;
		  int orient;
									
		  if (((4+index1-index0)&3) == 3) {
		    orient = (index1 + (((4+index0-fi)&3)==2))&1;
		  } else {
		    orient = (index1 + (((4+index0-fi)&3)==1))&1;
		  }

		  ch = add_cell(vh,orient,simplex);
// 		  ch->surf = surf;
									
		  vh[2] = index_13[rt_ch].V[edge_index[index0][index1]];
		  ch = add_cell(vh,1-orient,simplex);
// 		  ch->surf = surf;
									
		  vh[1] = index_03[rt_ch].V[index0];
		  ch = add_cell(vh,orient,simplex);
// 		  ch->surf = surf;
		} 
	      }
	    }
	  }
	}
      }
      // swap to the other cell
      Rt_Cell_handle ch_old = rt_ch;
      rt_ch = rt_ch->neighbor(index);
      index = rt_ch->index(ch_old);
    }

    CGAL_assertion(rt_ch == fit->first);
    CGAL_assertion(index == fit->second);
  }
}


// Constructs 3-cells of the mixed complex corresponding to cells
// of the regular triangulation
template <class Skin_Traits_>
void
Mixed_complex_builder_3<Skin_Traits_>::construct_3_cells() {
  Sc_Vertex_handle vh[4];
  Sc_Cell_handle ch;

//   Quadratic_surface *surf;
  for (Rt_Finite_cells_iterator cit = T.finite_cells_begin();
       cit != T.finite_cells_end();
       cit++) {
//     surf = new Sphere_surface(
//       traits.s2m_converter(
// 	Sc_Weighted_point(
// 	  orthocenter_obj(
// 	    traits.r2s_converter(cit->vertex(0)->point()),
// 	    traits.r2s_converter(cit->vertex(1)->point()),
// 	    traits.r2s_converter(cit->vertex(2)->point()),
// 	    traits.r2s_converter(cit->vertex(3)->point())),
// 	  orthoweight_obj(
// 	    traits.r2s_converter(cit->vertex(0)->point()),
// 	    traits.r2s_converter(cit->vertex(1)->point()),
// 	    traits.r2s_converter(cit->vertex(2)->point()),
// 	    traits.r2s_converter(cit->vertex(3)->point()))
// 	  )),
//       1-shrink_m, -1);
		
    // construct the tetrahedron:
    //   C[ch], C[Facet(ch,fi)], C[Edge(ch,ei,vi)], C[ch->vertex(vi)]
    Rt_Simplex simplex(cit);
    vh[0] = index_33[cit];
    for (int fi=0; fi<4; fi++) {
      vh[1] = index_23[cit].V[fi];
      if (vh[0] != vh[1]) {
	for (int ei=1; ei<4; ei++) {
	  int index1 = (fi+ei)&3;
	  for (int vi=1; vi<4; vi++) {
	    int index0 = (fi+vi)&3;
	    if (vi != ei) {
	      vh[2] = index_13[cit].V[edge_index[index0][index1]];
	      vh[3] = index_03[cit].V[index0];
	      if ((vh[1] != vh[2]) && (vh[2] != vh[3])) {
		int orient;
								
		if (((4+index1-index0)&3) == 1) {
		  orient = (index1 + (vi==2))&1;
		} else {
		  orient = (index1 + (vi==3))&1;
		}
		ch = add_cell(vh, orient, simplex);
// 		ch->surf = surf;
	      }
	    }
	  }
	}
      }
    }
  }
}

template <class Skin_Traits_>
void 
Mixed_complex_builder_3<Skin_Traits_>::compute_value_on_vertices() {
//   Sc_Cell_handle ch;
//   for (Sc_Finite_vertices_iterator vit=sc.finite_vertices_begin();
//        vit!=sc.finite_vertices_end(); vit++) {
//     ch = vit->cell();
//     if (sc.is_infinite(ch)) {
//       ch = ch->neighbor(ch->index(sc.infinite_vertex()));
//     }
//     vit->iValue = ch->surf->value(traits.s2m_converter(vit->point()));
//   }
	
//   for (Sc_Finite_cells_iterator cit = sc.finite_cells_begin();
//        cit != sc.finite_cells_end();
//        cit ++) {
//     cit->initialise();
//   }
}


// Adds a vertex to the simplicial complex
template <class Skin_Traits_>
typename Mixed_complex_builder_3<Skin_Traits_>::Sc_Vertex_handle
Mixed_complex_builder_3<Skin_Traits_>::add_vertex (
						   Rt_Simplex &sDel, Rt_Simplex &sVor)
{
  Sc_Vertex_handle vh = triangulation_incr_builder.add_vertex();
  vh->point() = get_anchor(sDel, sVor);

  return vh;
}

// Gets a vertex from the simplicial complex based on the anchors
template <class Skin_Traits_>
typename Mixed_complex_builder_3<Skin_Traits_>::Sc_Vertex_handle
Mixed_complex_builder_3<Skin_Traits_>::get_vertex (
						   Rt_Simplex &sDel, Rt_Simplex &sVor)
{
  Rt_Vertex_handle vh;
  Rt_Edge e;
  Rt_Facet f;
  Rt_Cell_handle ch;
  if (sDel.dimension() == 0) {
    vh = sDel;
    switch (sVor.dimension()) {
    case 0:
      assert(sDel == sVor);
      return index_00[vh];
    case 1:
      e=sVor;
      if (e.first->vertex(e.second) == vh) {
	return index_01[vh].V[e.first->vertex(e.third)];
      } else {
	return index_01[vh].V[e.first->vertex(e.second)];
      }
    case 2:
      f=sVor;
      return index_02[f.first].V[f.second][f.first->index(vh)];
    case 3:
      ch=sVor;
      return index_03[ch].V[ch->index(vh)];
    default:
      assert(false);
    }
  } else if (sDel.dimension() == 1) {
    e=sDel;
    switch (sVor.dimension()) {
    case 1:
      assert(sDel == sVor);
      return index_11[e.first->vertex(e.second)].V[e.first->vertex(e.third)];
    case 2:
      f=sVor;
      return index_12[f.first].V[f.second][6 - f.second -
					   f.first->index(e.first->vertex(e.second)) -
					   f.first->index(e.first->vertex(e.third))];
    case 3:
      ch=sVor;
      assert(index_13.is_defined(ch));
      return index_13[ch].V[edge_index
			    [ch->index(e.first->vertex(e.third))]
			    [ch->index(e.first->vertex(e.second))]];
    default:
      assert(false);
    }
  } else if (sDel.dimension() == 2) {
    f=sDel;
    switch (sVor.dimension()) {
    case 2:
      assert(sDel == sVor);
      return index_22[f.first].V[f.second];
    case 3:
      ch=sVor;
      if (f.first == ch) {
	return index_23[ch].V[f.second];
      } else {
	assert(f.first->neighbor(f.second) == ch);
	return index_23[ch].V[f.first->neighbor(f.second)->index(f.first)];
      }
    default:
      assert(false);
    }
  } else if (sDel.dimension() == 3) {
    assert(sDel == sVor);
    ch=sDel;
    return index_33[ch];
  }
  assert(false);
  return Sc_Vertex_handle();
}

// Adds a cell to the simplicial complex
template <class Skin_Traits_>
typename Mixed_complex_builder_3<Skin_Traits_>::Sc_Cell_handle
Mixed_complex_builder_3<Skin_Traits_>::
add_cell(Sc_Vertex_handle vh[], int orient, Rt_Simplex s) {
  assert((orient==0) || (orient==1));
  assert(vh[0] != Sc_Vertex_handle()); assert(vh[1] != Sc_Vertex_handle());
  assert(vh[2] != Sc_Vertex_handle()); assert(vh[3] != Sc_Vertex_handle());
  assert(vh[1] != vh[2]); assert(vh[1] != vh[3]); assert(vh[1] != vh[4]);
  assert(vh[2] != vh[3]); assert(vh[2] != vh[4]); assert(vh[3] != vh[4]);
  Sc_Cell_handle ch;
  if (orient) {
    ch = triangulation_incr_builder.add_cell(vh[0], vh[1], vh[2], vh[3]);
  } else {
    ch = triangulation_incr_builder.add_cell(vh[0], vh[1], vh[3], vh[2]);
  }
  backRef[ch] = s;

  return ch;
}

template <class Skin_Traits_>
typename Mixed_complex_builder_3<Skin_Traits_>::Sc_Point
Mixed_complex_builder_3<Skin_Traits_>::get_anchor(Rt_Simplex &sDel, Rt_Simplex &sVor)
{
  Rt_Vertex_handle vh;
  Rt_Edge           e;
  Rt_Facet          f;
  Rt_Cell_handle   ch;
	
  Sc_Point vfoc, dfoc;
  switch (sVor.dimension()) {
  case 0:
    vh=sVor;
    vfoc = traits.r2s_converter_object()(vh->point());
    break;
  case 1:
    e=sVor;
    vfoc = orthocenter_obj(
			   traits.r2s_converter_object()(e.first->vertex(e.second)->point()),
			   traits.r2s_converter_object()(e.first->vertex(e.third)->point()));
    break;
  case 2:
    f=sVor;
    vfoc = orthocenter_obj(
			   traits.r2s_converter_object()(f.first->vertex((f.second+1)&3)->point()),
			   traits.r2s_converter_object()(f.first->vertex((f.second+2)&3)->point()),
			   traits.r2s_converter_object()(f.first->vertex((f.second+3)&3)->point()));
    break;
  case 3:
    ch=sVor;
    vfoc = orthocenter_obj(
			   traits.r2s_converter_object()(ch->vertex(0)->point()),
			   traits.r2s_converter_object()(ch->vertex(1)->point()),
			   traits.r2s_converter_object()(ch->vertex(2)->point()),
			   traits.r2s_converter_object()(ch->vertex(3)->point()));
    break;
  }
  switch (sDel.dimension()) {
  case 0:
    vh=sDel;
    dfoc = traits.r2s_converter_object()(vh->point());
    break;
  case 1:
    e=sDel;
    dfoc = orthocenter_obj(
			   traits.r2s_converter_object()(e.first->vertex(e.second)->point()),
			   traits.r2s_converter_object()(e.first->vertex(e.third)->point()));
    break;
  case 2:
    f=sDel;
    dfoc = orthocenter_obj(
			   traits.r2s_converter_object()(f.first->vertex((f.second+1)&3)->point()),
			   traits.r2s_converter_object()(f.first->vertex((f.second+2)&3)->point()),
			   traits.r2s_converter_object()(f.first->vertex((f.second+3)&3)->point()));
    break;
  case 3:
    ch=sDel;
    dfoc = orthocenter_obj(
			   traits.r2s_converter_object()(ch->vertex(0)->point()),
			   traits.r2s_converter_object()(ch->vertex(1)->point()),
			   traits.r2s_converter_object()(ch->vertex(2)->point()),
			   traits.r2s_converter_object()(ch->vertex(3)->point()));
    break;
  }
	
  return dfoc + shrink*(vfoc - dfoc);
}

CGAL_END_NAMESPACE

#endif // MIXED_COMPLEX_BUILDER_H

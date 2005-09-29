#ifndef MARCHING_TETRAHEDRA_H
#define MARCHING_TETRAHEDRA_H

#define COARSE_PERTURBATION  1e-05

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Cartesian_converter.h>

CGAL_BEGIN_NAMESPACE 

template <class Simplicial_complex,
	  class HalfedgeDS,
	  class MarchingTetrahedraTraits_3 >
class Mesh_builder : public CGAL::Modifier_base<HalfedgeDS> {
 public:
  typedef Simplicial_complex                    Sc;
  typedef HalfedgeDS                            HDS;
  typedef MarchingTetrahedraTraits_3            MarchingTetrahedraTraits;
 private:
  typedef typename Sc::Vertex_handle            Sc_vertex_handle;
  typedef typename Sc::Edge                     Sc_edge;
  typedef typename Sc::Facet                    Sc_facet;
  typedef typename Sc::Cell_handle              Sc_cell_handle;
  
  typedef typename Sc::Finite_vertices_iterator Sc_finite_vertices_iterator;
  typedef typename Sc::Finite_edges_iterator    Sc_finite_edges_iterator;
  typedef typename Sc::Finite_facets_iterator   Sc_finite_facets_iterator;
  typedef typename Sc::All_cells_iterator       Sc_all_cells_iterator;
  typedef typename Sc::Finite_cells_iterator    Sc_finite_cells_iterator;
  
  typedef typename Sc::Cell_circulator          Sc_cell_circulator;
  typedef typename Sc::Geom_traits              Sc_traits;
  typedef typename Sc_traits::RT                Sc_rt;
  typedef typename Sc_traits::Point_3           Sc_point;
  typedef typename Sc_traits::Triangle_3        Sc_triangle;
  
  typedef typename HDS::Vertex_handle           HDS_vertex_handle;
  typedef typename HDS::Halfedge_handle         HDS_halfedge_handle;
  typedef typename HDS::Face_handle             HDS_face_handle;
  
  typedef typename HDS::Traits                  HDS_traits;
  typedef typename HDS_traits::Point_3          HDS_point;
  typedef typename HDS_traits::RT               HDS_rt;
  
  typedef std::pair<Sc_vertex_handle,Sc_vertex_handle> Vpair;
  typedef CGAL::Bounded_side                    Bounded_side;

 public:
  Mesh_builder(Sc &sc, const MarchingTetrahedraTraits &traits)
    : sc(sc), traits(traits), nVertices(0) {
  }
  
  void operator()( HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    
    B.begin_surface(0,0,0);
    
    // First generate vertices of the mesh:
    Sc_vertex_handle vh0, vh1;
    for (Sc_finite_edges_iterator eit=sc.finite_edges_begin();
	 eit!=sc.finite_edges_end(); eit++) {
      if ((traits.sign(eit->first->vertex(eit->second))==POSITIVE)!=
	  (traits.sign(eit->first->vertex(eit->third))==POSITIVE)) {
	add_vertex(B,eit);
      }
    }
    
    // Then triangulate each cell of the triangulation:
    int in[4], out[4], Nin;
    Vpair vpair;
    for (Sc_finite_cells_iterator cit= sc.finite_cells_begin();
         cit!=sc.finite_cells_end();
         cit++) {
      Nin = 0;
      
      for (int i=0; i<4; i++) { in[i]=-1; out[i]=-1; }

      for (int i=0; i<4; i++) {
        if (traits.sign(cit->vertex(i)) == POSITIVE) {
          in[Nin] = i; Nin++;
        } else {
          out[i-Nin] = i;
        }
      }
      if (Nin==1) {
        if ((in[0]%2) == 1) {
          add_facet(B,
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[0]))],
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[1]))],
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[2]))],
            cit);
	  assert(!B.error());
        } else {
          add_facet(B,
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[0]))],
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[2]))],
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[1]))],
            cit);
	  assert(!B.error());
        }
      } else if (Nin==2) {
        if (((out[0] == 0) && (out[1]==1)) || 
	    ((out[0] == 0) && (out[1]==3)) || 
	    ((out[0] == 1) && (out[1]==2)) || 
	    ((out[0] == 2) && (out[1]==3)))
	{
          add_facet(B,
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[0]))],
            vertices[makePair(cit->vertex(in[1]), cit->vertex(out[0]))],
            vertices[makePair(cit->vertex(in[1]), cit->vertex(out[1]))],
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[1]))],
            cit);
	  assert(!B.error());
        } else {
          add_facet(B,
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[0]))],
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[1]))],
            vertices[makePair(cit->vertex(in[1]), cit->vertex(out[1]))],
            vertices[makePair(cit->vertex(in[1]), cit->vertex(out[0]))],
            cit);
	  assert(!B.error());
        }
      } else if (Nin==3) {
        if ((out[0]%2) == 1) {
          add_facet(B,
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[0]))],
            vertices[makePair(cit->vertex(in[2]), cit->vertex(out[0]))],
            vertices[makePair(cit->vertex(in[1]), cit->vertex(out[0]))],
            cit);
	  assert(!B.error());
        } else {
          add_facet(B,
            vertices[makePair(cit->vertex(in[0]), cit->vertex(out[0]))],
            vertices[makePair(cit->vertex(in[1]), cit->vertex(out[0]))],
            vertices[makePair(cit->vertex(in[2]), cit->vertex(out[0]))],
            cit);
	  assert(!B.error());
        }
       }
      assert(!B.error());
    }

    B.end_surface();
  }
private:
  Sc &sc;
  const MarchingTetrahedraTraits &traits;
  std::map< Vpair, int > vertices;
  int nVertices;

  Vpair makePair(Sc_vertex_handle vh0, Sc_vertex_handle vh1) {
    if (vh0 < vh1) return Vpair(vh0,vh1); else return Vpair(vh1,vh0);
  }

  void add_vertex(CGAL::Polyhedron_incremental_builder_3<HDS> &B,
    Sc_finite_edges_iterator const &e) {
    Vpair vpair = makePair(
      e->first->vertex(e->second),
      e->first->vertex(e->third));
    Sc_cell_circulator ccir = sc.incident_cells(*e);
    while (sc.is_infinite(ccir)) ccir ++;
    assert(!sc.is_infinite(ccir));

    B.add_vertex(traits.intersection(Sc_edge(ccir, 
      ccir->index(e->first->vertex(e->second)),
      ccir->index(e->first->vertex(e->third)))));
    vertices[vpair] = nVertices;
    nVertices ++;
  }
  
  // Orientation is right
  void add_facet(CGAL::Polyhedron_incremental_builder_3<HDS> &B,
    int v0, int v1, int v2,
    Sc_cell_handle ch) {
    assert((v0!=v1) && (v0!=v2) && (v1!=v2));
    HDS_face_handle f = B.begin_facet();
    B.add_vertex_to_facet( v0 );
    B.add_vertex_to_facet( v1 );
    B.add_vertex_to_facet( v2 );
    B.end_facet();
  }

  // Add quadruple or two triangles (orientation is right)
  void add_facet(CGAL::Polyhedron_incremental_builder_3<HDS> &B,
    int v0, int v1, int v2, int v3,
    Sc_cell_handle ch) {
//     assert((v0!=v1) && (v0!=v2) && (v0!=v3));
//     assert((v1!=v2) && (v1!=v3) && (v2!=v3));
//     HDS_Face_handle f = B.begin_facet();
//     B.add_vertex_to_facet( v0 );
//     B.add_vertex_to_facet( v1 );
//     B.add_vertex_to_facet( v2 );
//     B.add_vertex_to_facet( v3 );
//     B.end_facet();
// OR:
    add_facet(B, v0, v1, v3, ch);
    add_facet(B, v1, v2, v3, ch);
  }
};

template <class Triangulation_3,
	  class Polyhedron_3,
	  class MarchingTetrahedraTraits_3 >
void marching_tetrahedra_3(
  Triangulation_3 &triangulation,
  Polyhedron_3 &polyhedron,
  const MarchingTetrahedraTraits_3 &marching_traits) {
  
  typedef typename Polyhedron_3::HalfedgeDS                    HDS;
  typedef MarchingTetrahedraTraits_3                           Marching_traits;
  typedef Mesh_builder<Triangulation_3,HDS,Marching_traits>    Mesh_builder;
  
  Mesh_builder builder(triangulation, marching_traits);
  polyhedron.delegate(builder);
}

CGAL_END_NAMESPACE 

#endif // MARCHING_TETRAHEDRA_H

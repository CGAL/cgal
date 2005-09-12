#ifndef MARCHING_TETRAHEDRA_H
#define MARCHING_TETRAHEDRA_H

#define COARSE_PERTURBATION  1e-05

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Cartesian_converter.h>

CGAL_BEGIN_NAMESPACE 

template <class T_, class HDS_, class Converter_ >
class Skin_surface_extractor {
public:
  typedef T_         T;
  typedef HDS_       HDS;
  typedef Converter_ Converter;

  typedef typename T::Vertex_handle            Sc_vertex_handle;
  typedef typename T::Edge                     Sc_edge;
  typedef typename T::Facet                    Sc_facet;
  typedef typename T::Cell_handle              Sc_cell_handle;
  typedef typename T::Geom_traits::Point_3     Sc_point;

  typedef typename HDS::Traits                 Mesh_K;
  typedef typename Mesh_K::RT                  Mesh_rt;
  typedef typename Mesh_K::Point_3             Mesh_point;

  Skin_surface_extractor(Mesh_rt iso_value=0) : iso_value(iso_value) {
  }
  
  Sign sign(Sc_vertex_handle const vh) {
    return CGAL::sign(
      vh->cell()->surf->value(converter(vh->point())) - iso_value);
  }
  Mesh_point intersection(Sc_edge const& e) {
    // Precondition: e.first is not an infinite cell: they have not surface set
    return e.first->surf->to_surface(
      converter(e.first->vertex(e.second)->point()),
      converter(e.first->vertex(e.third)->point()));
  }
  
  Converter converter;
  Mesh_rt iso_value;
};


template <class Simplicial_complex, class HalfedgeDS, class MeshExtractor >
class Mesh_builder : public CGAL::Modifier_base<HalfedgeDS> {
 public:
  typedef Simplicial_complex Sc;
  typedef HalfedgeDS         HDS;
  typedef MeshExtractor      Extractor;
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
  Mesh_builder(Sc &sc, Extractor &extr)
    : sc(sc), extr(extr), nVertices(0) {
  }
  
  void operator()( HDS& hds) {
    CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
    
    B.begin_surface(0,0,0);
    
    // First generate vertices of the mesh:
    Sc_vertex_handle vh0, vh1;
    for (Sc_finite_edges_iterator eit=sc.finite_edges_begin();
	 eit!=sc.finite_edges_end(); eit++) {
      if ((extr.sign(eit->first->vertex(eit->second))==POSITIVE)!=
	  (extr.sign(eit->first->vertex(eit->third))==POSITIVE)) {
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
        if (extr.sign(cit->vertex(i)) == POSITIVE) {
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
  Extractor &extr;
  std::map< Vpair, int > vertices;
  int nVertices;
//CGAL::Cartesian_converter<Sc_traits, HDS_traits, CGAL::To_double<Sc_RT> > sc2m_converter;
//CGAL::Cartesian_converter<HDS_traits, Sc_traits> m2sc_converter;

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

    B.add_vertex(extr.intersection(Sc_edge(ccir, 
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

template < 
  class Triangulation_3, class Mesh_3, 
  class S2M_converter = Cartesian_converter<
      typename Triangulation_3::Geom_traits,
      typename Mesh_3::Traits>,
  class MarchingSurfaceExtractor =
    Skin_surface_extractor<
      Triangulation_3, typename Mesh_3::HalfedgeDS,S2M_converter > >
class Marching_tetrahedra_3 {
public:
  typedef Triangulation_3                     Triang;
  typedef Mesh_3                              Mesh;
  typedef typename Mesh::HalfedgeDS           HDS;
  typedef S2M_converter                       Converter;
  typedef MarchingSurfaceExtractor            Extractor;
  typedef typename Mesh::Traits::RT           Mesh_rt;
  typedef Mesh_builder<Triang,HDS,Extractor>  Mesh_builder;
  
  Marching_tetrahedra_3() {}
  
  void operator()(Triang &t, Mesh &m, Mesh_rt iso=0) {
    Extractor extr(iso);
    Mesh_builder builder(t, extr);
    m.delegate(builder);
  }
  void operator()(Triang &t, Mesh &m, Extractor &extr, Mesh_rt iso=0) {
    Mesh_builder builder(t, extr);
    m.delegate(builder);
  }
};

CGAL_END_NAMESPACE 

#endif // MARCHING_TETRAHEDRA_H

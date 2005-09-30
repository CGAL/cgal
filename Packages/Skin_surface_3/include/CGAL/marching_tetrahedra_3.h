#ifndef MARCHING_TETRAHEDRA_H
#define MARCHING_TETRAHEDRA_H

#define COARSE_PERTURBATION  1e-05

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Marching_tetrahedra_observer_default_3.h>

CGAL_BEGIN_NAMESPACE 

template <class Triangulation_3,
	  class HalfedgeDS,
	  class MarchingTetrahedraTraits_3,
	  class MarchingTetrahedraObserver_3 >
class Marching_tetrahedra_builder : public Modifier_base<HalfedgeDS> {
 public:
  typedef Triangulation_3                       Triang;
  typedef HalfedgeDS                            HDS;
  typedef MarchingTetrahedraTraits_3            Traits;
  typedef MarchingTetrahedraObserver_3          Observer;
 private:
  typedef typename Triang::Vertex_handle            Triang_vertex_handle;
  typedef typename Triang::Edge                     Triang_edge;
  typedef typename Triang::Facet                    Triang_facet;
  typedef typename Triang::Cell_handle              Triang_cell_handle;
  
  typedef typename Triang::Finite_vertices_iterator Triang_finite_vertices_iterator;
  typedef typename Triang::Finite_edges_iterator    Triang_finite_edges_iterator;
  typedef typename Triang::Finite_facets_iterator   Triang_finite_facets_iterator;
  typedef typename Triang::All_cells_iterator       Triang_all_cells_iterator;
  typedef typename Triang::Finite_cells_iterator    Triang_finite_cells_iterator;
  
  typedef typename Triang::Cell_circulator          Triang_cell_circulator;
  typedef typename Triang::Geom_traits              Triang_traits;
  typedef typename Triang_traits::RT                Triang_rt;
  typedef typename Triang_traits::Point_3           Triang_point;
  typedef typename Triang_traits::Triangle_3        Triang_triangle;
  
  typedef typename HDS::Vertex_handle           HDS_vertex_handle;
  typedef typename HDS::Halfedge_handle         HDS_halfedge_handle;
  typedef typename HDS::Face_handle             HDS_face_handle;
  
  typedef typename HDS::Traits                  HDS_traits;
  typedef typename HDS_traits::Point_3          HDS_point;
  typedef typename HDS_traits::RT               HDS_rt;
  
  typedef std::pair<Triang_vertex_handle,Triang_vertex_handle> Vpair;

 public:
  Marching_tetrahedra_builder(const Triang &triang, const Traits &traits, Observer &observer)
    : triang(triang), traits(traits), observer(observer), nVertices(0) {
  }
  
  void operator()( HDS& hds) {
    Polyhedron_incremental_builder_3<HDS> B( hds, true);
    
    B.begin_surface(0,0,0);
    
    // First generate vertices of the mesh:
    Triang_vertex_handle vh0, vh1;
    for (Triang_finite_edges_iterator eit=triang.finite_edges_begin();
	 eit!=triang.finite_edges_end(); eit++) {
      if ((traits.sign(eit->first->vertex(eit->second))==POSITIVE)!=
	  (traits.sign(eit->first->vertex(eit->third))==POSITIVE)) {
	add_vertex(B,eit);
      }
    }
    
    // Then triangulate each cell of the triangulation:
    int in[4], out[4], Nin;
    Vpair vpair;
    for (Triang_finite_cells_iterator cit= triang.finite_cells_begin();
         cit!=triang.finite_cells_end();
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
  const Triang &triang;
  const Traits &traits;
  Observer &observer;
  std::map< Vpair, int > vertices;
  int nVertices;

  Vpair makePair(Triang_vertex_handle vh0, Triang_vertex_handle vh1) {
    if (vh0 < vh1) return Vpair(vh0,vh1); else return Vpair(vh1,vh0);
  }

  void add_vertex(Polyhedron_incremental_builder_3<HDS> &B,
    Triang_finite_edges_iterator const &e) {
    Vpair vpair = makePair(
      e->first->vertex(e->second),
      e->first->vertex(e->third));
    Triang_cell_circulator ccir = triang.incident_cells(*e);
    while (triang.is_infinite(ccir)) ccir ++;
    assert(!triang.is_infinite(ccir));

    B.add_vertex(traits.intersection(Triang_edge(ccir, 
      ccir->index(e->first->vertex(e->second)),
      ccir->index(e->first->vertex(e->third)))));
    vertices[vpair] = nVertices;
    nVertices ++;
  }
  
  // Orientation is right
  void add_facet(Polyhedron_incremental_builder_3<HDS> &B,
    int v0, int v1, int v2,
    Triang_cell_handle ch) {
    assert((v0!=v1) && (v0!=v2) && (v1!=v2));
    HDS_face_handle f = B.begin_facet();
    B.add_vertex_to_facet( v0 );
    B.add_vertex_to_facet( v1 );
    B.add_vertex_to_facet( v2 );
    B.end_facet();
  }

  // Add quadruple or two triangles (orientation is right)
  void add_facet(Polyhedron_incremental_builder_3<HDS> &B,
    int v0, int v1, int v2, int v3,
    Triang_cell_handle ch) {
    add_facet(B, v0, v1, v3, ch);
    add_facet(B, v1, v2, v3, ch);
  }
};

template <class Triangulation_3,
	  class Polyhedron_3,
	  class MarchingTetrahedraTraits_3,
	  class MarchingTetrahedraObserver_3>
void marching_tetrahedra_3(
  const Triangulation_3 &triangulation,
  Polyhedron_3 &polyhedron,
  const MarchingTetrahedraTraits_3 &traits,
  MarchingTetrahedraObserver_3 &observer) {
  
  typedef typename Polyhedron_3::HalfedgeDS                   HDS;
  typedef MarchingTetrahedraTraits_3                          Traits;
  typedef MarchingTetrahedraObserver_3                        Observer;
  typedef Marching_tetrahedra_builder<Triangulation_3,HDS, Traits, Observer>
                                                              Builder;
  
  Builder builder(triangulation, traits, observer);
  polyhedron.delegate(builder);
}

template <class Triangulation_3,
	  class Polyhedron_3,
	  class MarchingTetrahedraTraits_3>
void marching_tetrahedra_3(
  const Triangulation_3 &triangulation,
  Polyhedron_3 &polyhedron,
  const MarchingTetrahedraTraits_3 &traits) {
  
  typedef typename Polyhedron_3::HalfedgeDS                    HDS;
  typedef typename Triangulation_3::Cell_handle                Cell_handle;
  typedef MarchingTetrahedraTraits_3                           Traits;
  typedef Marching_tetrahedra_observer_default_3<Cell_handle, Polyhedron_3>
                                                               Observer; 
  typedef Marching_tetrahedra_builder<Triangulation_3,HDS,Traits,Observer>
                                                               Builder;

  Observer observer;
  Builder builder(triangulation, traits, observer);
  polyhedron.delegate(builder);
}

CGAL_END_NAMESPACE 

#endif // MARCHING_TETRAHEDRA_H

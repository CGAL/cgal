#ifndef SKIN_SURFACE_SQRT3_3_H
#define SKIN_SURFACE_SQRT3_3_H

#include <CGAL/Skin_surface_refinement_traits_3.h>

CGAL_BEGIN_NAMESPACE

// This code is based on the Polyhedron tutorial 

template <class Polyhedron_3,
	  class SkinSurfaceRefinementTraits_3>
class Skin_surface_sqrt3
{
  typedef Polyhedron_3                                        Polyhedron;
  typedef SkinSurfaceRefinementTraits_3                       Refinement_traits;
  typedef typename Refinement_traits::Triangulation           Triangulation;
  
  typedef typename Polyhedron::Traits                         Kernel;
  typedef typename Kernel::Point_3                            Point;
  typedef typename Kernel::Vector_3                           Vector;
	
  typedef typename Polyhedron::Vertex                         Vertex;
  typedef typename Polyhedron::Vertex_iterator                Vertex_iterator;
  typedef typename Polyhedron::Edge_iterator                  Edge_iterator;
  typedef typename Polyhedron::Halfedge_handle                Halfedge_handle;
  typedef typename Polyhedron::Halfedge_iterator              Halfedge_iterator;
  typedef typename Polyhedron::Halfedge_around_vertex_const_circulator  
                                                              HV_circulator;
  typedef typename Polyhedron::Halfedge_around_facet_circulator
                                                              HF_circulator;
  typedef typename Polyhedron::Facet                          Facet;
  typedef typename Polyhedron::Facet_handle                   Facet_handle;
  typedef typename Polyhedron::Facet_iterator                 Facet_iterator;
  typedef typename Kernel::FT                                 FT;

public:
  Skin_surface_sqrt3(Polyhedron &P, Refinement_traits &traits) 
    : P(P), traits(traits) {}

  //*********************************************
  // Subdivision
  //*********************************************
  int subdivide(int iter)
  {
    // check for valid polygon mesh
    if(P.size_of_facets() == 0)
      return false;

    while (iter > 0) {
      // normalize border
      P.normalize_border();

      if (iter >1) {
	subdivide_twice();
	iter = iter - 2;
      } else {
	subdivide_once();
	iter--;
      }
    }
    return true;
  }

private:
	
  //*********************************************
  // Subdivide
  //*********************************************
  void subdivide_once()
  {

    // We use that new vertices/halfedges/facets are appended at the end.
    Vertex_iterator last_v = P.vertices_end();
    -- last_v;  // the last of the old vertices
    Edge_iterator last_e = P.edges_end();
    -- last_e;  // the last of the old edges
    Facet_iterator last_f = P.facets_end();
    -- last_f;  // the last of the old facets

    // split edges
    Edge_iterator  e = P.edges_begin ();
    do {
      split_halfedge(e);
    } while ( e++ != last_e);


    Vertex_iterator v = P.vertices_begin();
    do {
      Halfedge_handle h_cir, h_start;
      h_cir = h_start = v->halfedge () ;
      do {
        P.split_facet (h_cir->prev(), h_cir->next());
        h_cir = h_cir->next()->opposite();
      } while (h_cir != h_start);
    } while (v++ != last_v);

    v = ++last_v; // First new vertex
    last_v = P.vertices_end();
    -- last_v;  // the last of the old vertices
    do {
      to_surface(v);
    } while (v++ != last_v);

    CGAL_postcondition( P.is_valid());
  }

  //*********************************************
  // Subdivide
  //*********************************************
  void subdivide_twice()
  {

    // We use that new vertices/halfedges/facets are appended at the end.
    Vertex_iterator last_v = P.vertices_end();
    -- last_v;  // the last of the old vertices
    Edge_iterator last_e = P.edges_end();
    -- last_e;  // the last of the old edges
    Facet_iterator last_f = P.facets_end();
    -- last_f;  // the last of the old facets

    // split edges
    Edge_iterator  e = P.edges_begin ();
    do {
      trisect_halfedge(P, e);
    } while ( e++ != last_e);

    Vertex_iterator v = P.vertices_begin();
    do {
      Halfedge_handle h_cir, h_start;
      h_cir = h_start = v->halfedge () ;
      do {
        P.split_facet (h_cir->prev(), h_cir->next());
        h_cir = h_cir->next()->opposite();
      } while (h_cir != h_start);
    } while (v++ != last_v);

    Facet_iterator f = P.facets_begin();
    do {
      create_center_vertex(P, f);
    } while (f++ != last_f);

    v = ++last_v; // First new vertex
    last_v = P.vertices_end();
    -- last_v;  // the last of the old vertices
    do {
      to_surface(v);
    } while (v++ != last_v);

    CGAL_postcondition( P.is_valid());
  }
  
  //*********************************************
  // Split halfedge
  //*********************************************
  void split_halfedge(Halfedge_handle e)
  {
    // Create two new vertices on e.
    Point p_new = e->vertex()->point();
    e = e->prev();
    p_new = p_new + .5*(e->vertex()->point()-p_new);

    P.split_vertex( e, e->next()->opposite());
    e->next()->vertex()->point() = p_new;
  }

  //*********************************************
  // Trisect halfedge
  //*********************************************
  void trisect_halfedge(Polyhedron& P,
			Halfedge_handle e)
  {
    // Create two new vertices on e.
    Point p2 = e->vertex()->point();
    e = e->prev();
    Point p1 = e->vertex()->point();
    Vector v = p2 - p1;

    P.split_vertex( e, e->next()->opposite());
    e->next()->vertex()->point() = p1 + .6666 * v;
    P.split_vertex( e, e->next()->opposite());
    e->next()->vertex()->point() = p1 + .3333 * v;
  }


  //*********************************************
  // Create center vertex
  //*********************************************
  void create_center_vertex(Polyhedron& P,
			    Facet_iterator f)
  {
    Vector vec( 0.0, 0.0, 0.0);
    std::size_t order = 0;
    HF_circulator h = f->facet_begin();
    do {
      vec = vec + ( h->vertex()->point() - CGAL::ORIGIN);
      ++ order;
    } while ( ++h != f->facet_begin());
    CGAL_assertion( order >= 3); // guaranteed by definition of Polyhedron
    Point center =  CGAL::ORIGIN + (vec / (FT)order);
    Halfedge_handle new_center = P.create_center_vertex( f->halfedge());
    new_center->vertex()->point() = center;
  }

  void to_surface(Vertex_iterator v) {
    v->point() = traits.to_surface_along_transversal_segment(v->point());
  }
  
  Polyhedron &P;
  Refinement_traits &traits;
};

template <class Polyhedron_3,
	  class SkinSurfaceRefinementTraits_3>
void skin_surface_sqrt3(
          Polyhedron_3 &p, 
          SkinSurfaceRefinementTraits_3 &traits,
          int nSubdiv = 2) {
            
  typedef Skin_surface_sqrt3<Polyhedron_3, SkinSurfaceRefinementTraits_3> Sqrt3;
  Sqrt3 subdivider(p, traits);
            
  subdivider.subdivide(nSubdiv);
}

CGAL_END_NAMESPACE

#endif // SKIN_SURFACE_SQRT3_3_H

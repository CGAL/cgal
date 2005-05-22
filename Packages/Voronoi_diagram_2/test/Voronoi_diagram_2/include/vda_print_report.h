#ifndef VDA_PRINT_REPORT_H
#define VDA_PRINT_REPORT_H 1

#include <CGAL/basic.h>
#include "vda_aux.h"
#include <iostream>

template<class VDA, class Projector, class Dual_primal_projector,
	 class Stream>
void print_report(const VDA& vda, const Projector& project,
		  const Dual_primal_projector& dp_project,
		  Stream& os = std::cout)
{
  std::cout << std::endl;
  std::cout << "is Delaunay graph valid? "
	    << (vda.dual().is_valid() ? "yes" : "no") << std::endl;
  std::cout << std::endl;


  if ( vda.dual().number_of_vertices() == 1 ) {
    std::cout << "The dual graph (Delaunay graph) has only 1 site."
	      << std::endl;
    std::cout << "This implies that the dual graph is 0-dimensional."
	      << std::endl;
    std::cout << "Cannot view the dual graph as an arrangement."
	      << std::endl;
    return;
  }

  typename VDA::Edge_iterator eit;
  for (eit = vda.edges_begin(); eit != vda.edges_end(); ++eit) {
    print_halfedge(vda, *eit, project, os);
  }

  os << std::endl << std::endl << std::endl;


  typename VDA::Halfedge_iterator heit;
  for (heit = vda.halfedges_begin(); heit != vda.halfedges_end(); ++heit) {
    print_halfedge(vda, typename VDA::Halfedge_handle(heit), project, os);
  }

  os << std::endl << std::endl;
  os << "=====================" << std::endl;
  os << std::endl << std::endl;

  typename VDA::Face_iterator fit = vda.faces_begin();
  for (; fit != vda.faces_end(); ++fit) {
    os << "Face of: ";
    print_dual_vertex(vda, fit->dual_vertex(), project, os);

    typename VDA::Ccb_halfedge_circulator ccb_start;
    // the two lines below are equivalent since Face_iterator is
    // convertible to Face_handle
    //    ccb_start = vda.ccb_halfedges(typename VDA::Face_handle(fit));
    ccb_start = vda.ccb_halfedges(fit);

    typename VDA::Ccb_halfedge_circulator ccb = ccb_start;

    
    os << std::endl;
    os << "TESTING INCREMENT OPERATORS" << std::endl;
    os << std::endl;

    print_halfedge(vda, *ccb, project, os);
    do {
      ++ccb;
      print_halfedge(vda, *ccb, project, os);
    } while ( ccb_start != ccb );
    os << std::endl << std::endl;

    print_halfedge(vda, *ccb, project, os);
    do {
      ccb++;
      print_halfedge(vda, *ccb, project, os);
    } while ( ccb_start != ccb );
    os << std::endl << std::endl;

    os << std::endl;
    os << "TESTING DECREMENT OPERATORS" << std::endl;
    os << std::endl;

    print_halfedge(vda, *ccb, project, os);
    do {
      --ccb;
      print_halfedge(vda, *ccb, project, os);
    } while ( ccb_start != ccb );
    os << std::endl << std::endl;

    print_halfedge(vda, *ccb, project, os);
    do {
      ccb--;
      print_halfedge(vda, *ccb, project, os);
    } while ( ccb_start != ccb );
    os << std::endl << std::endl;
  }

  os << std::endl << std::endl;
  os << "=====================" << std::endl;
  os << std::endl << std::endl;

  int n_all = 0, n_empty = 0, n_vert = 0;
  for (fit = vda.faces_begin(); fit != vda.faces_end(); ++fit) {
    CGAL_assertion( !vda.has_empty_Voronoi_cell_interior(fit->dual_vertex()) );
    n_all++;
  }

  typename VDA::Dual_vertices_iterator vit;
  for ( vit = vda.dual().finite_vertices_begin();
	vit != vda.dual().finite_vertices_end(); ++vit) {
    n_vert++;
    typename VDA::Dual_vertex_handle v(vit);
    if ( vda.has_empty_Voronoi_cell_interior(v) ) {
      n_empty++;
    }    
  }
  
  {
    unsigned int n_dual_edges = 0;
    unsigned int n_dual_fin_faces = 0;
    unsigned int n_dual_faces = 0;
    unsigned int n_edges = 0;
    unsigned int n_edge_inf = 0;
    unsigned int n_edge_fin = 0;
    unsigned int n_edge_degen = 0;
    unsigned int n_edge_non_degen = 0;
    unsigned int n_hedges = 0;
    unsigned int n_halfedges = 0;
    unsigned int n_faces = 0;
    unsigned int n_empty_faces = 0;
    unsigned int n_vertices = 0;
    unsigned int sum_deg = 0;
    unsigned int n_unbounded_faces = 0;
    unsigned int n_unbounded_faces2 = 0;

    // computing statistics in the dual (Delaunay graph)
    typename VDA::Dual_graph::All_edges_iterator deit;
    for (deit = vda.dual().all_edges_begin();
	 deit != vda.dual().all_edges_end(); ++deit) {
      n_dual_edges++;
      if ( vda.is_degenerate_edge(deit) ) {
	os << "degenerate edge: " << std::flush;
	print_dual_edge(vda, *deit, project, os);
	n_edge_degen++;
      } else {
	n_edge_non_degen++;
      }
    }
    os << std::endl;

    typename VDA::Dual_faces_iterator dfit;
    for (dfit = vda.dual().finite_faces_begin();
	 dfit != vda.dual().finite_faces_end(); ++dfit) {
      n_dual_fin_faces++;
    }
    os << std::endl;

    typename VDA::Dual_graph::All_faces_iterator dafit;
    for (dafit = vda.dual().all_faces_begin();
	 dafit != vda.dual().all_faces_end(); ++dafit) {
      n_dual_faces++;
    }
    os << std::endl;

    typename VDA::Dual_vertices_iterator dvit;
    for ( dvit = vda.dual().finite_vertices_begin();
	  dvit != vda.dual().finite_vertices_end(); ++dvit) {
      if ( vda.has_empty_Voronoi_cell_interior(dvit) ) {
	n_empty_faces++;
      }
    }

    // computing statistics in the primal (Voronoi diagram)
    typename VDA::Edge_iterator eit;
    for (eit = vda.edges_begin(); eit != vda.edges_end(); ++eit) {
      typename VDA::Dual_edge e = eit->dual_edge();
      if ( vda.dual().is_infinite(e) ) {
	n_edge_inf++;
      } else {
	n_edge_fin++;
      }
      n_edges++;
    }

    typename VDA::Halfedge_iterator heit;
    for (heit = vda.halfedges_begin(); heit != vda.halfedges_end(); ++heit) {
      typename VDA::Halfedge_handle hh(heit);
      CGAL_assertion( heit->opposite()->opposite() == hh );
      n_hedges++;
    }

    typename VDA::Face_iterator fit;
    for (fit = vda.faces_begin(); fit != vda.faces_end(); ++fit) {
      n_faces++;
    }

    typename VDA::Vertex_iterator vit;
    for (vit = vda.vertices_begin(); vit != vda.vertices_end(); ++vit) {
      typename VDA::Halfedge_around_vertex_circulator vc;
      // the two lines below are equivalent since Vertex_iterator is
      // convertible to Vertex_handle
      //      vc = vda.incident_halfedges(typename VDA::Vertex_handle(vit));
      vc = vda.incident_halfedges(vit);

      typename VDA::Halfedge_around_vertex_circulator vc_start = vc;
      typename VDA::size_type deg = 0;
      do {
	deg++;
	vc++;
      } while ( vc != vc_start );

      typename VDA::Dual_face_handle f = vit->dual_face();

      os << "vertex (degree = " << deg << "): " << std::flush;

      typename VDA::Dual_face_handle fvalid =
	CGAL_VORONOI_DIAGRAM_2_NS::Find_valid_vertex<VDA>()(&vda,f);
      os << dp_project(vda, fvalid) << std::endl;

      n_vertices++;
    }
    os << std::endl;

    // testing circulators
    for (typename VDA::Vertex_iterator vit = vda.vertices_begin();
	 vit != vda.vertices_end(); ++vit) {
      typename VDA::Halfedge_around_vertex_circulator hc, hc_start;
      hc = vit->incident_halfedges();
      hc_start = hc;
      do {
	hc++;
	CGAL_assertion( vit->is_incident_edge(*hc) );
	CGAL_assertion( vit->is_incident_face((*hc)->face()) );
	CGAL_assertion( vit->is_incident_face((*hc)->opposite()->face()) );
      } while ( hc != hc_start );
      sum_deg += vit->degree();
    }

    for (typename VDA::Face_iterator fit = vda.faces_begin();
	 fit != vda.faces_end(); ++fit) {
      typename VDA::Halfedge_handle he = fit->halfedge();
      typename VDA::Ccb_halfedge_circulator hc, hc_start;
      hc = he->ccb();
      hc_start = hc;
      do {
	hc++;
	n_halfedges++;
	CGAL_assertion( fit->is_halfedge_on_outer_ccb(*hc) );
      } while ( hc != hc_start );
    }

    // computing number of unbounded faces
    typename VDA::Unbounded_faces_iterator ufit;
    for (ufit = vda.unbounded_faces_begin();
	 ufit != vda.unbounded_faces_end(); ++ufit) {
      n_unbounded_faces++;
    }

    if ( vda.unbounded_faces_begin() != vda.unbounded_faces_end() ) {
      for (ufit = --vda.unbounded_faces_end();
	   ufit != vda.unbounded_faces_begin(); --ufit) {
	n_unbounded_faces2++;
      }
      n_unbounded_faces2++;
    }

    std::cout << "STATISTICS:" << std::endl;
    std::cout << "-----------" << std::endl;
    std::cout << "# of Voronoi edges: "	<< n_edges << std::endl;
    std::cout << "# of finite Voronoi edges: "
	      << n_edge_fin << std::endl;
    std::cout << "# of infinite Voronoi edges: " << n_edge_inf << std::endl;
    std::cout << std::endl;
    std::cout << "# of Voronoi halfedges: " << n_hedges << std::endl;
    std::cout << "# of Voronoi halfedges (2nd count): "
	      << n_halfedges << std::endl;
    std::cout << "# of Voronoi halfedges (3rd count): "
	      << vda.number_of_halfedges() << std::endl;
    std::cout << "# of Voronoi faces: " << n_faces << std::endl;
    std::cout << "# of Voronoi faces (2nd count): "
	      << vda.number_of_faces() << std::endl;
    std::cout << "# of unbounded faces: " << n_unbounded_faces
	      << std::endl;
    std::cout << "# of unbounded faces (2nd count): " << n_unbounded_faces2
	      << std::endl;
    //    std::cout << "# of unbounded faces (3rd count): "
    //	      << vda.dual().degree( vda.dual().infinite_vertex() )
    //	      << std::endl;
    std::cout << "# of Voronoi vertices: " << n_vertices << std::endl;
    std::cout << "# of Voronoi vertices (2nd count): "
	      << vda.number_of_vertices() << std::endl;
    std::cout << std::endl;

    std::cout << "sum of degrees of Voronoi vertices: "
	      << sum_deg << std::endl;
    std::cout << std::endl;

    std::cout << "# of Delaunay edges: " << n_dual_edges << std::endl;
    std::cout << "# of Voronoi edges in the dual: "
	      << n_dual_edges << std::endl;
    std::cout << "# of non-degenerate Voronoi edges in the dual: "
	      << n_edge_non_degen << std::endl;
    std::cout << "# of degenerate Voronoi edges in the dual: "
	      << n_edge_degen << std::endl;
    std::cout << std::endl;
    std::cout << "# of Delaunay faces (triangles): "
	      << n_dual_faces << std::endl;
    std::cout << "# of finite Delaunay faces (triangles): "
	      << n_dual_fin_faces << std::endl;
    std::cout << "# of Voronoi faces in the dual with empty interior: "
	      << n_empty_faces << std::endl;
    std::cout << std::endl;  
  }


  // now we instantiate the iterators for holes in each one of the
  // faces...
  unsigned int n_holes(0);
  for (typename VDA::Face_iterator fit = vda.faces_begin();
       fit != vda.faces_end(); ++fit) {
    for (typename VDA::Face::Holes_iterator hit = fit->holes_begin();
	 hit != fit->holes_end(); ++hit) {
      typename VDA::Halfedge_handle he = *hit;
      typename VDA::Halfedge_handle he_opp = (*hit)->opposite();
      n_holes++;
    }
  }
  std::cout << "# of holes inside the Voronoi faces: "
	      << n_holes << std::endl;
  std::cout << std::endl;  

#if 0
  std::cout << "# of dual vertices: " << n_vert << std::endl;
  std::cout << "# of Voronoi cells: " << n_all << std::endl;
  std::cout << "# of Voronoi cells with empty interior: "
	    << n_empty << std::endl;
#endif

  std::cout << "is Voronoi diagram valid? "
	    << (vda.is_valid() ? "yes" : "no") << std::endl;
}


#endif // VDA_PRINT_REPORT_H

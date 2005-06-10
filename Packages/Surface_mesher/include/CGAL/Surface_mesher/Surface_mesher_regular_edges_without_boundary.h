#ifndef CGAL_SURFACE_MESHER_REGULAR_EDGES_WITHOUT_BOUNDARY_H
#define CGAL_SURFACE_MESHER_REGULAR_EDGES_WITHOUT_BOUNDARY_H

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3_surface_mesh.h>
#include <CGAL/Surface_mesher/Surface_mesher.h>
#include <CGAL/Surface_mesher/Surface_mesher_regular_edges.h>

namespace CGAL {

  namespace Surface_mesher {

  template < typename Tr, 
    typename Surface,
    typename Criteria,
    typename SMREB = Surface_mesher_regular_edges_base <Tr, Surface, Criteria> >
  class Surface_mesher_regular_edges_without_boundary_base
    : public SMREB {
      
    public:
      typedef Complex_2_in_triangulation_3_surface_mesh<Tr> C2t3;
      typedef typename Tr::Point Point;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;
      typedef typename C2t3::Face_type Face_type;

    public:
      Surface_mesher_regular_edges_without_boundary_base (Tr& t, 
						      C2t3& co, 
						      Surface& s, 
						      Criteria& c)
      :SMREB(t, co, s, c){}

    // Initialization function
    void scan_triangulation_impl() {
      SMREB::scan_triangulation_impl(false);
    }

    void after_insertion_impl(const Vertex_handle v) {
      SMREB::after_insertion_impl(v, false);
    }

  };  // end Surface_mesher_regular_edges

  template <typename Tr,
    typename Surface,
    typename Criteria,
    typename Base = Surface_mesher_regular_edges_without_boundary_base<
      Tr, 
      Surface, 
      Criteria> >
  class Surface_mesher_regular_edges_without_boundary : 
      public Base, 
      public Mesher_level <
    Tr,
    Surface_mesher_regular_edges_without_boundary<Tr, Surface, Criteria, Base>,
    typename Tr::Facet,
    Null_mesher_level,
    Triangulation_mesher_level_traits_3<Tr>
  >
  {
  public:
    typedef Surface_mesher_regular_edges_without_boundary<
      Tr, 
      Surface, 
      Criteria, 
      Base
    > Self;
    typedef Mesher_level <
      Tr,
      Self,
      typename Tr::Facet,
      Null_mesher_level,
      Triangulation_mesher_level_traits_3<Tr>
    > Mesher;
    
    typedef Complex_2_in_triangulation_3_surface_mesh<Tr> C2t3;

    using Mesher::scan_triangulation;
    using Mesher::refine;
    using Mesher::is_algorithm_done;
    using Mesher::one_step;
    using Base::check_restricted_delaunay;
    

  protected:
    Null_mesher_level null_mesher_level;
    Null_mesh_visitor null_visitor;
    bool initialized;
    
  public:
    Surface_mesher_regular_edges_without_boundary(Tr& t, 
					      C2t3& co, 
					      Surface& s, 
					      Criteria& c): 
      Base(t, co, s, c), 
      Mesher(null_mesher_level),
      initialized(false)
      {}
    
    
    // Initialization
    void init(bool debug = false)
      {
	scan_triangulation();
	initialized = true;
	if (debug)
	  check_restricted_delaunay();
      }
    
    void refine_mesh (bool verbose = false, bool debug = false) {
      if(!initialized)
	init (debug);
      
      
      if (!verbose)
	refine (null_visitor);
      else {  
	std::cerr << "Refining...\n";
	int nbsteps = 0;
	std::cerr << "(" << nbsteps << "," 
		  << this->facets_to_refine.size() << ")";
	while (!is_algorithm_done()) {
	  one_step (null_visitor);
	  std::cerr << "\r             \r"
		    << "(" << ++nbsteps << "," 
		    << this->facets_to_refine.size() 
		    << "," 
		    << this->bad_edges.size() 
		    << ")";
	}
	std::cerr << "\ndone.\n";
      }
      
      if (debug)
	check_restricted_delaunay();
      
      initialized = false;
    }
    
  };  // end Surface_mesher_regular_edges_without_boundary
      
  }  // end namespace Surface_mesher
  
}  // end namespace CGAL

  
#endif // CGAL_SURFACE_MESHER_REGULAR_EDGES_H
  

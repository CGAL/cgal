#ifndef CGAL_SURFACE_MESHER_REGULAR_EDGES_H
#define CGAL_SURFACE_MESHER_REGULAR_EDGES_H

#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3_surface_mesh.h>
#include <CGAL/Surface_mesher/Surface_mesher.h>
#include <CGAL/utility.h>

namespace CGAL {

  namespace Surface_mesher {

  template < typename Tr, 
    typename Surface,
    typename Criteria,
    typename SMB = Surface_mesher_base <Tr, Surface, Criteria> >
  class Surface_mesher_regular_edges_base
    : public SMB {
      
    public:
      typedef Complex_2_in_triangulation_3_surface_mesh<Tr> C2t3;
      typedef typename Tr::Point Point;
    //typedef typename Tr::Edge Edge;
      typedef typename Tr::Facet Facet;
      typedef typename Tr::Edge Edge;
      typedef typename Tr::Cell Cell;
      typedef typename Tr::Cell_handle Cell_handle;
      typedef typename Tr::Vertex_handle Vertex_handle;
      typedef std::pair<Vertex_handle, Vertex_handle> EdgeVV;
      //      typedef typename Tr::Facet_circulator Facet_circulator;
      typedef typename Triangulation_mesher_level_traits_3<Tr>::Zone Zone;
      typedef typename C2t3::Face_type Face_type;
      typedef std::list<Edge> Edges;
      typedef std::list<Facet> Facets;
      typedef std::list<Cell> Cells;
      typedef Const_circulator_from_container<Facets> Facet_circulator;
      typedef typename Tr::Finite_edges_iterator Finite_edges_iterator;

    protected:
      set <pair <Vertex_handle, Vertex_handle> > bad_edges;

    protected:
      bool is_in_complex(const C2t3& c, const Facet& f) const {
	return (c.complex_subface_type(f) != C2t3::NOT_IN_COMPLEX);
      }

      // computes and return an ordered pair of Vertex
      EdgeVV make_edgevv(const Vertex_handle vh1, 
			 const Vertex_handle vh2) const {
	if (vh1 < vh2) {
	  return std::make_pair(vh1, vh2);
	}
	else {
	  return std::make_pair(vh2, vh1);
	}
      }

      // computes and returns the Edge type object from the EdgeVV object
      // use tr.is_edge(...)
      Edge edgevv_to_edge(const EdgeVV& arete) const {
    	Vertex_handle sommet1 = arete.first;
    	Vertex_handle sommet2 = arete.second;
    	Cell_handle c;
	int index1, index2;

	SMB::tr.is_edge(sommet1, sommet2, c, index1, index2);

    	return make_triple( c, index1, index2 );
      }

      // computes and returns the EdgeVV type object from the Edge object
      EdgeVV edge_to_edgevv(const Edge& arete) const {
	return make_edgevv(arete.first->vertex(arete.second), 
			   arete.first->vertex(arete.third) );
      }

      // computes and returns the opposite side of a given Facet side
      Facet other_side(const Facet& cote) const {
	return Facet
	  ((cote.first)->neighbor(cote.second),
	   (cote.first)->neighbor(cote.second)->
	   index(cote.first));
      }

      // computes and returns the list of Edges in a given Facet side
      Edges edges_in_facet(const Facet& cote) const {

	Edges loe;
	Cell_handle c = cote.first;
	int i = cote.second;
	
	for (int j = 0; j < 3; j++) {
	  for (int k = j + 1; k < 4; k++) {
	    if ( (i != j) && (i != k) ){
	      loe.push_back(make_triple(c, j, k));
	    }
	  }
	}
	return loe;
      }

      Facet biggest_incident_facet_in_complex(const Edge& arete) const {
	// Find the first facet in the incident facets 
        // of the edge which is in the Complex 
	// use the list of incident facets in the complex
	Facet_circulator fcirc = SMB::c2t3.incident_facets(arete);
	Facet first_facet = *fcirc;
	Facet biggest_facet = first_facet;
	
	for (++fcirc; 
	     (*fcirc != first_facet) && 
	       (*fcirc != other_side(first_facet)); 
	     ++fcirc) {
	  Vertex_handle fev = edge_to_edgevv(arete).first;
	  // is the current facet bigger than the current biggest one
	  if ( SMB::c2t3.compute_squared_distance(*fcirc, fev) > 
	       SMB::c2t3.compute_squared_distance(biggest_facet, fev) ) {
	    biggest_facet = *fcirc;
	  }
	  else {
	    Facet autre_cote = other_side(*fcirc);
	    // is the current facet bigger than the current biggest one
	    if ( SMB::c2t3.compute_squared_distance(autre_cote, fev) > 
		 SMB::c2t3.compute_squared_distance(biggest_facet, fev) ) {
	      biggest_facet = autre_cote;
	    }
	  }
	}
	return biggest_facet;
      }

      ///////////////////////
      // For before_insertion
    
      // Actions to perform on a facet inside the conflict zone
      void handle_facet_inside_conflict_zone (const Facet& f) {
	Facet cote = f;
	Facet autre_cote = other_side(cote);
	
	if ( SMB::c2t3.complex_subface_type(cote) != 
	     C2t3::NOT_IN_COMPLEX ) {
	  
	  Edges loe = edges_in_facet(cote);
	  // foreach edge of cote
	  for ( typename Edges::iterator eit = loe.begin();
		eit != loe.end(); 
		++eit) {
	    bad_edges.erase( edge_to_edgevv(*eit) );
	  }
	}
      }
    
      // Action to perform on a facet on the boundary of the conflict zone
      void handle_facet_on_boundary_of_conflict_zone (const Facet& f) {
	// perform the same operations as for an internal facet
	handle_facet_inside_conflict_zone (f);
      }

    public:
      Surface_mesher_regular_edges_base(Tr& t, C2t3& co, Surface& s, Criteria& c)
      :SMB(t, co, s, c){}

    // Initialization function
    void scan_triangulation_impl() {
      scan_triangulation_impl(true);
    }

    void scan_triangulation_impl(const bool withBoundary) {
      SMB::scan_triangulation_impl();
      std::cout << "scanning edges..." << std::endl;
      for (Finite_edges_iterator eit = SMB::tr.finite_edges_begin(); eit != 
	     SMB::tr.finite_edges_end(); ++eit) {
	if ( (SMB::c2t3.complex_subface_type(*eit) 
	      == C2t3::SINGULAR) || 
	     ( (!withBoundary) && 
	       (SMB::c2t3.complex_subface_type(*eit) 
		== C2t3::BOUNDARY) ) ) {
	  bad_edges.insert( edge_to_edgevv(*eit) );
	}
      }
    }

    // Tells whether there remain elements to refine
    bool no_longer_element_to_refine_impl() const {
      return (SMB::facets_to_refine.empty() && bad_edges.empty());
    }

    // Returns the next element to refine
    Facet get_next_element_impl() {

      if (!SMB::facets_to_refine.empty()) {
	return SMB::facets_to_refine.front()->second;
      }
      else {
	Edge first_bad_edge = edgevv_to_edge(*(bad_edges.begin()));
	return biggest_incident_facet_in_complex(first_bad_edge);
      }
    }

    void before_insertion_impl(const Facet&, const Point& s,
			       Zone& zone) {
      for (typename Zone::Facets_iterator fit = 
	     zone.internal_facets.begin(); 
	   fit != zone.internal_facets.end();
	   ++fit)
	handle_facet_inside_conflict_zone (*fit);
      
      for (typename Zone::Facets_iterator fit =
	     zone.boundary_facets.begin(); fit !=
	     zone.boundary_facets.end(); ++fit)
	handle_facet_on_boundary_of_conflict_zone (*fit);
      
      SMB::before_insertion_impl(Facet(), s, zone);
    }

    void after_insertion_impl(const Vertex_handle v, 
			      const bool withBoundary) {
      SMB::after_insertion_impl(v);

      // search for incident facets around v
      Facets facets;
      SMB::tr.incident_facets (v, std::back_inserter(facets));

      // foreach f in star(v)
      for (typename Facets::iterator fit = facets.begin(); 
	   fit != facets.end(); 
	   ++fit) { 
	Facet cote = *fit;
	Edges loe = edges_in_facet(cote);

	// foreach edge of cote
	for ( typename Edges::iterator eit = loe.begin();
	      eit != loe.end();
	      ++eit ) {
	  // test if edge is in Complex
	  if ( SMB::c2t3.complex_subface_type(*eit) 
	       != C2t3::NOT_IN_COMPLEX ) {
	    // test if edge is not regular to store it as a "bad_edge"
	    // e.g. more than or equal to 3 incident facets (SINGULAR) 
	    // or less than or equal to 1 
	    // (BOUNDARY only, because ISOLATED is NA)
	    // This test is not efficient because
	    // edges are tried to be inserted several times
	    // TODO one day: test if the edge is still singular
	    if ( (SMB::c2t3.complex_subface_type(*eit) 
		  == C2t3::SINGULAR) || 
		 ( (!withBoundary) && 
		   (SMB::c2t3.complex_subface_type(*eit) 
		    == C2t3::BOUNDARY) ) ) {
	      bad_edges.insert( edge_to_edgevv(*eit) );
	    }
	    else {
	      bad_edges.erase( edge_to_edgevv(*eit) );
	    }
	  }
	}
      }
    }
    
    void after_insertion_impl(const Vertex_handle v) {
      after_insertion_impl(v, true);
    }
  };  // end Surface_mesher_base

  template <typename Tr,
    typename Surface,
    typename Criteria,
    typename Base = Surface_mesher_regular_edges_base<Tr, Surface, Criteria> >
  class Surface_mesher_regular_edges : 
      public Base, 
      public Mesher_level <
    Tr,
    Surface_mesher_regular_edges<Tr, Surface, Criteria, Base>,
    typename Tr::Facet,
    Null_mesher_level,
    Triangulation_mesher_level_traits_3<Tr>
  >
  {
  public:
    typedef Surface_mesher_regular_edges<Tr, Surface, Criteria, Base> Self;
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
    Surface_mesher_regular_edges(Tr& t, C2t3& co, Surface& s, Criteria& c): 
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
    
  };  // end Surface_mesher_regular_edges
      
  }  // end namespace Surface_mesher
  
}  // end namespace CGAL

  
#endif // CGAL_SURFACE_MESHER_REGULAR_EDGES_H
  

// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Apollonius_graph_2.h
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================



#ifndef CGAL_APOLLONIUS_GRAPH_2_H
#define CGAL_APOLLONIUS_GRAPH_2_H

#include <vector>
#include <map>

#include <CGAL/Triangulation_2.h>
#include <CGAL/Apollonius_graph_data_structure_2.h>
#include <CGAL/Apollonius_graph_face_base_2.h>
#include <CGAL/Apollonius_graph_vertex_base_2.h>

#include <CGAL/in_place_edge_list.h>
#include <CGAL/Apollonius_graph_euclidean_traits_wrapper_2.h>


#include <CGAL/Iterator_project.h>
#include <CGAL/Nested_iterator.h>
#include <CGAL/Concatenate_iterator.h>


CGAL_BEGIN_NAMESPACE

template < class Node >
struct Project_site_2 {
  typedef Node                   argument_type;
  typedef typename Node::Site_2  Site;
  typedef Site                   result_type;
  Site&       operator()( Node& x) const { return x.site(); }
  const Site& operator()( const Node& x) const { return x.site(); }
};

template < class Gt,
	   class Agds = Apollonius_graph_data_structure_2 < 
               Apollonius_graph_vertex_base_2<Gt,true>,
               Apollonius_graph_face_base_2<Gt> > >
class Apollonius_graph_2
  : private Triangulation_2< Apollonius_graph_gt_wrapper<Gt>, Agds >
{
protected:
  // some local types
  typedef Apollonius_graph_gt_wrapper<Gt>        Modified_traits;
  typedef Triangulation_2<Modified_traits,Agds>  DG;

  typedef DG                                 Delaunay_graph;
  typedef typename DG::Vertex                Vertex;
  typedef typename DG::Face                  Face;

public:
  // TYPES
  //------
  typedef Agds                                   Data_structure;
  typedef Gt                                     Geom_traits;
  typedef typename Gt::Point_2                   Point_2;
  typedef typename Gt::Site_2                    Site_2;

  typedef typename DG::Edge                      Edge;
  typedef typename DG::Vertex_handle             Vertex_handle;
  typedef typename DG::Face_handle               Face_handle;

  typedef typename DG::Vertex_circulator         Vertex_circulator;
  typedef typename DG::Edge_circulator           Edge_circulator;
  typedef typename DG::Face_circulator           Face_circulator;

  typedef typename DG::All_faces_iterator        All_faces_iterator;
  typedef typename DG::Finite_faces_iterator     Finite_faces_iterator;
  typedef typename DG::All_vertices_iterator     All_vertices_iterator;
  typedef typename DG::Finite_vertices_iterator  Finite_vertices_iterator;
  typedef typename DG::All_edges_iterator        All_edges_iterator;
  typedef typename DG::Finite_edges_iterator     Finite_edges_iterator;


  // Auxiliary iterators for convenience
  // do not use default template argument to please VC++
  typedef Project_site_2<Vertex>                               Proj_site;
  typedef Iterator_project<Finite_vertices_iterator, 
                           Proj_site>
  /*                                           */ Visible_sites_iterator;

  typedef
  Apollonius_graph_vertex_base_nested_iterator_traits<
    Finite_vertices_iterator>  Hidden_sites_nested_iterator_traits;
  

  typedef Nested_iterator<Finite_vertices_iterator,
			  Hidden_sites_nested_iterator_traits>
  /*                                            */ Hidden_sites_iterator;

  typedef Concatenate_iterator<Visible_sites_iterator,
			       Hidden_sites_iterator>     Sites_iterator;

  typedef Site_2               value_type; // to have a back_inserter
  typedef const value_type&    const_reference; 
  typedef value_type&          reference;
private:
  typedef struct   Vertex_iterator {};
  typedef struct   Face_iterator {};
  typedef struct   Edge_iterator {};

protected:
  // some more local types
  typedef typename Agds::Vertex            Vertex;

  // point lists
  typedef std::vector<Site_2>              Site_list;
  typedef typename Site_list::iterator     Site_list_iterator;

  typedef std::map<Face_handle,bool>           Face_map;
  typedef std::map<Face_handle, Face_handle>   Face_face_map;
  typedef std::map<Vertex_handle,bool>         Vertex_map;
  typedef std::vector<Edge>                    Edge_list;

  typedef std::list<Vertex_handle>         Vertex_list;
  typedef typename Vertex_list::iterator   Vertex_list_iterator;
  typedef Vertex_handle                    Vh_triple[3];

  // the in place edge list
  typedef In_place_edge_list<Edge>          List;

  typedef enum { NO_CONFLICT = -1, INTERIOR, LEFT_VERTEX,
		 RIGHT_VERTEX, BOTH_VERTICES, ENTIRE_EDGE }
  Conflict_type;

  static Conflict_type opposite(const Conflict_type& ct) {
    if ( ct == RIGHT_VERTEX ) { return LEFT_VERTEX; }
    if ( ct == LEFT_VERTEX ) { return RIGHT_VERTEX; }
    return ct;
  }

protected:
  // Less_than comparator for weights of weighted points;
  // used to sort sites by decreasing weight when a sequence of sites
  // is inserted
  class Site_less_than_comparator
  {
  private:
    const Gt& gt;
  public:
    Site_less_than_comparator(const Gt& gt) : gt(gt) {}

    bool operator ()(const Site_2& p,
		     const Site_2& q) {
      Comparison_result result = gt.compare_weight_2_object()(p, q);
      return (result == LARGER);
    }
  };

public:
  // CREATION
  //---------
  Apollonius_graph_2(const Gt& gt=Gt()) :
    DG( Modified_traits(gt) ) {}

  template< class Input_iterator >
  Apollonius_graph_2(Input_iterator first, Input_iterator beyond,
		     const Gt& gt=Gt())
    : DG( Modified_traits(gt) )
  {
    insert(first, beyond);
  }

  Apollonius_graph_2(const Apollonius_graph_2 &ag)
    : DG(ag)
  {
    CGAL_postcondition( is_valid() );
  }

  Apollonius_graph_2&
  operator=(const Apollonius_graph_2& ag)
  {
    if ( this != &ag ) {
      DG::operator=(ag);
    }
    return (*this);
  }

public:
  // ACCESS METHODS
  // --------------
  inline const Geom_traits& geom_traits() const {
    return DG::geom_traits();
  }

  inline int dimension() const {
    return this->_tds.dimension();
  }

  inline int number_of_faces() const {
    return this->_tds.number_of_faces();
  }

  inline int number_of_vertices() const {
    return DG::number_of_vertices();
  }

  inline int number_of_hidden_vertices() const {
    //    if ( !Vertex::StoreHidden ) { return 0; }

    int n_hidden(0);
    for (Finite_vertices_iterator vit = finite_vertices_begin();
	 vit != finite_vertices_end(); ++vit) {
      n_hidden += vit->number_of_hidden_sites();
    }

    return n_hidden;
  }

  inline Vertex_handle infinite_vertex() const {
    return DG::infinite_vertex();
  }

  inline Face_handle infinite_face() const {
    return DG::infinite_face();
  }

  inline Vertex_handle finite_vertex() const {
    return DG::finite_vertex();
  }

public:
  // TRAVERSAL OF THE APOLLONIUS GRAPH
  //----------------------------------
  Finite_faces_iterator finite_faces_begin() const {
    return DG::finite_faces_begin();
  }

  Finite_faces_iterator finite_faces_end() const {
    return DG::finite_faces_end();
  }

  Finite_vertices_iterator finite_vertices_begin() const {
    return DG::finite_vertices_begin();
  }

  Finite_vertices_iterator finite_vertices_end() const {
    return DG::finite_vertices_end();
  }

  Finite_edges_iterator finite_edges_begin() const {
    return DG::finite_edges_begin();    
  }
  Finite_edges_iterator finite_edges_end() const {
    return DG::finite_edges_end();    
  }


  Sites_iterator sites_begin() const {
    return Sites_iterator(visible_sites_end(),
			  hidden_sites_begin(),
			  visible_sites_begin());
  }

  Sites_iterator sites_end() const {
    return Sites_iterator(visible_sites_end(),
			  hidden_sites_begin(),
			  hidden_sites_end(),0);
  }

  Visible_sites_iterator visible_sites_begin() const {
    return Visible_sites_iterator(finite_vertices_begin());
  }

  Visible_sites_iterator visible_sites_end() const {
    return Visible_sites_iterator(finite_vertices_end());
  }

  Hidden_sites_iterator hidden_sites_begin() const {
    return Hidden_sites_iterator(finite_vertices_begin(),
				 finite_vertices_end(),
				 finite_vertices_begin());
  }

  Hidden_sites_iterator hidden_sites_end() const {
    return Hidden_sites_iterator(finite_vertices_begin(),
				 finite_vertices_end(),
				 finite_vertices_end());
  }


  All_faces_iterator all_faces_begin() const {
    return DG::all_faces_begin();
  }

  All_faces_iterator all_faces_end() const {
    return DG::all_faces_end();
  }

  All_vertices_iterator all_vertices_begin() const {
    return DG::all_vertices_begin();
  }

  All_vertices_iterator all_vertices_end() const {
    return DG::all_vertices_end();
  }

  All_edges_iterator all_edges_begin() const {
    return DG::all_edges_begin();
  }

  All_edges_iterator all_edges_end() const {
    return DG::all_edges_end();
  }

public:
  // CIRCULATORS
  //------------
  // I had to add the initialization of Face_handle to NULL because
  // CGAL-2.5-I-82 was not working with Face_handle()
  inline Face_circulator
  incident_faces(Vertex_handle v,
		 Face_handle f = Face_handle(NULL)) const {
    return DG::incident_faces(v, f);
  }

  inline Vertex_circulator
  incident_vertices(Vertex_handle v,
		    Face_handle f = Face_handle(NULL)) const { 
    return DG::incident_vertices(v, f);
  }

  inline Edge_circulator
  incident_edges(Vertex_handle v,
		 Face_handle f = Face_handle(NULL)) const {
    return DG::incident_edges(v, f);
  }
 
public:
  // PREDICATES
  //-----------
  inline bool is_infinite(const Vertex_handle& v) const {
    return DG::is_infinite(v);
  }

  inline bool is_infinite(const Face_handle& f) const {
    return DG::is_infinite(f);
  }

  inline bool is_infinite(const Face_handle& f, int i) const {
    return DG::is_infinite(f, i);
  }

  inline bool is_infinite(const Edge& e) const {
    return is_infinite(e.first, e.second);
  }

  inline bool is_infinite(const Edge_circulator& ec) const {
    return DG::is_infinite(ec);
  }

public:
  // INSERTION
  //----------
  template< class Input_iterator >
  void insert(Input_iterator first, Input_iterator beyond) {
    // copy to a local container
    Site_list wp_list;
    for (Input_iterator it = first; it != beyond; ++it) {
      wp_list.push_back(*it);
    }

    // sort by decreasing weight
    Site_less_than_comparator less_than(geom_traits());
    std::sort(wp_list.begin(), wp_list.end(), less_than);

    // now insert
    Site_list_iterator lit;
    for (lit = wp_list.begin(); lit != wp_list.end(); ++lit) {
      insert(*lit);
    }

    // clear the local container
    wp_list.clear();
  }

  Vertex_handle  insert(const Site_2& p);
  Vertex_handle  insert(const Site_2& p, Vertex_handle vnear);

public:
  // REMOVAL
  //--------
  void remove(Vertex_handle v);

public:
  // NEAREST NEIGHBOR LOCATION
  //--------------------------
  Vertex_handle  nearest_neighbor(const Point_2& p) const;
  Vertex_handle  nearest_neighbor(const Point_2& p,
				  Vertex_handle vnear) const;

public:
  // ACCESS TO THE DUAL
  //-------------------
  typename Gt::Object_2 dual(const Face_handle& f) const;

private:
  typename Gt::Object_2 dual(const Edge e) const;

  typename Gt::Object_2 dual(const Edge_circulator& ec) const {
    return dual(*ec);
  }

  typename Gt::Object_2 dual(const Finite_edges_iterator& ei) const {
    return dual(*ei);
  }

public:
  // I/O
  //----
  template< class Stream >
  inline
  Stream& draw_primal(Stream &str) const
  {
    if ( number_of_vertices() < 2 ) {
      // do nothing
    } else if ( number_of_vertices() == 2 ) {
      Vertex_handle v1(finite_vertices_begin());
      Vertex_handle v2(++finite_vertices_begin());
      Site_2 p1 = v1->site();
      Site_2 p2 = v2->site();
      typename Geom_traits::Segment_2 seg =
	geom_traits().construct_Apollonius_primal_segment_2_object()(p1,p2);
      typename Geom_traits::Ray_2 ray1 =
	geom_traits().construct_Apollonius_primal_ray_2_object()(p1,p2,p2);
      typename Geom_traits::Ray_2 ray2 =
	geom_traits().construct_Apollonius_primal_ray_2_object()(p2,p1,p1);

      str << seg;
      str << ray1;
      str << ray2;
    } else {
      All_edges_iterator eit = all_edges_begin();
      for (; eit != all_edges_end(); ++eit) {
	draw_primal_edge< Stream >(*eit, str);
      }
    }
    return str;
  }

  template < class Stream > 
  Stream& draw_dual(Stream &str) const
  {
    Finite_edges_iterator eit = finite_edges_begin();
    for (; eit != finite_edges_end(); ++eit) {
      typename Gt::Object_2 o = dual(eit);
      typename Geom_traits::Line_2     l;
      typename Geom_traits::Segment_2  s;
      typename Geom_traits::Ray_2      r;
      typename Geom_traits::Hyperbola_2          h;
      typename Geom_traits::Hyperbola_segment_2  hs;
      typename Geom_traits::Hyperbola_ray_2      hr;
      if (assign(hs, o))  str << hs;
      if (assign(s, o))   str << s; 
      if (assign(hr, o))  str << hr;
      if (assign(r, o))   str << r;
      if (assign(h, o))   str << h;
      if (assign(l, o))   str << l;
    }
    return str;
  }

protected:
  template < class Stream > 
  Stream& draw_dual_weighted_points(Stream &str) const
  {
    All_faces_iterator fit = all_faces_begin();
    for (; fit != all_faces_end(); ++fit) {
      Face_handle f(fit);
      if ( is_infinite(f) ) {
	if (  is_infinite(f->vertex(0))  ) {
	  str << circumcircle( f->vertex(1)->site(),
			       f->vertex(2)->site() );
	} else if (  is_infinite(f->vertex(1))  ){
	  str << circumcircle( f->vertex(2)->site(),
			       f->vertex(0)->site() );
	} else {
	  str << circumcircle( f->vertex(0)->site(),
			       f->vertex(1)->site() );	  
	}
      } else {
	Site_2 wp = circumcircle(f);
	typename Gt::Rep::Circle_2 c(wp.point(),
				     CGAL_NTS square(wp.weight()));
	str << c;
      }
    }
    return str;
  }

public:
  // VALIDITY CHECK
  //---------------
  bool is_valid(bool verbose = false, int level = 1) const;

public:
  // MISCELLANEOUS
  //--------------
  inline void clear() {
    DG::clear();
  }

  inline void swap(Apollonius_graph_2& ag) {
    DG::swap(ag);
  }

public:
  // MK: THE FOLLOWING ARE NOT IN THE SPEC
  //======================================
  // Primal
  //  Weighted_point primal(const Face_handle& f) const;
  typename Gt::Object_2 primal(const Edge e) const;
  typename Gt::Object_2 primal(const Edge_circulator& ec) const {
    return primal(*ec);
  }
  typename Gt::Object_2 primal(const Finite_edges_iterator& ei) const {
    return primal(*ei);
  }

protected:
  // wrappers for the geometric predicates

  // checks is q is contained inside p
  bool is_hidden(const Site_2 &p,
		 const Site_2 &q) const;

  // returns:
  //   ON_POSITIVE_SIDE if q is closer to p1
  //   ON_NEGATIVE_SIDE if q is closer to p2
  //   ON_ORIENTED_BOUNDARY if q is on the bisector of p1 and p2
  Oriented_side side_of_bisector(const Site_2 &p1,
				 const Site_2 &p2,
				 const Point_2 &q) const;

  Sign incircle(const Site_2 &p1, const Site_2 &p2,
		const Site_2 &p3, const Site_2 &q) const;

  Sign incircle(const Site_2 &p1, const Site_2 &p2,
		const Site_2 &q) const;


  Sign incircle(const Face_handle& f, const Site_2& q) const;


  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v) const;

  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Vertex_handle& v) const;


  
  bool finite_edge_interior(const Site_2& p1,
			    const Site_2& p2,
			    const Site_2& p3,
			    const Site_2& p4,
			    const Site_2& q,
			    bool endpoints_in_conflict) const;

  bool finite_edge_interior(const Face_handle& f, int i,
			    const Site_2& q,
			    bool endpoints_in_conflict) const;

  bool finite_edge_interior(const Vertex_handle& v1,
			    const Vertex_handle& v2,
			    const Vertex_handle& v3,
			    const Vertex_handle& v4,
			    const Vertex_handle& v,
			    bool endpoints_in_conflict) const;

  bool finite_edge_interior_degenerated(const Site_2& p1,
					const Site_2& p2,
					const Site_2& p3,
					const Site_2& q,
					bool endpoints_in_conflict) const;


  bool finite_edge_interior_degenerated(const Site_2& p1,
					const Site_2& p2,
					const Site_2& q,
					bool endpoints_in_conflict) const;

  bool finite_edge_interior_degenerated(const Face_handle& f, int i,
					const Site_2& p,
					bool endpoints_in_conflict) const;

  bool finite_edge_interior_degenerated(const Vertex_handle& v1,
					const Vertex_handle& v2,
					const Vertex_handle& v3,
					const Vertex_handle& v4,
					const Vertex_handle& v,
					bool endpoints_in_conflict) const;
  bool infinite_edge_interior(const Site_2& p2,
			      const Site_2& p3,
			      const Site_2& p4,
			      const Site_2& q,
			      bool endpoints_in_conflict) const;


  bool infinite_edge_interior(const Face_handle& f, int i,
			      const Site_2& p,
			      bool endpoints_in_conflict) const;

  bool infinite_edge_interior(const Vertex_handle& v1,
			      const Vertex_handle& v2,
			      const Vertex_handle& v3,
			      const Vertex_handle& v4,
			      const Vertex_handle& v,
			      bool endpoints_in_conflict) const;

  Conflict_type
  finite_edge_conflict_type_degenerated(const Site_2& p1,
					const Site_2& p2,
					const Site_2& q) const;

  bool edge_interior(const Face_handle& f, int i,
		     const Site_2& p, bool b) const;


  inline bool edge_interior(const Edge& e,
			    const Site_2& p, bool b) const {
    return edge_interior(e.first, e.second, p, b);
  }

  bool edge_interior(const Vertex_handle& v1,
		     const Vertex_handle& v2,
		     const Vertex_handle& v3,
		     const Vertex_handle& v4,
		     const Vertex_handle& v,
		     bool endpoints_in_conflict) const;

  inline bool is_degenerate_edge(const Site_2& p1,
				 const Site_2& p2,
				 const Site_2& p3,
				 const Site_2& p4) const {
    return geom_traits().is_degenerate_edge_2_object()
      (p1, p2, p3, p4);
  }

  inline bool is_degenerate_edge(const Vertex_handle& v1,
				 const Vertex_handle& v2,
				 const Vertex_handle& v3,
				 const Vertex_handle& v4) const {
    CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) &&
		       !is_infinite(v3) && !is_infinite(v4) );

    return is_degenerate_edge(v1->site(), v2->site(),
			      v3->site(), v4->site());
  }

  inline bool is_degenerate_edge(const Face_handle& f, int i) const {
    Vertex_handle v1 = f->vertex( ccw(i) );
    Vertex_handle v2 = f->vertex(  cw(i) );
    Vertex_handle v3 = f->vertex(     i  );
    Vertex_handle v4 = f->mirror_vertex(i);

    return is_degenerate_edge(v1, v2, v3, v4);
  }

  inline bool is_degenerate_edge(const Edge& e) const {
    return is_degenerate_edge(e.first, e.second);
  }


protected:
  // wrappers for constructions
  Point_2 circumcenter(const Face_handle& f) const;
  Point_2 circumcenter(const Site_2& p0, 
		       const Site_2& p1, 
		       const Site_2& p2) const;

  Site_2 circumcircle(const Face_handle& f) const;
  Site_2 circumcircle(const Site_2& p0, 
		      const Site_2& p1, 
		      const Site_2& p2) const;

  typename Gt::Line_2 circumcircle(const Site_2& p0,
				   const Site_2& p1) const;

protected:
  // wrappers for combinatorial operations on the data structure

  // getting the symmetric edge
  inline Edge sym_edge(const Edge e) const {
    return sym_edge(e.first, e.second);
  }

  inline Edge sym_edge(const Face_handle& f, int i) const {
    Face_handle f_sym = f->neighbor(i);
    return Edge(  f_sym, f_sym->index( f->mirror_vertex(i) )  );
  }

  Edge flip(Face_handle& f, int i);
  Edge flip(Edge e);

  Vertex_handle insert_in_face(Face_handle& f, const Site_2& p);

  bool          is_degree_2(const Vertex_handle& v) const;

  Vertex_handle insert_degree_2(Edge e);
  Vertex_handle insert_degree_2(Edge e, const Site_2& p);
  void          remove_degree_2(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v, Face* f);

  // this was defined because the hierarchy needs it
  inline Vertex_handle create_vertex() {
    return this->_tds.create_vertex();
  }

protected:
  // insertion of the first three sites
  Vertex_handle  insert_first(const Site_2& p);
  Vertex_handle  insert_second(const Site_2& p);
  Vertex_handle  insert_third(const Site_2& p);

  // methods for insertion
  void initialize_conflict_region(const Face_handle& f, List& l);
  bool check_edge_for_hidden_sites(const Face_handle& f, int i,
				   const Site_2& p, Vertex_map& vm);
  void expand_conflict_region(const Face_handle& f,
			      const Site_2& p,
			      List& l, Face_map& fm, Vertex_map& vm,
			      std::vector<Vh_triple*>* fe);

  Vertex_handle add_bogus_vertex(Edge e, List& l);
  Vertex_list   add_bogus_vertices(List& l);
  void          remove_bogus_vertices(Vertex_list& vl);

  void move_hidden_sites(Vertex_handle& vold, Vertex_handle& vnew);

  // MK: this is not currently used
  inline  std::vector<Face*> get_faces_for_recycling(Face_map& fm,
					     unsigned int n_wanted);
  void remove_hidden_vertices(Vertex_handle&v, Vertex_map& vm,
			      Face_map& fm);
  Vertex_handle retriangulate_conflict_region(const Site_2& p,
					      List& l,
					      Face_map& fm,
					      Vertex_map& vm);

protected:
  // methods for removal
  void  remove_first(Vertex_handle v);
  void  remove_second(Vertex_handle v);
  void  remove_third(Vertex_handle v);
  void  remove_degree_d_vertex(Vertex_handle v);
  void  minimize_degree(Vertex_handle v);

  void find_conflict_region_remove(const Vertex_handle& v,
				   const Vertex_handle& vnearest,
				   List& l, Face_map& fm,
				   Vertex_map& vm,
				   std::vector<Vh_triple*>* fe);

protected:
  // methods for I/O

  template<class T>
  bool assign(T& t2, const typename Gt::Object_2& o2) const
  {
    return geom_traits().assign_2_object()(t2, o2);
  }

  template< class Stream >
  Stream& draw_primal_edge(Edge e, Stream &str) const
  {
    Object o = primal(e);
    //      typename Geom_traits::Line_2     l;
    typename Geom_traits::Segment_2  s;
    typename Geom_traits::Ray_2      r;
    //      typename Geom_traits::Hyperbola_2          h;
    typename Geom_traits::Hyperbola_segment_2  hs;
    //      typename Geom_traits::Hyperbola_ray_2      hr;
    typename Geom_traits::Parabola_segment_2   ps;
    if (assign(hs, o))  str << hs;
    if (assign(s, o))   str << s; 
    if (assign(ps, o))  str << ps;
    if (assign(r, o))   str << r;
    //      if (assign(hr, o))  str << hr;
    //      if (assign(h, o))   str << h;
    //      if (assign(l, o))   str << l;
    return str;
  }

  template < class Stream > 
  Stream& draw_dual_edge(Edge e, Stream &str) const
  {
    if ( is_infinite(e) ) { return str; }
    Object o = dual(e);
    typename Geom_traits::Line_2     l;
    typename Geom_traits::Segment_2  s;
    typename Geom_traits::Ray_2      r;
    typename Geom_traits::Hyperbola_2          h;
    typename Geom_traits::Hyperbola_segment_2  hs;
    typename Geom_traits::Hyperbola_ray_2      hr;
    if (assign(hs, o))  str << hs;
    if (assign(s, o))   str << s; 
    if (assign(hr, o))  str << hr;
    if (assign(r, o))   str << r;
    if (assign(h, o))   str << h;
    if (assign(l, o))   str << l;

    return str;
  }

  template< class Stream >
  inline
  Stream& draw_face(const Face_handle& f, Stream &str) const
  {
    for (int i = 0; i < 3; i++) {
      draw_primal_edge< Stream >(Edge(f, i), str);
    }
    return str;
  }

}; // Apollonius_graph_2


CGAL_END_NAMESPACE


#include <CGAL/Apollonius_graph_2.C>



#endif // CGAL_APOLLONIUS_GRAPH_2_H

#ifndef CGAL_APOLLONIUS_GRAPH_2_H
#define CGAL_APOLLONIUS_GRAPH_2_H

#include <vector>
#include <map>

#include <CGAL/Triangulation_2.h>
//#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Apollonius_graph_data_structure_2.h>
#include <CGAL/Apollonius_graph_face_base_2.h>
#include <CGAL/Apollonius_graph_vertex_base_2.h>

#include <CGAL/in_place_edge_list.h>
//#include <CGAL/sorted_vertex_triple.h>
#include <CGAL/Apollonius_graph_euclidean_traits_wrapper_2.h>


CGAL_BEGIN_NAMESPACE


template < class Gt, bool StoreHidden = true,
  class Agds = Apollonius_graph_data_structure_2 < 
               Apollonius_graph_vertex_base_2<Gt,StoreHidden>,
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
  typedef typename Gt::Bare_point                Point;
  typedef typename Gt::Weighted_point            Weighted_point;
  typedef typename Gt::Weight                    Weight;

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

  typedef double   Vertex_iterator;
  typedef double   Face_iterator;
  typedef double   Edge_iterator;

protected:
  // some more local types
  typedef typename Agds::Vertex_base       Vertex_base;

  // point lists
  typedef std::vector<Weighted_point>      Weighted_point_list;
  typedef typename Weighted_point_list::iterator Weighted_point_list_iterator;

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

  static inline Conflict_type opposite(const Conflict_type& ct) {
    if ( ct == RIGHT_VERTEX ) { return LEFT_VERTEX; }
    if ( ct == LEFT_VERTEX ) { return RIGHT_VERTEX; }
    return ct;
  }

protected:
  // Less_than comparator for weights of weighted points;
  // used to sort sites by decreasing weight when a sequence of sites
  // is inserted
  class Weighted_point_less_than_comparator
  {
  private:
    const Gt& gt;
  public:
    Weighted_point_less_than_comparator(const Gt& gt) : gt(gt) {}

    bool operator ()(const Weighted_point& p,
		     const Weighted_point& q) {
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
    DG::operator=(ag);
    return (*this);
  }

public:
  // ACCESS METHODS
  // --------------
  inline const Geom_traits& geom_traits() const {
    return DG::geom_traits();
  }

  inline int number_of_vertices() const {
    return DG::number_of_vertices();
  }

  inline int number_of_hidden_vertices() const {
    if ( !StoreHidden ) { return 0; }

    int n_hidden(0);
    for (Finite_vertices_iterator vit = finite_vertices_begin();
	 vit != finite_vertices_end(); ++vit) {
      n_hidden += vit->number_of_hidden_weighted_points();
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
  inline Finite_faces_iterator finite_faces_begin() const {
    return DG::finite_faces_begin();
  }

  inline Finite_faces_iterator finite_faces_end() const {
    return DG::finite_faces_end();
  }

  inline Finite_vertices_iterator finite_vertices_begin() const {
    return DG::finite_vertices_begin();
  }

  inline Finite_vertices_iterator finite_vertices_end() const {
    return DG::finite_vertices_end();
  }

  inline Finite_edges_iterator finite_edges_begin() const {
    return DG::finite_edges_begin();    
  }
  inline Finite_edges_iterator finite_edges_end() const {
    return DG::finite_edges_end();    
  }

  //  Point_iterator points_begin() const;
  //  Point_iterator points_end() const;

  inline All_faces_iterator all_faces_begin() const {
    return DG::all_faces_begin();
  }

  inline All_faces_iterator all_faces_end() const {
    return DG::all_faces_end();
  }

  inline All_vertices_iterator all_vertices_begin() const {
    return DG::all_vertices_begin();
  }

  inline All_vertices_iterator all_vertices_end() const {
    return DG::all_vertices_end();
  }

  inline All_edges_iterator all_edges_begin() const {
    return DG::all_edges_begin();
  }

  inline All_edges_iterator all_edges_end() const {
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
    Weighted_point_list wp_list;
    for (Input_iterator it = first; it != beyond; ++it) {
      wp_list.push_back(*it);
    }

    // sort by decreasing weight
    Weighted_point_less_than_comparator less_than(geom_traits());
    std::sort(wp_list.begin(), wp_list.end(), less_than);

    // now insert
    Weighted_point_list_iterator lit;
    for (lit = wp_list.begin(); lit != wp_list.end(); ++lit) {
      insert(*lit);
    }

    // clear the local container
    wp_list.clear();
  }

  Vertex_handle  insert(const Weighted_point& p);
  Vertex_handle  insert(const Weighted_point& p, Vertex_handle vnear);

public:
  // REMOVAL
  //--------
  void remove(Vertex_handle v);

public:
  // NEAREST NEIGHBOR LOCATION
  //--------------------------
  Vertex_handle  nearest_neighbor(const Point& p) const;
  Vertex_handle  nearest_neighbor(const Point& p,
				  Vertex_handle vnear) const;

public:
  // ACCESS TO THE DUAL
  //-------------------
  Point  dual(const Face_handle& f) const;
  Object dual(const Edge e) const;

  inline Object dual(const Edge_circulator& ec) const {
    return dual(*ec);
  }

  inline Object dual(const Finite_edges_iterator& ei) const {
    return dual(*ei);
  }

public:
  // I/O
  //----
  template < class Stream >
  Stream& write_non_hidden_weighted_points(Stream& str) const
  {
    str << number_of_vertices() << std::endl;

    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      Weighted_point wp = vit->point();
      str << wp << std::endl;
    }
    return str;
  }

  template < class Stream >
  Stream& write_all_weighted_points(Stream& str) const
  {
    int n_total = number_of_vertices() + number_of_hidden_vertices();

    str << n_total << std::endl;

    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      Weighted_point wp = vit->point();
      str << wp << std::endl;
      // now write the hidden vertices
      if ( StoreHidden ) {
	typename Vertex_base::Hidden_weighted_point_iterator wpit;
	for (wpit = vit->hidden_weighted_points_begin();
	     wpit != vit->hidden_weighted_points_end(); ++wpit) {
	  str << (*wpit) << std::endl;
	}
      }
    }
    return str;
  }

  template < class Stream >
  Stream& draw_hidden_weighted_points(Stream &str) const
  {
    if ( !StoreHidden ) { return str; }

    Finite_vertices_iterator vit = finite_vertices_begin();

    for (; vit != finite_vertices_end(); ++vit) {
      if ( vit->number_of_hidden_weighted_points() > 0 ) {
	typename Vertex_base::Hidden_weighted_point_iterator wpit;
	for (wpit = vit->hidden_weighted_points_begin();
	     wpit != vit->hidden_weighted_points_end(); ++wpit) {
	  Weighted_point wp = *wpit;
	  typename Gt::Rep::Circle_2 c(wp.point(),
				       CGAL_NTS square(wp.weight()));
	  str << c;
	}
      }
    }
    return str;
  }

  template < class Stream >
  Stream& draw_hidden_weighted_point_centers(Stream &str) const
  {
    if ( !StoreHidden ) { return str; }

    Finite_vertices_iterator vit = finite_vertices_begin();

    for (; vit != finite_vertices_end(); ++vit) {
      if ( vit->number_of_hidden_weighted_points() > 0 ) {
	typename Vertex_base::Hidden_weighted_point_iterator wpit;
	for (wpit = vit->hidden_weighted_points_begin();
	     wpit != vit->hidden_weighted_points_end(); ++wpit) {
	  Weighted_point wp = *wpit;
	  str << wp.point();
	}
      }
    }
    return str;
  }

  template < class Stream >
  Stream& draw_non_hidden_weighted_points(Stream &str) const
  {
    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      Weighted_point wp = vit->point();
      typename Gt::Rep::Circle_2 c(wp.point(),
				   CGAL_NTS square(wp.weight()));
      str << c;
    }
    return str;
  }

  template < class Stream >
  Stream& draw_non_hidden_weighted_point_centers(Stream &str) const
  {
    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      Weighted_point wp = vit->point();
      str << wp.point();
    }
    return str;
  }

  template< class Stream >
  inline
  Stream& draw_primal(Stream &str) const
  {
    if ( number_of_vertices() < 2 ) {
      // do nothing
    } else if ( number_of_vertices() == 2 ) {
      Vertex_handle v1(finite_vertices_begin());
      Vertex_handle v2(++finite_vertices_begin());
      Weighted_point p1 = v1->point();
      Weighted_point p2 = v2->point();
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
      Object o = dual(eit);
      typename Geom_traits::Line_2     l;
      typename Geom_traits::Segment_2  s;
      typename Geom_traits::Ray_2      r;
      typename Geom_traits::Hyperbola_2          h;
      typename Geom_traits::Hyperbola_segment_2  hs;
      typename Geom_traits::Hyperbola_ray_2      hr;
      if (CGAL::assign(hs, o))  str << hs;
      if (CGAL::assign(s, o))   str << s; 
      if (CGAL::assign(hr, o))  str << hr;
      if (CGAL::assign(r, o))   str << r;
      if (CGAL::assign(h, o))   str << h;
      if (CGAL::assign(l, o))   str << l;
    }
    return str;
  }

  template < class Stream > 
  Stream& draw_dual_weighted_points(Stream &str) const
  {
    All_faces_iterator fit = all_faces_begin();
    for (; fit != all_faces_end(); ++fit) {
      Face_handle f(fit);
      if ( is_infinite(f) ) {
	if (  is_infinite(f->vertex(0))  ) {
	  str << circumcircle( f->vertex(1)->point(),
			       f->vertex(2)->point() );
	} else if (  is_infinite(f->vertex(1))  ){
	  str << circumcircle( f->vertex(2)->point(),
			       f->vertex(0)->point() );
	} else {
	  str << circumcircle( f->vertex(0)->point(),
			       f->vertex(1)->point() );	  
	}
      } else {
	Weighted_point wp = circumcircle(f);
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
    return DG::clear();
  }

  inline void swap(Apollonius_graph_2& ag) {
    DG::swap(ag);
  }

public:
  // MK: THE FOLLOWING ARE NOT IN THE SPEC
  //======================================
  // Primal
  //  Weighted_point primal(const Face_handle& f) const;
  Object primal(const Edge e) const;
  inline Object primal(const Edge_circulator& ec) const {
    return primal(*ec);
  }
  inline Object primal(const Finite_edges_iterator& ei) const {
    return primal(*ei);
  }


protected:
  // wrappers for the geometric predicates

  // checks is q is contained inside p
  bool is_hidden(const Weighted_point &p,
		 const Weighted_point &q) const;

  // returns:
  //   ON_POSITIVE_SIDE if q is closer to p1
  //   ON_NEGATIVE_SIDE if q is closer to p2
  //   ON_ORIENTED_BOUNDARY if q is on the bisector of p1 and p2
  Oriented_side side_of_bisector(const Weighted_point &p1,
				 const Weighted_point &p2,
				 const Point &q) const;

  Sign incircle(const Weighted_point &p1, const Weighted_point &p2,
		const Weighted_point &p3, const Weighted_point &q) const;

  Sign incircle(const Weighted_point &p1, const Weighted_point &p2,
		const Weighted_point &q) const;


  Sign incircle(const Face_handle& f, const Weighted_point& q) const;


  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v) const;

  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Vertex_handle& v) const;


  
  bool finite_edge_interior(const Weighted_point& p1,
			    const Weighted_point& p2,
			    const Weighted_point& p3,
			    const Weighted_point& p4,
			    const Weighted_point& q,
			    bool endpoints_in_conflict) const;

  bool finite_edge_interior(const Face_handle& f, int i,
			    const Weighted_point& q,
			    bool endpoints_in_conflict) const;

  bool finite_edge_interior(const Vertex_handle& v1,
			    const Vertex_handle& v2,
			    const Vertex_handle& v3,
			    const Vertex_handle& v4,
			    const Vertex_handle& v,
			    bool endpoints_in_conflict) const;

  bool finite_edge_interior_degenerated(const Weighted_point& p1,
					const Weighted_point& p2,
					const Weighted_point& p3,
					const Weighted_point& q,
					bool endpoints_in_conflict) const;


  bool finite_edge_interior_degenerated(const Weighted_point& p1,
					const Weighted_point& p2,
					const Weighted_point& q,
					bool endpoints_in_conflict) const;

  bool finite_edge_interior_degenerated(const Face_handle& f, int i,
					const Weighted_point& p,
					bool endpoints_in_conflict) const;

  bool finite_edge_interior_degenerated(const Vertex_handle& v1,
					const Vertex_handle& v2,
					const Vertex_handle& v3,
					const Vertex_handle& v4,
					const Vertex_handle& v,
					bool endpoints_in_conflict) const;
  bool infinite_edge_interior(const Weighted_point& p2,
			      const Weighted_point& p3,
			      const Weighted_point& p4,
			      const Weighted_point& q,
			      bool endpoints_in_conflict) const;


  bool infinite_edge_interior(const Face_handle& f, int i,
			      const Weighted_point& p,
			      bool endpoints_in_conflict) const;

  bool infinite_edge_interior(const Vertex_handle& v1,
			      const Vertex_handle& v2,
			      const Vertex_handle& v3,
			      const Vertex_handle& v4,
			      const Vertex_handle& v,
			      bool endpoints_in_conflict) const;

  Conflict_type
  finite_edge_conflict_type_degenerated(const Weighted_point& p1,
					const Weighted_point& p2,
					const Weighted_point& q) const;

  bool edge_interior(const Face_handle& f, int i,
		     const Weighted_point& p, bool b) const;


  inline bool edge_interior(const Edge& e,
			    const Weighted_point& p, bool b) const {
    return edge_interior(e.first, e.second, p, b);
  }

  bool edge_interior(const Vertex_handle& v1,
		     const Vertex_handle& v2,
		     const Vertex_handle& v3,
		     const Vertex_handle& v4,
		     const Vertex_handle& v,
		     bool endpoints_in_conflict) const;

  inline bool is_degenerate_edge(const Weighted_point& p1,
				 const Weighted_point& p2,
				 const Weighted_point& p3,
				 const Weighted_point& p4) const {
    return geom_traits().is_degenerate_edge_2_object()
      (p1, p2, p3, p4);
  }

  inline bool is_degenerate_edge(const Vertex_handle& v1,
				 const Vertex_handle& v2,
				 const Vertex_handle& v3,
				 const Vertex_handle& v4) const {
    CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) &&
		       !is_infinite(v3) && !is_infinite(v4) );

    return is_degenerate_edge(v1->point(), v2->point(),
			      v3->point(), v4->point());
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
  Point circumcenter(const Face_handle& f) const;
  Point circumcenter(const Weighted_point& p0, 
		     const Weighted_point& p1, 
		     const Weighted_point& p2) const;

  Weighted_point circumcircle(const Face_handle& f) const;
  Weighted_point circumcircle(const Weighted_point& p0, 
			      const Weighted_point& p1, 
			      const Weighted_point& p2) const;

  typename Gt::Line_2 circumcircle(const Weighted_point& p0,
				   const Weighted_point& p1) const;

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

  Vertex_handle insert_in_face(Face_handle& f, const Weighted_point& p);

  bool          is_degree_2(const Vertex_handle& v) const;

  Vertex_handle insert_degree_2(Edge e);
  Vertex_handle insert_degree_2(Edge e, const Weighted_point& p);
  void          remove_degree_2(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v, Face* f);

  // this was defined because the hierarchy needs it
  inline Vertex_handle create_vertex() {
    return _tds.create_vertex();
  }


protected:
  // insertion of the first three sites
  Vertex_handle  insert_first(const Weighted_point& p);
  Vertex_handle  insert_second(const Weighted_point& p);
  Vertex_handle  insert_third(const Weighted_point& p);

  // methods for insertion
  void initialize_conflict_region(const Face_handle& f, List& l);
  bool check_edge_for_hidden_weighted_points(const Face_handle& f, int i,
					     const Weighted_point& p,
					     Vertex_map& vm);
  void expand_conflict_region(const Face_handle& f, const Weighted_point& p,
			      List& l, Face_map& fm, Vertex_map& vm,
			      std::vector<Vh_triple*>* fe);

  Vertex_handle add_bogus_vertex(Edge e, List& l);
  Vertex_list   add_bogus_vertices(List& l);
  void          remove_bogus_vertices(Vertex_list& vl);

  void move_hidden_weighted_points(Vertex_handle& vold,
				   Vertex_handle& vnew);

  // MK: this is not currently used
  inline  std::vector<Face*> get_faces_for_recycling(Face_map& fm,
					     unsigned int n_wanted);
  void remove_hidden_vertices(Vertex_handle&v, Vertex_map& vm,
			      Face_map& fm);
  Vertex_handle retriangulate_conflict_region(const Weighted_point& p,
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
    if (CGAL::assign(hs, o))  str << hs;
    if (CGAL::assign(s, o))   str << s; 
    if (CGAL::assign(ps, o))  str << ps;
    if (CGAL::assign(r, o))   str << r;
    //      if (CGAL::assign(hr, o))  str << hr;
    //      if (CGAL::assign(h, o))   str << h;
    //      if (CGAL::assign(l, o))   str << l;
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
    if (CGAL::assign(hs, o))  str << hs;
    if (CGAL::assign(s, o))   str << s; 
    if (CGAL::assign(hr, o))  str << hr;
    if (CGAL::assign(r, o))   str << r;
    if (CGAL::assign(h, o))   str << h;
    if (CGAL::assign(l, o))   str << l;

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

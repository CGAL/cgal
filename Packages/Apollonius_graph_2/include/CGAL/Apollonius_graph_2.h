#ifndef CGAL_ADDITIVELY_WEIGHTED_VORONOI_DIAGRAM_2_H
#define CGAL_ADDITIVELY_WEIGHTED_VORONOI_DIAGRAM_2_H

#include <vector>
#include <CGAL/Point_2.h>
#include <CGAL/Weighted_point.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Hyperbola_2.h>
#include <CGAL/Hyperbola_segment_2.h>
#include <CGAL/Hyperbola_ray_2.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Additively_weighted_Voronoi_diagram_2_enum.h>
#include <CGAL/Additively_weighted_Voronoi_diagram_data_structure_2.h>
#include <CGAL/Additively_weighted_Voronoi_diagram_face_base_2.h>
#include <CGAL/Additively_weighted_Voronoi_diagram_vertex_base_2.h>

CGAL_BEGIN_NAMESPACE

//#define DO_PRINTOUTS 1
#if !defined DO_PRINTOUTS
#define DO_PRINTOUTS 0
#endif

template < class Gt,
  class Tds = Additively_weighted_Voronoi_diagram_data_structure_2 < 
              Additively_weighted_Voronoi_diagram_vertex_base_2<Gt>,
              Additively_weighted_Voronoi_diagram_face_base_2<Gt> > >
class Additively_weighted_Voronoi_diagram_2
  : public Triangulation_2< Gt, Tds >
{
public:
  // type definitions
  typedef Triangulation_2< Gt, Tds >   Delaunay_graph;
  typedef Gt                           Geom_traits;
  typedef typename Gt::Bare_point      Point;
  typedef typename Gt::Weighted_point  Weighted_point;
  typedef typename Gt::Weight          Weight;
  typedef typename Gt::Line_2          Line;

  typedef typename Delaunay_graph::Vertex         Vertex;
  typedef typename Delaunay_graph::Face           Face;

  typedef typename Delaunay_graph::Face_handle    Face_handle;
  typedef typename Delaunay_graph::Vertex_handle  Vertex_handle;
  typedef typename Delaunay_graph::Edge           Edge;
  typedef typename Delaunay_graph::Face_circulator       Face_circulator;
  typedef typename Delaunay_graph::Edge_circulator       Edge_circulator;
  typedef typename Delaunay_graph::All_edges_iterator    All_edges_iterator;
  typedef typename Delaunay_graph::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Delaunay_graph::All_faces_iterator    All_faces_iterator;
  typedef typename Delaunay_graph::Finite_faces_iterator Finite_faces_iterator;
  typedef typename Delaunay_graph::Finite_vertices_iterator 
                                                     Finite_vertices_iterator;

private:
  typedef typename Tds::Vertex_base        Vertex_base;
  // point lists
  typedef std::vector<Weighted_point>      Weighted_point_list;
  typedef Weighted_point_list::iterator    Weighted_point_list_iterator;

  typedef std::map<Face_handle,bool>           Face_map;
  typedef std::map<Face_handle, Face_handle>   Face_face_map;
  typedef std::map<Vertex_handle,bool>         Vertex_map;
  typedef std::vector<Edge>                    Edge_list;

  typedef std::list<Vertex_handle>         Vertex_list;
  typedef Vertex_list::iterator            Vertex_list_iterator;
  typedef Vertex_handle                    Vh_triple[3];

  // the conflict test results
  typedef enum {DO_NOT_FLIP = -1, FLIP}  Conflict_test_result;

private:
  // helper classes
  class Weighted_point_less_than_comparator
  {
  private:
    const Gt& gt;
  public:
    Weighted_point_less_than_comparator(const Gt& gt) : gt(gt) {}

    bool operator ()(const Weighted_point& p,
		     const Weighted_point& q) {
      Comparison_result result =
	gt.awvd_compare_weight_test_2_object()(p, q);
      return (result == LARGER);
    }
  };

  // declare the class Queue as a friend
  friend class Queue;

  // the queue
  class Queue {
  private:
    Edge _front;
    Additively_weighted_Voronoi_diagram_2& awvd_2;

  public:
    bool is_valid() const {}

    inline void push_first(const Edge& e) {
      _front = e;
      set_next(e, e);
      set_prev(e, e);
    }

    inline Edge next(const Edge& e) const {
      CGAL_triangulation_precondition( is_in_queue(e) );
      std::pair<void*, int> _next = e.first->next_in_queue(e.second);
      Face* fptr = static_cast<Face*>(_next.first);
      return Edge(Face_handle(fptr), _next.second);
    }

    inline Edge prev(const Edge& e) const {
      CGAL_triangulation_precondition( is_in_queue(e) );
      std::pair<void*, int> _prev = e.first->prev_in_queue(e.second);
      Face* fptr = static_cast<Face*>(_prev.first);
      return Edge(Face_handle(fptr), _prev.second);
    }

    inline void set_next(const Edge& e, const Edge& next) {
      pair<void*,int> _next(next.first.ptr(), next.second);
      Edge e_sym = awvd_2.sym_edge(e);
      e.first->set_next_in_queue(e.second, _next);
      e_sym.first->set_next_in_queue(e_sym.second, _next);
    }

    inline void set_prev(const Edge& e, const Edge& prev) {
      pair<void*,int> _prev(prev.first.ptr(), prev.second);
      Edge e_sym = awvd_2.sym_edge(e);
      e.first->set_prev_in_queue(e.second, _prev);
      e_sym.first->set_prev_in_queue(e_sym.second, _prev);
    }

    inline bool is_first(const Edge& e) const {
      Edge e_sym = awvd_2.sym_edge(e);
      return ( (e.first == _front.first &&
		e.second == _front.second) ||
	       (e_sym.first == _front.first &&
		e_sym.second == _front.second) );
    }

  public:
    inline bool is_singleton() const {
      CGAL_triangulation_precondition( !is_empty() );
      Edge last_edge = back();
      Edge sym_last_edge = awvd_2.sym_edge(last_edge);
      return ( (_front.first == last_edge.first &&
		_front.second == last_edge.second) ||
	       (_front.first == sym_last_edge.first &&
		_front.second == sym_last_edge.second) );
    }

  public:
    Queue(Additively_weighted_Voronoi_diagram_2& awvd_2,
	  const Edge& e = Edge(Face_handle(NULL),-1) )
      : awvd_2(awvd_2) {
      _front = e;
    }

    inline Edge front() const {
      CGAL_triangulation_precondition( !is_empty() );
      return _front;
    }

    inline Edge back() const {
      CGAL_triangulation_precondition( !is_empty() );
      return prev(_front);
    }

    inline bool is_empty() const {
      return ( _front.first == NULL );
    }

    inline void pop() {
      CGAL_triangulation_precondition( !is_empty() );
      remove(front()); // it is important here that I do not pass the
      // variable _front but rather a copy of it...
    }

    inline void push_front(const Edge& e) {
      CGAL_triangulation_precondition( !is_in_queue(e) );
      push(e);
      _front = e;
    }

    inline void push_back(const Edge& e) {
      push(e);
    }

    void push(const Edge& e) {
      CGAL_triangulation_precondition( !is_in_queue(e) );

      if ( is_empty() ) {
	push_first(e);
	return;
      }
      Edge last_edge = back();
      set_next(last_edge, e);
      set_next(e, _front);
      set_prev(e, last_edge);
      set_prev(_front, e);
    }
    
    void remove(const Edge& e) {
      CGAL_triangulation_precondition( is_in_queue(e) );
      static Edge SENTINEL_QUEUE_EDGE = Edge(Face_handle(NULL), -1);

      if ( is_singleton() ) {
	_front = SENTINEL_QUEUE_EDGE;
	set_next(e, SENTINEL_QUEUE_EDGE);
	set_prev(e, SENTINEL_QUEUE_EDGE);
	return;
      }

      Edge _next = next(e);
      Edge _prev = prev(e);

      if ( is_first(e) ) {
	_front = _next;
      }

      set_next(e, SENTINEL_QUEUE_EDGE);
      set_prev(e, SENTINEL_QUEUE_EDGE);

      set_next(_prev, _next);
      set_prev(_next, _prev);
    }

    bool is_in_queue(const Edge& e) const {
      Edge sym = awvd_2.sym_edge(e);
      CGAL_triangulation_precondition
	( e.first->is_in_queue(e.second) ==
	  sym.first->is_in_queue(sym.second) );

      return e.first->is_in_queue(e.second);
    }

    void clear() {
      while ( !is_empty() ) {
	pop();
      }
    }
  };

  friend class List;

  class List {
    friend class List_iterator;
  private:
    Edge _front;
    Additively_weighted_Voronoi_diagram_2& awvd_2;
    unsigned int _size;

  private:
    inline void increase_size() {
      _size++;
    }

    inline void decrease_size() {
      _size--;
    }

  public:
    bool is_valid() const {}

    inline unsigned int size() const {
      return _size;
    }

    inline void push_first(const Edge& e) {
      _front = e;
      set_next(e, e);
      set_prev(e, e);
      increase_size();
    }

    inline Edge next(const Edge& e) const {
      CGAL_triangulation_precondition( is_in_list(e) );
      std::pair<void*, int> _next = e.first->next_in_queue(e.second);
      Face* fptr = static_cast<Face*>(_next.first);
      return Edge(Face_handle(fptr), _next.second);
    }

    inline Edge prev(const Edge& e) const {
      CGAL_triangulation_precondition( is_in_list(e) );
      std::pair<void*, int> _prev = e.first->prev_in_queue(e.second);
      Face* fptr = static_cast<Face*>(_prev.first);
      return Edge(Face_handle(fptr), _prev.second);
    }

    inline void set_next(const Edge& e, const Edge& next) {
      pair<void*,int> _next(next.first.ptr(), next.second);
      e.first->set_next_in_queue(e.second, _next);
    }

    inline void set_prev(const Edge& e, const Edge& prev) {
      pair<void*,int> _prev(prev.first.ptr(), prev.second);
      e.first->set_prev_in_queue(e.second, _prev);
    }

    inline bool is_first(const Edge& e) const {
      return ( (e.first == _front.first &&
		e.second == _front.second) );
    }

  public:
    inline bool is_singleton() const {
      CGAL_triangulation_precondition( !is_empty() );
      return (size() == 1);
    }

  public:
    List(Additively_weighted_Voronoi_diagram_2& awvd_2,
	 const Edge& e = Edge(Face_handle(NULL),-1) )
      : awvd_2(awvd_2), _size(0) {
      _front = e;
    }

    inline Edge front() const {
      CGAL_triangulation_precondition( !is_empty() );
      return _front;
    }

    inline Edge back() const {
      CGAL_triangulation_precondition( !is_empty() );
      return prev(_front);
    }

    inline bool is_empty() const {
      return ( _front.first == NULL );
    }

    inline void pop() {
      CGAL_triangulation_precondition( !is_empty() );
      remove(front()); // it is important here that I do not pass the
      // variable _front but rather a copy of it...
    }

    inline void push_front(const Edge& e) {
      CGAL_triangulation_precondition( !is_in_list(e) );
      push(e);
      _front = e;
    }

    inline void push_back(const Edge& e) {
      push(e);
    }

    void push(const Edge& e) {
      CGAL_triangulation_precondition( !is_in_list(e) );

      if ( is_empty() ) {
	push_first(e);
	return;
      }
      Edge last_edge = back();
      set_next(last_edge, e);
      set_next(e, _front);
      set_prev(e, last_edge);
      set_prev(_front, e);

      increase_size();
    }

    inline void insert_after(const Edge& e, const Edge& new_e) {
      CGAL_triangulation_precondition( is_in_list(e) );
      Edge old_front = _front;
      _front = next(e);
      push_front(new_e);
      _front = old_front;
    }

    inline void insert_before(const Edge& e, const Edge& new_e) {
      CGAL_triangulation_precondition( is_in_list(e) );
      Edge old_front = _front;
      _front = e;
      push(new_e);
      _front = old_front;
    }

    inline void replace(const Edge& e, const Edge& new_e) {
      insert_before(e, new_e);
      remove(e);
    }

    void remove(const Edge& e) {
      CGAL_triangulation_precondition( is_in_list(e) );
      static Edge SENTINEL_QUEUE_EDGE = Edge(Face_handle(NULL), -1);

      if ( is_singleton() ) {
	_front = SENTINEL_QUEUE_EDGE;
	set_next(e, SENTINEL_QUEUE_EDGE);
	set_prev(e, SENTINEL_QUEUE_EDGE);
	return;
      }

      Edge _next = next(e);
      Edge _prev = prev(e);

      if ( is_first(e) ) {
	_front = _next;
      }

      set_next(e, SENTINEL_QUEUE_EDGE);
      set_prev(e, SENTINEL_QUEUE_EDGE);

      set_next(_prev, _next);
      set_prev(_next, _prev);
      decrease_size();
    }

    bool is_in_list(const Edge& e) const {
      return e.first->is_in_queue(e.second);
    }

    void clear() {
      while ( !is_empty() ) {
	pop();
      }
    }
  };


  class Vertex_triple {
  private:
    Vertex* _v[3];

    void init(Vertex* v1, Vertex* v2, Vertex* v3) {
      if ( v1 <= v2 && v1 <= v3 ) {
	_v[0] = v1;
	_v[1] = v2;
	_v[2] = v3;
      } else if ( v2 <= v1 && v2 <= v3 ) {
	_v[0] = v2;
	_v[1] = v3;
	_v[2] = v1;
      } else {
	_v[0] = v3;
	_v[1] = v1;
	_v[2] = v2;
      }
    }

  public:
    Vertex_triple() { _v[0] = _v[1] = _v[2] = NULL; }
    Vertex_triple(const Face_handle& f) {
      init(f->vertex(0).ptr(),f->vertex(1).ptr(),f->vertex(2).ptr());
    }
    Vertex_triple(const Face* f) {
      init(f->vertex(0).ptr(),f->vertex(1).ptr(),f->vertex(2).ptr());
    }
    Vertex_triple(const Vertex_handle& v1, const Vertex_handle& v2,
		  const Vertex_handle& v3) {
      init(v1.ptr(), v2.ptr(), v3.ptr());
    }
    Vertex_triple(Vertex* v1, Vertex* v2, Vertex* v3) {
      init(v1, v2, v3);
    }

    inline bool operator==(const Vertex_triple& vt) const {
      return (_v[0] == vt._v[0] && _v[1] == vt._v[1] &&_v[2] == vt._v[2]);
    }

    inline bool operator<(const Vertex_triple& vt) const {
      for (int i = 0; i < 3; i++) {
	if ( _v[i] > vt._v[i] ) return false;
	if ( _v[i] < vt._v[i] ) return true;
      }
      return false;
    }

    Vertex* operator[](int i) { return _v[i]; }

  };



public:
  Additively_weighted_Voronoi_diagram_2(const Gt& gt=Gt()) :
    Delaunay_graph(gt) {}

  template< class Input_iterator >
  Additively_weighted_Voronoi_diagram_2(Input_iterator first,
				       Input_iterator beyond,
				       const Gt& gt=Gt()) :
    Delaunay_graph(gt)
  {
    Input_iterator it;
    Weighted_point_less_than_comparator less_than(gt);
    Weighted_point_list wp_list;
    for (it = first; it != beyond; ++it)  wp_list.push_back(*it);
    std::sort(wp_list.begin(), wp_list.end(), less_than);
    Weighted_point_list_iterator lit;
    for (lit = wp_list.begin(); lit != wp_list.end(); ++lit) {
      insert(*lit, true);
    }
  }

  Additively_weighted_Voronoi_diagram_2
  (const Additively_weighted_Voronoi_diagram_2 &awvd)
    : Delaunay_graph(awvd)
  {
    CGAL_triangulation_postcondition( is_valid() );
  }

  Additively_weighted_Voronoi_diagram_2
  operator=(const Additively_weighted_Voronoi_diagram_2& awvd)
  {
    Delaunay_graph::operator=(awvd);
    return (*this);
  }

public:
  // predicates
  Sign
  awvd_distance2_test(const Weighted_point &p,
		      const Weighted_point &q) const;

  Sign awvd_test(const Weighted_point &p, const Weighted_point &q,
		 const Weighted_point &r, const Weighted_point &s) const;

  Sign awvd_test_degenerated(const Weighted_point &p,
			     const Weighted_point &q,
			     const Weighted_point &r) const;

  Sign awvd_test(const Face_handle &f, 
		 const Weighted_point &p) const;

  Sign awvd_test(const Face_handle &f, int i,
		 const Weighted_point &p) const;

  Sign awvd_Voronoi_radii_difference_test(const Weighted_point &p,
					  const Weighted_point &q,
					  const Weighted_point &r, 
					  const Weighted_point &s) const;

  Sign awvd_Voronoi_radii_difference_test(const Face_handle& f, int i,
					  const Weighted_point &p) const;

  Orientation awvd_orientation_test(const Weighted_point& p1,
				    const Weighted_point& p2,
				    const Weighted_point& p3,
				    const Weighted_point& q1,
				    const Weighted_point& q2) const;

  Orientation awvd_orientation_test(const Face_handle& f,
				    const Weighted_point& p,
				    const Weighted_point& q) const;

  Orientation awvd_orientation_left_test(const Weighted_point& p1,
					 const Weighted_point& p2,
					 const Point& q1,
					 const Point& q2) const;

  Orientation awvd_orientation_right_test(const Weighted_point& p1,
					  const Weighted_point& p2,
					  const Point& q1,
					  const Point& q2) const;

  Bounded_side
  awvd_circle_existence_test(const Weighted_point &p,
			     const Weighted_point &q,
			     const Weighted_point &r) const;

  Bounded_side
  awvd_circle_existence_test(const Face_handle& f) const;

  typedef enum { FIRST = 1, SECOND, THIRD } Order;
  Order awvd_ordering_test(const Weighted_point &p,
			   const Weighted_point &q,
			   const Weighted_point &r) const;

  Comparison_result awvd_compare_distances_test(const Weighted_point &p1,
						const Weighted_point &p2,
						const Point &p) const;

  Sign
  awvd_vertex_test(const Weighted_point &p1, const Weighted_point &p2,
		   const Weighted_point &p3, const Weighted_point &q) const;
  Sign
  awvd_vertex_test_degenerated(const Weighted_point &p1,
			       const Weighted_point &p2,
			       const Weighted_point &q) const;

  Sign
  awvd_vertex_test(const Face_handle& f, const Weighted_point& p) const;

  //****************************************************************

  Comparison_result
  awvd_compare_weights(const Weighted_point& p,
		       const Weighted_point& q) const;

  Conflict_test_result
  conflict_test_for_Voronoi_center(const Face_handle& f,
				   const Weighted_point& p) const;

  Edge_conflict_test_result
  awvd_finite_edge_test(const Weighted_point& p1,
			const Weighted_point& p2,
			const Weighted_point& p3,
			const Weighted_point& p4,
			const Weighted_point& q) const;

  Edge_conflict_test_result
  awvd_finite_edge_test(const Face_handle& f, int i,
			const Weighted_point& p) const;

  Edge_conflict_test_result
  awvd_finite_edge_test_degenerated(const Weighted_point& p1,
				    const Weighted_point& p2,
				    const Weighted_point& p3,
				    const Weighted_point& q) const;


  Edge_conflict_test_result
  awvd_finite_edge_test_degenerated(const Weighted_point& p1,
				    const Weighted_point& p2,
				    const Weighted_point& q) const;

  Edge_conflict_test_result
  awvd_finite_edge_test_degenerated(const Face_handle& f, int i,
				    const Weighted_point& p) const;

  Edge_conflict_test_result
  awvd_infinite_edge_test(const Weighted_point& p2,
			  const Weighted_point& p3,
			  const Weighted_point& p4,
			  const Weighted_point& q) const;


  Edge_conflict_test_result
  awvd_infinite_edge_test(const Face_handle& f, int i,
			  const Weighted_point& p) const;

private:
  Edge_conflict_test_result
  edge_conflict_test(const Face_handle& f, int i,
		     const Weighted_point& p) const;

  Edge_conflict_test_result
  edge_conflict_test(const Face_handle& f, int i,
		     const Vertex_handle& v4,
		     const Weighted_point& p) const;

  inline
  Edge_conflict_test_result
  edge_conflict_test(const Edge& e,
		     const Weighted_point& p) const
  {
    return edge_conflict_test(e.first, e.second, p);
  }

public:
  // Primal
  Weighted_point primal(const Face_handle& f) const;
  Object primal(const Edge e) const;
  Object primal(const Edge_circulator& ec) const;
  Object primal(const Finite_edges_iterator& ei) const;

  // Dual
  Point  dual(const Face_handle& f) const;
  Object dual(const Edge e) const;
  Object dual(const Edge_circulator& ec) const;
  Object dual(const Finite_edges_iterator& ei) const;

  // constructions
  Point awvd_circumcenter(const Face_handle& f) const;
  Point awvd_circumcenter(const Weighted_point& p0, 
			  const Weighted_point& p1, 
			  const Weighted_point& p2) const;

  Weighted_point awvd_circumcircle(const Face_handle& f) const;
  Weighted_point awvd_circumcircle(const Weighted_point& p0, 
				   const Weighted_point& p1, 
				   const Weighted_point& p2) const;
  Line awvd_circumcircle(const Weighted_point& p0,
			 const Weighted_point& p1) const;

private:
  // combinatorial operations
  Edge flip(Face_handle& f, int i);
  Edge flip(Edge e);

  Vertex_handle insert_in_face(Face_handle& f, const Weighted_point& p);

  bool          is_degree_2(const Vertex_handle& v) const;

  Vertex_handle insert_degree_2(Edge e);
  Vertex_handle insert_degree_2(Edge e, const Weighted_point& p);
  void          remove_degree_2(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v, Face* f);

public:
  // insertion of weighted point

  template< class Input_iterator >
  inline
  void insert(Input_iterator first, Input_iterator beyond)
  {
    Input_iterator it;
    Weighted_point_less_than_comparator less_than(geom_traits());
    Weighted_point_list wp_list;
    for (it = first; it != beyond; ++it)  wp_list.push_back(*it);
    std::sort(wp_list.begin(), wp_list.end(), less_than);
    Weighted_point_list_iterator lit;
    for (lit = wp_list.begin(); lit != wp_list.end(); ++lit) {
      insert(*lit);
    }
    wp_list.clear();
  }

  Vertex_handle  insert(const Weighted_point& p);
  Vertex_handle  insert(const Weighted_point& p, Vertex_handle vnearest);

private:
  Vertex_handle  insert_first(const Weighted_point& p);
  Vertex_handle  insert_second(const Weighted_point& p);
  Vertex_handle  insert_third(const Weighted_point& p,
			      std::vector<Vh_triple*>* flipped_edges);

  Vertex_handle  insert(const Weighted_point& p, Vertex_handle vnearest,
			std::vector<Vh_triple*>* flipped_edges);



public:
  inline
  std::vector<Edge>
  find_conflict_region(const Weighted_point& p)
    {
      CGAL_precondition ( number_of_vertices() > 2 );
      Vertex_handle vnearest = find_nearest_neighbor(p);
      return find_conflict_region(p, vnearest);
    }

  std::vector<Edge>
  find_conflict_region(const Weighted_point& p,	Vertex_handle vnearest);

private:
  inline
  std::vector<Edge> find_conflict_region(const Vertex_handle& vnearest,
					 const Weighted_point& p,
					 Vertex_handle& v,
					 bool insert = false)
  {
    
    return find_conflict_region(vnearest, p, NULL, v, insert);
  }

  void find_conflict_region(const Weighted_point& p,
			    const Vertex_handle& vnearest,
			    List& l, Face_map& fm,
			    Vertex_map& vm,
			    std::vector<Vh_triple*>* fe);

  void initialize_conflict_region(const Face_handle& f, List& l);
  void check_edge_for_trivial_weighted_points(const Face_handle& f, int i,
					      const Weighted_point& p,
					      Vertex_map& vm);
  void expand_conflict_region(const Face_handle& f, const Weighted_point& p,
			      List& l, Face_map& fm, Vertex_map& vm,
			      std::vector<Vh_triple*>* fe);

private:
  Vertex_handle add_bogus_vertex(Edge e, List& l);
  Vertex_list   add_bogus_vertices(List& l);
  void          remove_bogus_vertices(Vertex_list& vl);

  void move_trivial_weighted_points(Vertex_handle& vold,
				    Vertex_handle& vnew);

  inline  std::vector<Face*> get_faces_for_recycling(Face_map& fm,
					     unsigned int n_wanted);
  void remove_trivial_vertices(Vertex_handle&v, Vertex_map& vm,
			       Face_map& fm);
  void retriangulate_conflict_region(Vertex_handle &v, List& l,
				     Face_map& fm, Vertex_map& vm);
public:
  // removal of weighted point
  void           remove(Vertex_handle v);
  void           remove_using_conflict_region(Vertex_handle v);

private:
  bool  is_degree_2_vertex(const Vertex_handle& v) const;
  void  remove_first(Vertex_handle v);
  void  remove_second(Vertex_handle v);
  void  remove_third(Vertex_handle v);
  void  remove_degree_d_vertex(Vertex_handle v);
  void  remove_degree_d_vertex_using_conflict_region(Vertex_handle v);
  void  remove_degree_d_vertex_using_conflict_region2(Vertex_handle v);

public:
  // point location methods
  inline
  std::list<Vertex_handle>
  path_to_nearest_neighbor(const Weighted_point& p) const
    {
      std::list<Vertex_handle> vertex_list;
      std::back_insert_iterator<std::list<Vertex_handle> >  it(vertex_list);
      path_to_nearest_neighbor(p, it);
      return vertex_list;
    }

  Vertex_handle  find_nearest_neighbor(const Point& p) const;
  Vertex_handle  find_nearest_neighbor(const Point& p,
				       Vertex_handle start_v) const;

private:
  template< class Back_insert_iterator>
  void
  path_to_nearest_neighbor(const Weighted_point& p,
			   Back_insert_iterator& it) const
    {
      if ( number_of_vertices() == 0 ) { return; }

      Vertex_handle start_vertex;
      if ( start_vertex.ptr() == NULL ) {
	Vertex_iterator vit = vertices_begin();
	for (; vit != vertices_end(); ++vit) {
	  if ( !is_infinite(vit) ) {
	    start_vertex = Vertex_handle(vit);
	    break;
	  }
	}
      }

      Vertex_handle vclosest;
      Vertex_handle v = start_vertex;

      if ( number_of_vertices() < 3 ) {
	vclosest = v;
	*it++ = vclosest;
	Vertex_iterator vit = vertices_begin();
	for (; vit != vertices_end(); ++vit) {
	  Vertex_handle v1(vit);
	  if ( v1 != vclosest && !is_infinite(v1) ) {
	    Weighted_point p1 = vclosest->point();
	    Weighted_point p2 = v1->point();
	    if ( awvd_compare_distances_test(p1, p2, p) == LARGER ) {
	      vclosest = v1;
	      *it++ = vclosest;
	    }
	  }
	}
	return;
      }

      vclosest = v;
      *it++ = vclosest;
      Weighted_point p1 = v->point();
      Face_circulator fc(v);
      Face_circulator fc_start = fc;
      Face_handle ff;
      do {
	ff = Face_handle(fc);
	Vertex_handle v1 = ff->vertex( ccw(ff->index(v)) );
	if ( !is_infinite(v1) ) {
	  Weighted_point p2 = v1->point();
	  if ( awvd_compare_distances_test(p1, p2, p) == LARGER ) {
	    v = v1;
	    break;
	  }
	}
	++fc;
      } while ( fc != fc_start );
      if ( vclosest == v ) { return; }

      do {
	vclosest = v;
	*it++ = vclosest;
	Weighted_point p1 = v->point();
	fc = Face_circulator(v, ff);
	fc_start = fc;
	do {
	  ff = Face_handle(fc);
	  Vertex_handle v1 = ff->vertex( ccw(ff->index(v)) );
	  if ( !is_infinite(v1) ) {
	    Weighted_point p2 = v1->point();
	    if ( awvd_compare_distances_test(p1, p2, p) == LARGER ) {
	      v = v1;
	      break;
	    }
	  }
	  ++fc;
	} while ( fc != fc_start );
      } while ( vclosest != v );
    }


public:
  // validity check
  bool is_valid(bool verbose = false, int level = 0) const;

private:
  // some access methods
  inline Edge sym_edge(const Edge e) const
  {
    return sym_edge(e.first, e.second);
  }

  inline Edge sym_edge(const Face_handle& f, int i) const
  {
    Face_handle f_sym = f->neighbor(i);
    return Edge(  f_sym, f_sym->index( f->mirror_vertex(i) )  );
  }


public:
  // I/O via streams

  template < class Stream >
  Stream& write_non_trivial_weighted_points(Stream& str) const
  {
    str << number_of_vertices() << std::endl;

    Vertex_iterator vit = vertices_begin();
    for (; vit != vertices_end(); ++vit) {
      Weighted_point wp = vit->point();
      str << wp.x() << " " << wp.y() << " " << wp.weight() << std::endl;
    }
    return str;
  }

  template < class Stream >
  Stream& draw_trivial_weighted_points(Stream &str) const
  {
    Vertex_iterator vit = vertices_begin();

    for (; vit != vertices_end(); ++vit) {
      if ( vit->number_of_weighted_points() > 0 ) {
	typename Vertex_base::Weighted_point_list_iterator wpit;
	for (wpit = vit->weighted_points_begin();
	     wpit != vit->weighted_points_end(); ++wpit) {
	  Weighted_point wp = *wpit;
	  Circle_2< typename Gt::Rep > c(wp.point(),
					 CGAL_NTS square(wp.weight()));
	  //	  str << wp.point();
	  str << c;
	}
      }
    }
    return str;
  }

  template < class Stream >
  Stream& draw_non_trivial_weighted_points(Stream &str) const
  {
    Vertex_iterator vit = vertices_begin();
    for (; vit != vertices_end(); ++vit) {
      Weighted_point wp = vit->point();
      Circle_2< typename Gt::Rep > c(wp.point(),
				     CGAL_NTS square(wp.weight()));
      //      str << wp.point();
      str << c;
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
	  str << awvd_circumcircle( f->vertex(1)->point(),
				    f->vertex(2)->point() );
	} else if (  is_infinite(f->vertex(1))  ){
	  str << awvd_circumcircle( f->vertex(2)->point(),
				    f->vertex(0)->point() );
	} else {
	  str << awvd_circumcircle( f->vertex(0)->point(),
				    f->vertex(1)->point() );	  
	}
      } else {
	Weighted_point wp = awvd_circumcircle(f);
	Circle_2< typename Gt::Rep > c(wp.point(),
				       CGAL_NTS square(wp.weight()));
	str << c;
      }
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
      Vertex_handle v1(vertices_begin());
      Vertex_handle v2(++vertices_begin());
      str << Line(v1->point(), v2->point());
    } else if ( number_of_vertices() < 3 ) {
      Finite_edges_iterator eit = finite_edges_begin();
      for (; eit != finite_edges_end(); ++eit) {
	draw_edge< Stream >(*eit, str);
      }
    } else {
      All_edges_iterator eit = all_edges_begin();
      for (; eit != all_edges_end(); ++eit) {
	draw_edge< Stream >(*eit, str);
      }
    }
    return str;
  }

  template< class Stream >
  Stream& draw_edge(Edge e, Stream &str) const
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
  Stream& draw_dual(Edge e, Stream &str) const
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
      draw_edge< Stream >(Edge(f, i), str);
    }
    return str;
  }

}; // Additively_weighted_Voronoi_diagram_2


//--------------------------------------------------------------------
// test method
//--------------------------------------------------------------------
template < class Gt, class Tds >
bool
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
is_valid(bool verbose, int level) const
{
  if (number_of_vertices() <= 1)    return true; 

  bool result = _tds.is_valid(verbose, level);

  CGAL_assertion( result );

  if (level <= 0) return true;

  if (number_of_vertices() < 3)  return true;

  Weighted_point p, q, r, s;

  for (All_faces_iterator fit = all_faces_begin();
       fit != all_faces_end(); ++fit) {
    Face_handle f(fit);
    if ( !is_infinite(f) ) {
      p = f->vertex(0)->point();
      q = f->vertex(1)->point();
      r = f->vertex(2)->point();
      result = result &&
	(awvd_circle_existence_test(p, q, r) != ON_BOUNDED_SIDE);
    }
  }

  CGAL_triangulation_assertion(result);

  for (All_edges_iterator eit = all_edges_begin(); 
       eit != all_edges_end(); ++eit) {
    Edge e = *eit;
    Vertex_handle vp = e.first->vertex( ccw(e.second) );
    Vertex_handle vq = e.first->vertex(  cw(e.second) );
    Vertex_handle vr = e.first->vertex(   e.second    );
    Vertex_handle vs = e.first->mirror_vertex(e.second);

    if ( is_infinite(vp) ) {
      q = vq->point();
      r = vr->point();
      s = vs->point();
      result = result &&
	(awvd_test_degenerated(s, q, r) != NEGATIVE) &&
	(awvd_test_degenerated(q, r, s) != NEGATIVE);
      CGAL_triangulation_assertion(result);
      continue;
    }
    if ( is_infinite(vq) ) {
      p = vp->point();
      r = vr->point();
      s = vs->point();
      result = result &&
      	(awvd_test_degenerated(r, p, s) != NEGATIVE) &&
	(awvd_test_degenerated(p, s, r) != NEGATIVE);
      CGAL_triangulation_assertion(result);
      continue;
    }
    if ( is_infinite(vr) ) {
      if ( !is_infinite(vs) ) {
	p = vp->point();
	q = vq->point();
	s = vs->point();
	result = result &&
	  (awvd_test_degenerated(p, q, s) != NEGATIVE);
	CGAL_triangulation_assertion(result);
      }
      continue;
    }
    if ( is_infinite(vs) ) {
      if ( !is_infinite(vr) ) {
	p = vp->point();
	q = vq->point();
	r = vr->point();
	result = result &&
	  (awvd_test_degenerated(q, p, r) != NEGATIVE);
	CGAL_triangulation_assertion(result);
      }
      continue;
    }
    CGAL_triangulation_assertion( !is_infinite(vp) &&
				  !is_infinite(vq) &&
				  !is_infinite(vr) &&
				  !is_infinite(vs) );
    p = vp->point();
    q = vq->point();
    r = vr->point();
    s = vs->point();

    result = result &&
      (awvd_test(p, q, r, s) != NEGATIVE) &&
      (awvd_test(q, p, s, r) != NEGATIVE);
    CGAL_triangulation_assertion(result);
  }

  if ( level <= 1 ) return true;

  for (All_faces_iterator fit = all_faces_begin();
       fit != all_faces_end(); ++fit) {
    Face_handle f(fit);
    if ( !is_infinite(f) ) {
      p = f->vertex(0)->point();
      q = f->vertex(1)->point();
      r = f->vertex(2)->point();
      for (Vertex_iterator vit = vertices_begin();
	   vit != vertices_end(); ++vit) {
	Vertex_handle vh(vit);
	if ( !is_infinite(vh) ) {
	  s = vh->point();
	  result = result &&
	    (awvd_test(p, q, r, s) != NEGATIVE);
	  CGAL_triangulation_assertion(result);
	}
      }
    } else { // f is infinite
      for (int i = 0; i < 3; i++) {
	if ( is_infinite(f->vertex(i)) ) {
	  p = f->vertex( ccw(i) )->point();
	  q = f->vertex(  cw(i) )->point();
	  break;
	}
      } // end for loop
      for (Vertex_iterator vit = vertices_begin();
	       vit != vertices_end(); ++vit) {
	Vertex_handle vh(vit);
	if ( !is_infinite(vh) ) {
	  r = vh->point();
	  result = result &&
	    (awvd_test_degenerated(p, q, r) != NEGATIVE);
	  CGAL_triangulation_assertion(result);
	} // endif
      } // end for loop
    } // end outer if
  } // and for all faces

  //**************************************************************
  // I ALSO WANTTO ADD A TEST FOR THE TRIVIAL CIRCLES WHICH WILL
  // VERIFY THAT THE TRIVIAL CIRCLES ARE INDEED TRIVIAL

  CGAL_triangulation_assertion(result);
  return result;
}



//--------------------------------------------------------------------
// embedding and visualization methods and constructions for primal
// and dual
//--------------------------------------------------------------------

// circumcenter
template < class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Point
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_circumcenter(const Face_handle& f) const
{
  CGAL_triangulation_precondition (dimension()==2 || !is_infinite(f));
  return awvd_circumcenter(f->vertex(0)->point(),
			   f->vertex(1)->point(),
			   f->vertex(2)->point());
}

template < class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Point
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_circumcenter(const Weighted_point& p0, 
		  const Weighted_point& p1, 
		  const Weighted_point& p2) const
{
  return geom_traits().construct_awvd_circumcenter_2_object()(p0, p1, p2);
}

// circumcircle
template < class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Weighted_point
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_circumcircle(const Face_handle& f) const
{
  CGAL_triangulation_precondition (dimension()==2 || !is_infinite(f));
  return awvd_circumcircle(f->vertex(0)->point(),
			   f->vertex(1)->point(),
			   f->vertex(2)->point());
}

template < class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Weighted_point
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_circumcircle(const Weighted_point& p0, 
		  const Weighted_point& p1, 
		  const Weighted_point& p2) const
{
  return geom_traits().construct_awvd_circumcircle_2_object()(p0, p1, p2);
}


template < class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Line
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_circumcircle(const Weighted_point& p0,
		  const Weighted_point& p1) const
{
  return geom_traits().construct_awvd_left_bitangent_line_2_object()
    (p0, p1);
}


// dual
template < class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Point
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
dual(const Face_handle& f) const
{
  return awvd_circumcenter(f);
}


template < class Gt, class Tds >
inline
Object
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
dual(const Edge e) const
{
  typedef typename Geom_traits::Line_2     Line;
  typedef typename Geom_traits::Segment_2  Segment;
  typedef typename Geom_traits::Ray_2      Ray;
  typedef typename Geom_traits::Hyperbola_2          Hyperbola;
  typedef typename Geom_traits::Hyperbola_segment_2  Hyperbola_segment;
  typedef typename Geom_traits::Hyperbola_ray_2      Hyperbola_ray;

  CGAL_triangulation_precondition( !is_infinite(e) );

  if ( dimension() == 1 ) {
    Weighted_point p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point q = (e.first)->vertex(ccw(e.second))->point();
    if ( CGAL_NTS compare(p.weight(), q.weight()) == EQUAL ) {
      Line l1 = geom_traits().construct_awvd_line_2_object()(p.point(),
							     q.point());
      //      Line l1(p1, p2);
      Point m((p.x() + q.x()) / 2, (p.y() + q.y()) / 2);
      Point m1(m.x() + l1.a(), m.y() + l1.b());
      Line l = geom_traits().construct_awvd_line_2_object()(m, m1);
      return CGAL::make_object(l);
    }
    Hyperbola h = geom_traits().construct_awvd_hyperbola_2_object()(p, q);
    return CGAL::Object(new CGAL::Wrapper< Hyperbola >(h));
  }

  // dimension == 2
  // noe of the two adjacent faces is infinite
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    Weighted_point p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point q = (e.first)->vertex(ccw(e.second))->point();
    if ( CGAL_NTS compare(p.weight(), q.weight()) == EQUAL ) {
      Segment s = geom_traits().construct_awvd_segment_2_object()
	(dual(e.first),dual(e.first->neighbor(e.second)));
      return CGAL::Object(new CGAL::Wrapper< Segment >(s));
    }
    Hyperbola_segment hs =
      geom_traits().construct_awvd_hyperbola_segment_2_object()
      (p, q, dual(e.first),dual(e.first->neighbor(e.second)));
    return CGAL::Object(new CGAL::Wrapper< Hyperbola_segment >(hs));
  }

  // both of the adjacent faces are infinite
  if ( is_infinite(e.first) &&
       is_infinite(e.first->neighbor(e.second)) )  {
    Weighted_point p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point q = (e.first)->vertex(ccw(e.second))->point();
    if ( CGAL_NTS compare(p.weight(), q.weight()) == EQUAL ) {
      Line l = geom_traits().construct_awvd_line_2_object()(p.point(),
							    q.point());
      return CGAL::Object(new CGAL::Wrapper< Line >(l));
    }
    Hyperbola h = geom_traits().construct_awvd_hyperbola_2_object()(p, q);
    return CGAL::Object(new CGAL::Wrapper< Hyperbola >(h));
  }

  // only one of the adjacent faces is infinite
  CGAL_triangulation_assertion( is_infinite( e.first ) ||
				is_infinite( e.first->neighbor(e.second) )
				);

  CGAL_triangulation_assertion( !(is_infinite( e.first ) &&
				  is_infinite( e.first->neighbor(e.second) )
				  )
				);

  CGAL_triangulation_assertion
    (  is_infinite( e.first->vertex(e.second) ) ||
       is_infinite( e.first->mirror_vertex(e.second) )  );

  Face_handle f;
  Sign direction;
  if (  is_infinite( e.first->vertex(e.second) )  ) {
    f = e.first->neighbor(e.second);
    direction = POSITIVE;
  } else {
    f = e.first;
    direction = NEGATIVE;
  }

  Weighted_point p = e.first->vertex( ccw(e.second) )->point();
  Weighted_point q = e.first->vertex( cw(e.second) )->point();
  //**********************************************************
  // here i have to add the ray_2 case...
  Hyperbola_ray hr =
    geom_traits().construct_awvd_hyperbola_ray_2_object()
    (p, q, dual(f), direction);
  return CGAL::Object(new CGAL::Wrapper< Hyperbola_ray >(hr));
}


template< class Gt, class Tds >
inline
Object
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
dual(const Edge_circulator& ec) const
{
  return dual(*ec);
}


template< class Gt, class Tds >
inline
Object
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
dual(const Finite_edges_iterator& ei) const
{
  return dual(*ei);
}


// primal
template < class Gt, class Tds >
inline
Object
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
primal(const Edge e) const
{
  typedef typename Geom_traits::Point::RT  RT;
  typedef typename Geom_traits::Segment_2  Segment;
  typedef typename Geom_traits::Ray_2      Ray;
  typedef typename Geom_traits::Hyperbola_segment_2  Hyperbola_segment;
  typedef typename Geom_traits::Parabola_segment_2   Parabola_segment;

  //  CGAL_triangulation_precondition( !is_infinite(e) );

  if ( number_of_vertices() != 2 ) {
    if ( is_infinite(e) ) {
      if (  is_infinite( e.first->vertex(cw(e.second)) )  ) {
	Weighted_point p = e.first->vertex( ccw(e.second) )->point();
	Weighted_point r = e.first->vertex( e.second )->point();
	Weighted_point s = e.first->mirror_vertex( e.second )->point();
	Line l1 = awvd_circumcircle(r, p);
	Line l2 = awvd_circumcircle(p, s);
	RT d1 = CGAL_NTS sqrt( CGAL_NTS square(l1.a()) +
			       CGAL_NTS square(l1.b()) );
	RT d2 = CGAL_NTS sqrt( CGAL_NTS square(l2.a()) +
			       CGAL_NTS square(l2.b()) );
	RT a = l1.a() / d1 - l2.a() / d2;
	RT b = l1.b() / d1 - l2.b() / d2;
	Point c(p.x() + b, p.y() - a);
	Ray ray(p, c);
	return CGAL::Object(new CGAL::Wrapper< Ray >(ray));
      } else {
	CGAL_triangulation_assertion
	  (   is_infinite( e.first->vertex(ccw(e.second)) )   );
	Weighted_point q = e.first->vertex( cw(e.second) )->point();
	Weighted_point r = e.first->vertex( e.second )->point();
	Weighted_point s = e.first->mirror_vertex( e.second )->point();
	Line l1 = awvd_circumcircle(s, q);
	Line l2 = awvd_circumcircle(q, r);
	RT d1 = CGAL_NTS sqrt( CGAL_NTS square(l1.a()) +
			       CGAL_NTS square(l1.b()) );
	RT d2 = CGAL_NTS sqrt( CGAL_NTS square(l2.a()) +
			       CGAL_NTS square(l2.b()) );
	RT a = l1.a() / d1 - l2.a() / d2;
	RT b = l1.b() / d1 - l2.b() / d2;
	Point c(q.x() + b, q.y() - a);
	Ray ray(q, c);
	return CGAL::Object(new CGAL::Wrapper< Ray >(ray));
      }
    }
  }

  if ( dimension() == 1 ) {
    Weighted_point p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point q = (e.first)->vertex(ccw(e.second))->point();
    Segment s = geom_traits().construct_awvd_segment_2_object()
      (p.point(), q.point());
    return CGAL::Object(new CGAL::Wrapper< Segment >(s));
  }

  // dimension == 2
  if( (!is_infinite(e.first)) &&
      (!is_infinite(e.first->neighbor(e.second))) ) {
    Weighted_point p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point q = (e.first)->vertex(ccw(e.second))->point();
    Weighted_point d1 = awvd_circumcircle(e.first);
    Weighted_point d2 = awvd_circumcircle(e.first->neighbor(e.second));
    if ( CGAL_NTS compare(d1.weight(), d2.weight()) == EQUAL ) {
      Segment s = geom_traits().construct_awvd_segment_2_object()
	(p.point(), q.point());
      return CGAL::Object(new CGAL::Wrapper< Segment >(s));
    }
    Hyperbola_segment hs =
      geom_traits().construct_awvd_hyperbola_segment_2_object()
      (d1, d2, p.point(), q.point());
    return CGAL::Object(new CGAL::Wrapper< Hyperbola_segment >(hs));
  }

  // both of the adjacent faces are infinite
  if ( is_infinite(e.first) &&
       is_infinite(e.first->neighbor(e.second)) )  {
    Weighted_point p = (e.first)->vertex(cw(e.second))->point();
    Weighted_point q = (e.first)->vertex(ccw(e.second))->point();
    Segment s = geom_traits().construct_awvd_segment_2_object()
      (p.point(), q.point());
    return CGAL::Object(new CGAL::Wrapper< Segment >(s));
  }

  // only one of the adjacent faces is infinite
  Weighted_point p = (e.first)->vertex(ccw(e.second))->point();
  Weighted_point q = (e.first)->vertex(cw(e.second))->point();
  if ( is_infinite(e.first) ) {
    Weighted_point c = awvd_circumcircle(e.first->neighbor(e.second));
    Line l = awvd_circumcircle(p, q);
    Parabola_segment ps =
      geom_traits().construct_awvd_parabola_segment_2_object()
      (c, l, p.point(), q.point());
    return CGAL::Object(new CGAL::Wrapper< Parabola_segment >(ps));
  }
  Weighted_point c = awvd_circumcircle(e.first);
  Line l = awvd_circumcircle(q, p);
  Parabola_segment ps =
    geom_traits().construct_awvd_parabola_segment_2_object()
    (c, l, p.point(), q.point());
  return CGAL::Object(new CGAL::Wrapper< Parabola_segment >(ps));
}


template< class Gt, class Tds >
inline
Object
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
primal(const Edge_circulator& ec) const
{
  return primal(*ec);
}


template< class Gt, class Tds >
inline
Object
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
primal(const Finite_edges_iterator& ei) const
{
  return primal(*ei);
}

//--------------------------------------------------------------------
// combinatorial operations
//--------------------------------------------------------------------
template < class Gt, class Tds >
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Edge
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
flip(Face_handle& f, int i)
{
  CGAL_triangulation_precondition ( ! f.is_null() );
  CGAL_triangulation_precondition (i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition( dimension()==2 ); 

  CGAL_triangulation_precondition( f->vertex(i) != f->mirror_vertex(i) );

  //  CGAL_triangulation_assertion( is_valid() );

  _tds.flip( &(*f), i);

  //  CGAL_triangulation_assertion( is_valid() );

  return Edge(f, ccw(i));
}

template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Edge
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
flip(Edge e)
{
  return flip(e.first, e.second);
}

template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
insert_in_face(Face_handle& f, const Weighted_point& p)
{
  Vertex_handle v = static_cast<Vertex*>(_tds.insert_in_face( &(*f) ));

  v->set_point(p);
  return v;
}

template< class Gt, class Tds >
inline
bool
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
is_degree_2(const Vertex_handle& v) const
{
  Face_circulator fc = v->incident_faces();
  Face_circulator fc1 = fc;
  ++(++fc1);
  return ( fc == fc1 );
}

template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
insert_degree_2(Edge e)
{
  return
    static_cast<Vertex*>(_tds.insert_degree_2(e.first.ptr(),e.second));
}

template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
insert_degree_2(Edge e, const Weighted_point& p)
{
  Vertex_handle v = insert_degree_2(e);

  v->set_point(p);
  return v;
}


template< class Gt, class Tds >
inline
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_degree_2(Vertex_handle v)
{
  CGAL_triangulation_precondition( is_degree_2(v) );

  _tds.remove_degree_2(v.ptr());
}


template< class Gt, class Tds >
inline
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_degree_3(Vertex_handle v)
{
  remove_degree_3(v, NULL);
}


template< class Gt, class Tds >
inline
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_degree_3(Vertex_handle v, Face* f)
{
  CGAL_triangulation_precondition( v->degree() == 3 );
  _tds.remove_degree_3(v.ptr(), f);
}

//--------------------------------------------------------------------
// insertion of weighted point
//--------------------------------------------------------------------

template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
insert_first(const Weighted_point& p)
{
  CGAL_triangulation_precondition( number_of_vertices() == 0 );
  return Delaunay_graph::insert(p);
}

template< class Gt, class Tds >
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
insert_second(const Weighted_point& p)
{
  CGAL_triangulation_precondition( number_of_vertices() == 1 );
  Vertex_handle vnew;
  Vertex_handle v(vertices_begin());
  if ( awvd_distance2_test(v->point(), p) != POSITIVE ) {
    if ( awvd_compare_weights(v->point(), p) != SMALLER ) {
      v->add_weighted_point(p);
      vnew = Vertex_handle(NULL);  
    } else {
      vnew = Delaunay_graph::insert(p);
      vnew->add_weighted_point(v->point());
      move_trivial_weighted_points(v, vnew);
      remove_first(v);
    }
  } else {
    vnew = Delaunay_graph::insert(p);
  }

  return vnew;
}

template< class Gt, class Tds >
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
insert_third(const Weighted_point& p,
	     std::vector<Vh_triple*>* flipped_edges)
{
  CGAL_triangulation_precondition( number_of_vertices() == 2 );

  Vertex_handle v1(vertices_begin());
  Vertex_handle v2(++vertices_begin());

  Sign s1 = awvd_distance2_test(v1->point(), p);
  Sign s2 = awvd_distance2_test(v2->point(), p);

  if ( s1 != POSITIVE &&
       awvd_compare_weights(v1->point(), p) != SMALLER ) {
    v1->add_weighted_point(p);
    return Vertex_handle(NULL);
  }
  if ( s2 != POSITIVE &&
       awvd_compare_weights(v2->point(), p) != SMALLER ) {
    v2->add_weighted_point(p);
    return Vertex_handle(NULL);
  }
  if ( s1 != POSITIVE && s2 == POSITIVE ) {
    Weighted_point p1 = v1->point();
    Vertex_handle v = Delaunay_graph::insert(p);
    v->add_weighted_point(p1);
    move_trivial_weighted_points(v1, v);
    remove_second(v1);
    return v;
  } else if ( s1 == POSITIVE && s2 != POSITIVE ) {
    Weighted_point p2 = v2->point();
    Vertex_handle v = Delaunay_graph::insert(p);
    v->add_weighted_point(p2);
    move_trivial_weighted_points(v2, v);
    remove_second(v2);
    return v;
  } else if ( s1 != POSITIVE && s2 != POSITIVE ) {
    Weighted_point p1 = v1->point();
    Weighted_point p2 = v2->point();
    Vertex_handle v = Delaunay_graph::insert(p);
    v->add_weighted_point(p1);
    v->add_weighted_point(p2);
    move_trivial_weighted_points(v1, v);
    move_trivial_weighted_points(v2, v);
    remove_first(v1);
    remove_second(v2);
    return v;
  }

  CGAL_triangulation_assertion( is_valid() );

  Vertex_handle v = Delaunay_graph::insert(p);

  Edge_conflict_test_result r =
    awvd_finite_edge_test_degenerated(v1->point(), v2->point(), p);

  if ( r == NO_CONFLICT ) {
    Comparison_result cr =
      awvd_compare_distances_test(v1->point(), v2->point(), p);

    CGAL_assertion( cr != EQUAL );
    Vertex_handle vv = ( cr == LARGER ) ? v1 : v2;

    Face_circulator fc = incident_faces(v);
    while ( true ) {
      Face_handle f(fc);
      int k = f->index(v);
      Vertex_handle vh = f->vertex(ccw(k));
      if ( vh == vv ) {
	flip(f, cw(k));
	break;
      }
      ++fc;
    }
  } else if ( r == INTERIOR ) {
    Edge_circulator ec = incident_edges(v);

    while ( true ) {
      if ( is_infinite(ec) ) {
	flip(*ec);
	break;
      }
      ec++;
    }
  } else if ( r == ENTIRE_EDGE ) {
    Face_circulator fc = incident_faces(v);

    while ( true ) {
      Face_handle f(fc);
      if ( !is_infinite(f) ) {
	flip(f, f->index(v));
	break;
      }
      ++fc;
    }
  } else if ( r == BOTH_VERTICES ) {
    Comparison_result cr =
      awvd_compare_distances_test(v1->point(), v2->point(), p);

    CGAL_assertion( cr != EQUAL );

    Edge_circulator ec;
    ec = ( cr == SMALLER ) ? incident_edges(v1) : incident_edges(v2);
    while ( true ) {
      if ( is_infinite(ec) ) {
	flip(*ec);
	break;
      }
      ec++;
    }
  } else {
    CGAL_assertion( r == RIGHT_VERTEX || r == LEFT_VERTEX );
    // do nothing here
  }

  CGAL_triangulation_assertion( is_valid() );

  return v;
}



template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
insert(const Weighted_point& p)
{
  return insert(p, NULL, NULL);
}



template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
insert(const Weighted_point& p, Vertex_handle vnearest)
{
  return insert(p, vnearest.ptr(), NULL);
}

template< class Gt, class Tds >
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
insert(const Weighted_point& p, Vertex_handle vnearest,
       std::vector<Vh_triple*>* flipped_edges)
{
  if ( number_of_vertices() == 0 ) {
    return insert_first(p);
  }
  if ( number_of_vertices() == 1 ) {
    return insert_second(p);
  }
  if ( number_of_vertices() == 2 ) {
    return insert_third(p, flipped_edges);
  }

  // first find the nearest neighbor
  if ( vnearest.ptr() == NULL ) {
    vnearest = find_nearest_neighbor(p);
  }

  CGAL_assertion( vnearest.ptr() != NULL );


  // check if it is trivial
  Weighted_point wp_nearest = vnearest->point();
  if ( awvd_distance2_test(wp_nearest, p) != POSITIVE &&
       (CGAL_NTS compare(wp_nearest.weight(), p.weight()) != SMALLER) ) {
    vnearest->add_weighted_point(p);
    return Vertex_handle(NULL);
  }

  // find the first conflict
  Edge_circulator ec_start = vnearest->incident_edges();
  Edge_circulator ec = ec_start;

  Edge_conflict_test_result r;
  Edge e;
  do {
    e = *ec;
    r = edge_conflict_test(e, p);

    if ( r != NO_CONFLICT ) { break; }
    ++ec;
  } while ( ec != ec_start );

  CGAL_assertion( r != NO_CONFLICT );

  // no conflict with Voronoi vertices, but conflict with the interior 
  // of a Voronoi edge
  if ( r == INTERIOR ) {
    return insert_degree_2(e, p);
  }

  List l(*this);
  Face_map fm;
  Vertex_map vm;

  Vertex_handle v = static_cast<Vertex*>(_tds.create_vertex());
  v->set_point(p);

  find_conflict_region(p, vnearest, l, fm, vm, flipped_edges);
  retriangulate_conflict_region(v, l, fm, vm);

  fm.clear();
  vm.clear();

  return v;
}


//--------------------------------------------------------------------
// find conflict region
//--------------------------------------------------------------------

template< class Gt, class Tds >
std::vector<Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Edge>
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
find_conflict_region(const Weighted_point& p, Vertex_handle vnearest)
{
  CGAL_precondition ( number_of_vertices() > 2 );

  List l(*this);
  Face_map fm;
  Vertex_map vm;

  find_conflict_region(p, vnearest, l, fm, vm, NULL);

  std::vector<Edge> edge_list;

  while ( !l.is_empty() ) {
    Edge e = l.front();
    l.pop();
    edge_list.push_back(e);
  }

  fm.clear();
  vm.clear();

  return edge_list;
}


template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
find_conflict_region(const Weighted_point& p,
		     const Vertex_handle& vnearest,
		     List& l, Face_map& fm,
		     Vertex_map& vm,
		     std::vector<Vh_triple*>* fe)
{
  // check if it is trivial
  Weighted_point wp_nearest = vnearest->point();
  if ( awvd_distance2_test(wp_nearest, p) != POSITIVE &&
       (CGAL_NTS compare(wp_nearest.weight(), p.weight()) != SMALLER) ) {
    vnearest->add_weighted_point(p);
    return;
  }

  CGAL_precondition( vnearest.ptr() != NULL );

  Edge_circulator ec_start = vnearest->incident_edges();
  Edge_circulator ec = ec_start;

  Edge_conflict_test_result r;
  Edge e;
  do {
    e = *ec;
    r = edge_conflict_test(e, p);

    if ( r != NO_CONFLICT ) { break; }
    ++ec;
  } while ( ec != ec_start );

  CGAL_assertion( r != NO_CONFLICT );

  if ( r == INTERIOR ) {
    l.push(e);
    l.push(sym_edge(e));
    return;
  }


  Face_handle start_f;
  if ( r != RIGHT_VERTEX ) {
    start_f = e.first;
  } else {
    start_f = e.first->neighbor(e.second);
  }

  initialize_conflict_region(start_f, l);
  expand_conflict_region(start_f, p, l, fm, vm, fe);
}

template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
initialize_conflict_region(const Face_handle& f, List& l)
{
  l.clear();
  for (int i = 0; i < 3; i++) {
    l.push(sym_edge(f, i));
  }
}

template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
check_edge_for_trivial_weighted_points(const Face_handle& f, int i,
				       const Weighted_point& p,
				       Vertex_map& vm)
{
  Vertex_handle v1 = f->vertex(ccw(i));
  if ( vm.find(v1) == vm.end() && !is_infinite(v1) ) {
    if ( awvd_distance2_test(v1->point(), p) != POSITIVE ) {
      Comparison_result r;
      r = CGAL_NTS compare(v1->point().weight(), p.weight());
      CGAL_assertion ( r == SMALLER );
      vm[v1] = true;
    }
  }

  Vertex_handle v2 = f->vertex(cw(i));
  if ( vm.find(v2) == vm.end() && !is_infinite(v2) ) {
    if ( awvd_distance2_test(v2->point(), p) != POSITIVE ) {
      Comparison_result r;
      r = CGAL_NTS compare(v2->point().weight(), p.weight());
      CGAL_assertion ( r == SMALLER );
      vm[v2] = true;
    }
  }
}

template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
expand_conflict_region(const Face_handle& f, const Weighted_point& p,
		       List& l, Face_map& fm, Vertex_map& vm,
		       std::vector<Vh_triple*>* fe)
{
  // setting fm[f] to true means that the face has been reached and
  // that the face is available for recycling. If we do not want the
  // face to be available for recycling we must set this flag to
  // false.
  fm[f] = true;

  //  CGAL_assertion( fm.find(f) != fm.end() );
  for (int i = 0; i < 3; i++) {
    check_edge_for_trivial_weighted_points(f, i, p, vm);

    Face_handle n = f->neighbor(i);

    Edge_conflict_test_result r = edge_conflict_test(f, i, p);

    if ( fm.find(n) != fm.end() ) {
      if ( r == ENTIRE_EDGE ) {
	Edge e = sym_edge(f, i);
	if ( l.is_in_list(sym_edge(e)) ) {
	  l.remove(e);
	  l.remove(sym_edge(e));
	}
      }
      continue;
    }

    if ( r == NO_CONFLICT || r == INTERIOR ||
	 r == RIGHT_VERTEX ) {
      std::cout << "r: " << int(r) << std::endl;
      Edge e(f,i);
      std::cout << "is infinite: " << is_infinite(e)
		<< std::endl;
      if ( is_infinite(e) ) {
	Vertex_handle v = find_nearest_neighbor(p);
	Weighted_point nn = v->point();
	Weighted_point p2, p3, p4;
	if ( is_infinite(f->vertex(ccw(i))) ) {
	  p2 = f->vertex(cw(i))->point();
	  p3 = f->vertex(i)->point();
	  p4 = f->mirror_vertex(i)->point();
	} else {
	  p2 = f->vertex(ccw(i))->point();
	  p3 = f->mirror_vertex(i)->point();
	  p4 = f->vertex(i)->point();
	}
	std::cout << "p2: " << p2.point() << "," << p2.weight() << std::endl;
	std::cout << "p3: " << p3.point() << "," << p3.weight() << std::endl;
	std::cout << "p4: " << p4.point() << "," << p4.weight() << std::endl;
	std::cout << " q: " << p.point() << "," << p.weight() << std::endl;
	std::cout << "nn: " << nn.point() << "," << nn.weight() << std::endl;
	std::cout << "trivial list size: " << vm.size() << std::endl;

	Sign s1 = awvd_distance2_test(p3, p);
	Sign s2 = awvd_distance2_test(p4, p);
	Sign s3 = awvd_distance2_test(p4, p2);
	Sign s4 = awvd_distance2_test(p3, p2);
	std::cout << "s1, s2: " << int(s1) << " " << int(s2) << std::endl;
	std::cout << "s3, s4: " << int(s3) << " " << int(s4) << std::endl;
      }

    }
    CGAL_assertion( r != NO_CONFLICT && r != INTERIOR &&
		    r != RIGHT_VERTEX );
    if ( r == ENTIRE_EDGE ) {
      Edge e = sym_edge(f, i);

      CGAL_assertion( l.is_in_list(e) );
      int j = f->mirror_index(i);
      Edge e_before = sym_edge(n, ccw(j));
      Edge e_after = sym_edge(n, cw(j));
      if ( !l.is_in_list(e_before) ) {
	l.insert_before(e, e_before);
      }
      if ( !l.is_in_list(e_after) ) {
	l.insert_after(e, e_after);
      }
      l.remove(e);

      if ( fe != NULL ) {
	Vh_triple* vhq = new Vh_triple[1];

	(*vhq)[0] = NULL;
	(*vhq)[1] = n->vertex(     j  );
	(*vhq)[2] = n->vertex( ccw(j) );

	fe->push_back(vhq);
      }
      expand_conflict_region(n, p, l, fm, vm, fe);
    }

  }
}

//--------------------------------------------------------------------
// retriangulate conflict region
//--------------------------------------------------------------------

template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
add_bogus_vertex(Edge e, List& l)
{
  Edge esym = sym_edge(e);
  Face_handle g1 = e.first;
  Face_handle g2 = esym.first;

  Vertex_handle v = insert_degree_2(e);
  Face_circulator fc(v);
  Face_handle f1(fc);
  Face_handle f2(++fc);
  int i1 = f1->index(v);
  int i2 = f2->index(v);

  CGAL_assertion( ((f1->neighbor(i1) == g1) && (f2->neighbor(i2) == g2)) ||
		  ((f1->neighbor(i1) == g2) && (f2->neighbor(i2) == g1)) );

  Edge ee, eesym;
  if ( f1->neighbor(i1) == g1 ) {
    ee = Edge(f2, i2);
    eesym = Edge(f1, i1);
  } else {
    ee = Edge(f1, i1);
    eesym = Edge(f2, i2);
  }

  l.replace(e, ee);
  l.replace(esym, eesym);

  return v;
}

template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_list
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
add_bogus_vertices(List& l)
{
  Vertex_list vertex_list;
  Edge_list edge_list;

  Edge e_start = l.front();
  Edge e = e_start;
  do {
    Edge esym = sym_edge(e);
    if ( l.is_in_list(esym) ) {
      bool found(false);
      Edge_list::iterator it;
      for (it = edge_list.begin(); it != edge_list.end(); ++it) {
	if ( *it == esym ) {
	  found = true;
	  break;
	}
      }
      if ( !found ) { edge_list.push_back(e); }
    }
    e = l.next(e);
  } while ( e != e_start );

  Edge_list::iterator it;

  for (it = edge_list.begin();  it != edge_list.end(); ++it) {
    e = *it;
    Vertex_handle v = add_bogus_vertex(e, l);
    vertex_list.push_back(v);
  }

  return vertex_list;
}

template< class Gt, class Tds >
inline
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_bogus_vertices(Vertex_list& vl)
{
  while ( vl.size() > 0 ) {
    Vertex_handle v = vl.front();
    vl.pop_front();
    remove_degree_2(v);
  }
}


template< class Gt, class Tds >
inline
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
move_trivial_weighted_points(Vertex_handle& vold, Vertex_handle& vnew)
{
  typename Vertex_base::Weighted_point_list_iterator wpit;

  for (wpit = vold->weighted_points_begin();
       wpit != vold->weighted_points_end(); ++wpit) {
    vnew->add_weighted_point(*wpit);
  }

  vold->clear_weighted_point_list();
}




template< class Gt, class Tds >
std::vector<Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Face*>
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
get_faces_for_recycling(Face_map& fm, unsigned int n_wanted)
{
  std::vector<Face*> vf;

  Face_map::iterator fmit;
  for (fmit = fm.begin(); fmit != fm.end(); ++fmit) {
    Face_handle f = (*fmit).first;
    if ( fm[f] == true ) { vf.push_back(f.ptr()); }
  }

  while ( vf.size() < n_wanted ) {
    Face* fp = static_cast<Face*>(_tds.create_face());
    vf.push_back(fp);
  }

  while ( vf.size() > n_wanted ) {
    Face* fp = vf.back();
    vf.pop_back();
    _tds.delete_face(fp);
  }
  
  return vf;
}

template< class Gt, class Tds>
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_trivial_vertices(Vertex_handle&v, Vertex_map& vm, Face_map& fm)
{
  Vertex_map::iterator it;

  for (it = vm.begin(); it != vm.end(); ++it) {
    Vertex_handle vtrivial = (*it).first;
    _tds.delete_vertex(vtrivial.ptr());
  }
  vm.clear();
}


template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
retriangulate_conflict_region(Vertex_handle& v,	List& l, 
			      Face_map& fm, Vertex_map& vm)
{
  //  CGAL_assertion( is_valid() );

  // 1. move all the trivial weighted points to the new one
  Vertex_map::iterator vmit;
  for (vmit = vm.begin(); vmit != vm.end(); ++vmit) {
    Vertex_handle vtrivial = (*vmit).first;
    move_trivial_weighted_points(vtrivial, v);
    v->add_weighted_point(vtrivial->point());
  }

  int vmsize = vm.size();
  int num_vert = number_of_vertices();
  CGAL_precondition( number_of_vertices() - vm.size() >= 2 );

  // 2. add the bogus vetrices
  Vertex_list dummy_vertices = add_bogus_vertices(l);

  // 3. repair the face pointers...
  Edge e_start = l.front();
  Edge eit = e_start;
  do {
    Edge esym = sym_edge(eit);
    Face_handle f = eit.first;
    int k = eit.second;
    CGAL_assertion( !l.is_in_list(esym) );
    CGAL_assertion( fm.find(f) == fm.end() );
    f->vertex(ccw(k))->set_face(f);
    f->vertex( cw(k))->set_face(f);
    eit = l.next(eit);
  } while ( eit != e_start );

  //  std::vector<Face*> vf = get_faces_for_recycling(fm, l.size());
  std::list<Face*> vf;

  // 4. copy the edge list to a vector of edges and clear the in place 
  //    list
  typedef typename Tds::Edge Tds_edge;
  std::vector<Tds_edge> ve;

  Edge efront = l.front();
  Edge e = efront;
  do {
    ve.push_back(Tds_edge(e.first.ptr(), e.second));
    e = l.next(e);
  } while ( e != efront );

  l.clear();

  // 5. remove the trivial vertices
  remove_trivial_vertices(v, vm, fm);

  // 6. retriangulate the hole
  //  _tds.star_hole(v.ptr(), ve.begin(), ve.end(), vf.begin(), vf.end());
  _tds.star_hole(v.ptr(), ve.begin(), ve.end());

  // 7. remove the bogus vertices
  remove_bogus_vertices(dummy_vertices);

  // 8. remove the unused faces
  Face_map::iterator it;
  for (it = fm.begin(); it != fm.end(); ++it) {
    Face_handle fh = (*it).first;
    _tds.delete_face(fh.ptr());
  }

  CGAL_assertion( number_of_vertices() == num_vert - vmsize );

  // 9. DONE!!!!
}



//--------------------------------------------------------------------
// point location
//--------------------------------------------------------------------
template< class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
find_nearest_neighbor(const Point& p) const
{
#if 0
  if ( number_of_vertices() == 0 ) {
    return Vertex_handle(NULL);
  }

  Vertex_iterator vit = vertices_begin();
  Vertex_handle start_vertex;
  for (; vit != vertices_end(); ++vit) {
    if ( !is_infinite(vit) ) {
      start_vertex = Vertex_handle(vit);
      break;
    }
  }
  return find_nearest_neighbor(p, start_vertex);
#endif
  return find_nearest_neighbor(p, Vertex_handle(NULL));
}


template< class Gt, class Tds >
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Vertex_handle
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
find_nearest_neighbor(const Point& p,
		      Vertex_handle start_vertex) const
{
  if ( number_of_vertices() == 0 ) {
    return Vertex_handle(NULL);
  }

  if ( start_vertex.ptr() == NULL ) {
    Vertex_iterator vit = vertices_begin();
    for (; vit != vertices_end(); ++vit) {
      if ( !is_infinite(vit) ) {
	start_vertex = Vertex_handle(vit);
	break;
      }
    }
  }

  //  if ( start_vertex.ptr() == NULL ) { return start_vertex; }

  Vertex_handle vclosest;
  Vertex_handle v = start_vertex;

  if ( number_of_vertices() < 3 ) {
    vclosest = v;
    Vertex_iterator vit = vertices_begin();
    for (; vit != vertices_end(); ++vit) {
      Vertex_handle v1(vit);
      if ( v1 != vclosest && !is_infinite(v1) ) {
	Weighted_point p1 = vclosest->point();
	Weighted_point p2 = v1->point();
	if ( awvd_compare_distances_test(p1, p2, p) == LARGER ) {
	  vclosest = v1;
	}
      }
    }
    return vclosest;
  }

#if 1
  do {
    vclosest = v;
    Weighted_point p1 = v->point();
    Vertex_circulator vc_start = incident_vertices(v);
    Vertex_circulator vc = vc_start;
    do {
      if ( !is_infinite(vc) ) {
	Vertex_handle v1(vc);
	Weighted_point p2 = v1->point();
	if ( awvd_compare_distances_test(p1, p2, p) == LARGER ) {
	  v = v1;
	  break;
	}
      }
      ++vc;
    } while ( vc != vc_start );
  } while ( vclosest != v );
#else
  vclosest = v;
  Weighted_point p1 = v->point();
  Face_circulator fc(v);
  Face_circulator fc_start = fc;
  Face_handle ff;
  do {
    ff = Face_handle(fc);
    Vertex_handle v1 = ff->vertex( ccw(ff->index(v)) );
    if ( !is_infinite(v1) ) {
      Weighted_point p2 = v1->point();
      if ( awvd_compare_distances_test(p1, p2, p) == LARGER ) {
	v = v1;
	break;
      }
    }
    ++fc;
  } while ( fc != fc_start );
  if ( vclosest == v ) { return vclosest; }

  do {
    vclosest = v;
    Weighted_point p1 = v->point();
    fc = Face_circulator(v, ff);
    fc_start = fc;
    do {
      ff = Face_handle(fc);
      Vertex_handle v1 = ff->vertex( ccw(ff->index(v)) );
      if ( !is_infinite(v1) ) {
	Weighted_point p2 = v1->point();
	if ( awvd_compare_distances_test(p1, p2, p) == LARGER ) {
	  v = v1;
	  break;
	}
      }
      ++fc;
    } while ( fc != fc_start );
  } while ( vclosest != v );
#endif

  return vclosest;
}


//----------------------------------------------------------------------
// methods for the predicates
//----------------------------------------------------------------------
template < class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_distance2_test(const Weighted_point &p, const Weighted_point &q) const
{
  return geom_traits().awvd_distance2_test_2_object()(p, q);
}


template < class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_test(const Weighted_point &p, const Weighted_point &q,
	  const Weighted_point &r, const Weighted_point &s) const
{
  return geom_traits().awvd_test_2_object()(p, q, r, s);
}

template < class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_test_degenerated(const Weighted_point &p, const Weighted_point &q,
		      const Weighted_point &r) const
{
  return geom_traits().awvd_test_degenerated_2_object()(p, q, r);
}

template < class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_Voronoi_radii_difference_test(const Weighted_point &p,
				   const Weighted_point &q,
				   const Weighted_point &r, 
				   const Weighted_point &s) const
{
  return
    geom_traits().awvd_sign_of_difference_of_Voronoi_radii_test_2_object()
    (p, q, r, s);
}

template < class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_test(const Face_handle& f,	const Weighted_point &p) const
{
  CGAL_precondition( !is_infinite(f) );
  return awvd_test(f->vertex(0)->point(),
		   f->vertex(1)->point(),
		   f->vertex(2)->point(), p);
}

template < class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_test(const Face_handle& f,	int i, const Weighted_point &p) const
{
  CGAL_triangulation_precondition( is_infinite(f) );
  CGAL_triangulation_precondition( is_infinite(f->vertex(i)) );
  CGAL_triangulation_precondition
    ( !is_infinite(  f->vertex(ccw(i)) ) &&
      !is_infinite(  f->vertex( cw(i)) ) );

  return awvd_test_degenerated( f->vertex( ccw(i) )->point(),
				f->vertex(  cw(i) )->point(),
				p );
}

template < class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_Voronoi_radii_difference_test(const Face_handle& f, int i,
				   const Weighted_point &p) const
{
  CGAL_triangulation_precondition( !is_infinite(f) );
  return awvd_Voronoi_radii_difference_test(f->vertex(ccw(i))->point(),
					    f->vertex( cw(i))->point(),
					    f->vertex(i)->point(), p);
}


template < class Gt, class Tds >
inline
Orientation
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_orientation_test(const Weighted_point &p1,
		      const Weighted_point &p2,
		      const Weighted_point &p3,
		      const Weighted_point &q1,
		      const Weighted_point &q2) const
{
  return geom_traits().awvd_orientation_test_2_object()
    (p1, p2, p3, q1.point(), q2.point());
}

template < class Gt, class Tds >
inline
Orientation
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_orientation_test(const Face_handle& f,
		      const Weighted_point& p,
		      const Weighted_point& q) const
{
  CGAL_triangulation_precondition( !is_infinite(f) );
  return awvd_orientation_test(f->vertex(0)->point(),
			       f->vertex(1)->point(),
			       f->vertex(2)->point(),
			       p, q);
}

template < class Gt, class Tds >
inline
Orientation
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_orientation_left_test(const Weighted_point &p1,
			   const Weighted_point &p2,
			   const Point &q1, const Point &q2) const
{
  return
    geom_traits().awvd_orientation_using_left_tangent_point_test_2_object()
    (p1, p2, q1, q2);
}

template < class Gt, class Tds >
inline
Orientation
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_orientation_right_test(const Weighted_point &p1,
			    const Weighted_point &p2,
			    const Point &q1, const Point &q2) const
{
  return
    geom_traits().awvd_orientation_using_right_tangent_point_test_2_object()
    (p1, p2, q1, q2);
}



template < class Gt, class Tds >
inline
Bounded_side
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_circle_existence_test(const Weighted_point &p,
			   const Weighted_point &q,
			   const Weighted_point &r) const
{
  return geom_traits().awvd_circle_existence_test_2_object()(p, q, r);
}

template < class Gt, class Tds >
inline
Bounded_side
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_circle_existence_test(const Face_handle& f) const
{
  CGAL_triangulation_precondition( !is_infinite(f) );
  return awvd_circle_existence_test(f->vertex(0)->point(),
				    f->vertex(1)->point(),
				    f->vertex(2)->point());
}

template < class Gt, class Tds >
inline
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Order
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_ordering_test(const Weighted_point &p,
		   const Weighted_point &q,
		   const Weighted_point &r) const
{
  int o = geom_traits().awvd_ordering_on_left_bitangent_test_2_object()
    (p, q, r);
  if ( o == 1 ) return FIRST;
  if ( o == 2 ) return SECOND;
  CGAL_assertion( o == 3 );
  return THIRD;
}

template < class Gt, class Tds >
inline
Comparison_result
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_compare_distances_test(const Weighted_point &p1,
			    const Weighted_point &p2,
			    const Point &p) const
{
  return geom_traits().awvd_compare_distances_2_object()(p1, p2, p);
}


//*********************************************

template< class Gt, class Tds >
inline
Comparison_result
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_compare_weights(const Weighted_point& p,
		     const Weighted_point& q) const
{
  return geom_traits().awvd_compare_weight_test_2_object()(p, q);
}

template< class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_vertex_test(const Weighted_point &p1, const Weighted_point &p2,
		 const Weighted_point &p3, const Weighted_point &q) const

{
  return geom_traits().awvd_test_2_object()(p1, p2, p3, q);
}

template < class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_vertex_test_degenerated(const Weighted_point &p1,
			     const Weighted_point &p2,
			     const Weighted_point &q) const
{
  return geom_traits().awvd_test_degenerated_2_object()(p1, p2, q);
}

template < class Gt, class Tds >
inline
Sign
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_vertex_test(const Face_handle& f, const Weighted_point &p) const
{
  if ( !is_infinite(f) ) {
    return awvd_vertex_test(f->vertex(0)->point(),
			    f->vertex(1)->point(),
			    f->vertex(2)->point(), p);
  }
  int inf_i;
  for (int i = 0; i < 3; i++) {
    if ( is_infinite(f->vertex(i)) ) {
      inf_i = i;
      break;
    }
  }
  return awvd_vertex_test_degenerated( f->vertex( ccw(inf_i) )->point(),
				       f->vertex(  cw(inf_i) )->point(),
				       p );
}

template< class Gt, class Tds >
inline
Edge_conflict_test_result
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_finite_edge_test(const Weighted_point& p1,
		      const Weighted_point& p2,
		      const Weighted_point& p3,
		      const Weighted_point& p4,
		      const Weighted_point& q) const
{
  return geom_traits().aw_Voronoi_diagram_finite_edge_test_2_object()
    (p1, p2, p3, p4, q);
}

template< class Gt, class Tds >
inline
Edge_conflict_test_result 
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_finite_edge_test(const Face_handle& f, int i,
		      const Weighted_point& p) const
{
  CGAL_precondition( !is_infinite(f) &&
		     !is_infinite(f->neighbor(i)) );
  return awvd_finite_edge_test( f->vertex( ccw(i) )->point(),
			        f->vertex(  cw(i) )->point(),
				f->vertex(     i  )->point(),
				f->mirror_vertex(i)->point(), p);
}

template< class Gt, class Tds >
inline
Edge_conflict_test_result
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_finite_edge_test_degenerated(const Weighted_point& p1,
				  const Weighted_point& p2,
				  const Weighted_point& p3,
				  const Weighted_point& q) const
{
  return
    geom_traits().aw_Voronoi_diagram_finite_edge_test_degenerated_2_object()
    (p1, p2, p3, q);
}

template< class Gt, class Tds >
inline
Edge_conflict_test_result
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_finite_edge_test_degenerated(const Weighted_point& p1,
				  const Weighted_point& p2,
				  const Weighted_point& q) const
{
  return
    geom_traits().aw_Voronoi_diagram_finite_edge_test_degenerated_2_object()
    (p1, p2, q);
}


template< class Gt, class Tds >
Edge_conflict_test_result 
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_finite_edge_test_degenerated(const Face_handle& f, int i,
				  const Weighted_point& p) const
{
  if ( !is_infinite( f->mirror_vertex(i) ) ) {
    CGAL_precondition( is_infinite(f->vertex(i)) );

    Face_handle g = f->neighbor(i);
    int j = f->mirror_index(i);

    return opposite(awvd_finite_edge_test_degenerated(g, j, p));
  }

  CGAL_precondition( is_infinite( f->mirror_vertex(i) ) );

  Weighted_point p1 = f->vertex( ccw(i) )->point();
  Weighted_point p2 = f->vertex(  cw(i) )->point();

  if ( is_infinite(f->vertex(i)) ) {
    return awvd_finite_edge_test_degenerated(p1, p2, p);
  }

  Weighted_point p3 = f->vertex(i)->point();
  return awvd_finite_edge_test_degenerated(p1, p2, p3, p);
}

template< class Gt, class Tds >
Edge_conflict_test_result 
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_infinite_edge_test(const Weighted_point& p2,
			const Weighted_point& p3,
			const Weighted_point& p4,
			const Weighted_point& q) const
{
  return
    geom_traits().aw_Voronoi_diagram_infinite_edge_test_2_object()
    (p2, p3, p4, q);
}

template< class Gt, class Tds >
Edge_conflict_test_result 
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
awvd_infinite_edge_test(const Face_handle& f, int i,
			const Weighted_point& p) const
{
  if ( !is_infinite( f->vertex(ccw(i)) ) ) {
    CGAL_precondition( is_infinite( f->vertex(cw(i)) ) );
    Face_handle g = f->neighbor(i);
    int j = f->mirror_index(i);

    return opposite(awvd_infinite_edge_test(g, j, p));
  }

  CGAL_precondition( is_infinite( f->vertex(ccw(i)) ) );

  Weighted_point p2 = f->vertex(  cw(i) )->point();
  Weighted_point p3 = f->vertex(     i  )->point();
  Weighted_point p4 = f->mirror_vertex(i)->point();

  return awvd_infinite_edge_test(p2, p3, p4, p);
}

//*********************************************

template< class Gt, class Tds >
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::Conflict_test_result
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
conflict_test_for_Voronoi_center(const Face_handle& f,
				 const Weighted_point& p) const
{
  bool is_inf = is_infinite(f);
  if ( !is_inf && awvd_circle_existence_test(f) != ON_UNBOUNDED_SIDE ) {
    return FLIP;
  }

  if ( !is_inf ) {
    Sign s = awvd_test(f, p);
    return (s == NEGATIVE) ? FLIP : DO_NOT_FLIP;
  }

  int idx;
  for (int i = 0; i < 3; i++) {
    if (  is_infinite( f->vertex(i) )  ) {
      idx = i;
      break;
    }
  }

  Sign s = awvd_test(f, idx, p);
  return (s == NEGATIVE) ? FLIP : DO_NOT_FLIP;
}



template< class Gt, class Tds >
Edge_conflict_test_result
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
edge_conflict_test(const Face_handle& f, int i,
		   const Weighted_point& p) const
{
  Vertex_handle v1 = f->vertex(ccw(i));
  if ( !is_infinite(v1) ) {
    if ( awvd_distance2_test(v1->point(), p) != POSITIVE ) {
      return ENTIRE_EDGE;
    }
  }

  Vertex_handle v2 = f->vertex(cw(i));
  if ( !is_infinite(v2) ) {
    if ( awvd_distance2_test(v2->point(), p) != POSITIVE ) {
      return ENTIRE_EDGE;
    }
  }

  Face_handle g = f->neighbor(i);

  bool is_inf_f = is_infinite(f);
  bool is_inf_g = is_infinite(g);

  Edge_conflict_test_result result;

  if ( !is_inf_f && !is_inf_g ) {
    result = awvd_finite_edge_test(f, i, p);
  } else if ( !is_inf_f || !is_inf_g ) {
    result = awvd_finite_edge_test_degenerated(f, i, p);
  } else {
    //    Edge e(f, i);
    if ( !is_infinite(f, i) ) {
      result = awvd_finite_edge_test_degenerated(f, i, p);
    } else {
      result = awvd_infinite_edge_test(f, i, p);
    }
  }

  return result;
}

template< class Gt, class Tds >
Edge_conflict_test_result
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
edge_conflict_test(const Face_handle& f, int i,
		   const Vertex_handle& v4,
		   const Weighted_point& p) const
{
  CGAL_triangulation_precondition( v4 != f->vertex(ccw(i)) );
  CGAL_triangulation_precondition( v4 != f->vertex( cw(i)) );

  Vertex_handle v1 = f->vertex(ccw(i));
  if ( !is_infinite(v1) ) {
    if ( awvd_distance2_test(v1->point(), p) != POSITIVE ) {
      return ENTIRE_EDGE;
    }
  }

  Vertex_handle v2 = f->vertex(cw(i));
  if ( !is_infinite(v2) ) {
    if ( awvd_distance2_test(v2->point(), p) != POSITIVE ) {
      return ENTIRE_EDGE;
    }
  }

  if ( !is_infinite(v1) && !is_infinite(v2) && !is_infinite(v4) ) {
    Weighted_point p1 = v1->point();
    Weighted_point p2 = v2->point();
    Weighted_point p4 = v4->point();
    if ( awvd_circle_existence_test(p1, p4, p2) != ON_UNBOUNDED_SIDE ) 
      {
	return ENTIRE_EDGE;
      }
  }

  bool is_inf_f = is_infinite(f);
  bool is_inf_v4 = is_infinite(v4);

  Edge_conflict_test_result result;

  if ( !is_inf_f && !is_inf_v4 ) {
    result = awvd_finite_edge_test(f->vertex( ccw(i) )->point(),
				   f->vertex(  cw(i) )->point(),
				   f->vertex(     i  )->point(),
				   v4->point(), p);
  } else if ( !is_inf_f && is_inf_v4 ) {
    result =
      awvd_finite_edge_test_degenerated(f->vertex( ccw(i) )->point(),
					f->vertex(  cw(i) )->point(),
					f->vertex(     i  )->point(), p);
  } else if ( is_inf_f && is_inf_v4 ) {
    CGAL_triangulation_precondition( !is_infinite(f->vertex(ccw(i))) );
    CGAL_triangulation_precondition( !is_infinite(f->vertex( cw(i))) );
    result =
      awvd_finite_edge_test_degenerated(f->vertex( ccw(i) )->point(),
					f->vertex(  cw(i) )->point(), p);
  } else {
    CGAL_assertion( is_inf_f && !is_inf_v4 );
    if ( is_infinite(f->vertex(i)) ) {
      result =
	awvd_finite_edge_test_degenerated(f->vertex(  cw(i) )->point(),
					  f->vertex( ccw(i) )->point(),
					  v4->point(), p);
    } else if ( is_infinite(f->vertex(ccw(i))) ) {
      result = awvd_infinite_edge_test(f->vertex(  cw(i) )->point(),
				       f->vertex(     i  )->point(),
				       v4->point(), p);
    } else {
      result = awvd_infinite_edge_test(f->vertex(  ccw(i) )->point(),
				       v4->point(),
				       f->vertex(      i  )->point(), p);
    }
  }

  return result;
}


//----------------------------------------------------------------------
// methods for disk removal
//----------------------------------------------------------------------


template< class Gt, class Tds >
inline
bool
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
is_degree_2_vertex(const Vertex_handle& v) const
{
  Face_circulator fc = v->incident_faces();
  Face_circulator fc1 = v->incident_faces();
  ++(++fc1);
  Face_handle fh(fc), fh1(fc1);
  return ( fh == fh1 );
}


template< class Gt, class Tds >
inline
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_first(Vertex_handle v)
{
  Delaunay_graph::remove(v);
}

template< class Gt, class Tds >
inline
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_second(Vertex_handle v)
{
  Delaunay_graph::remove(v);
}

template< class Gt, class Tds >
inline
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_third(Vertex_handle v)
{
  if ( is_degree_2_vertex(v) ) {
    Face_handle fh(v->incident_faces());
    int i = fh->index(v);
    flip(fh, i);
  } else if ( v->degree() == 4 ) {
    Edge_circulator ec = v->incident_edges();
    for (int i = 0; i < 4; i++) {
      Edge e = *ec;
      Edge sym = sym_edge(e);
      if ( e.first->vertex(e.second) !=	sym.first->vertex(sym.second) ) {
	flip(e);
	break;
      }
      ++ec;
    }
  }

  _tds.remove_dim_down(v.ptr());
}




template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));

  Weighted_point_list wp_list;
  typename Vertex_base::Weighted_point_list_iterator wpit;

  for (wpit = v->weighted_points_begin();
       wpit != v->weighted_points_end(); ++wpit) {
    wp_list.push_back(*wpit);
  }

  if ( wp_list.size() > 1 ) {
    Weighted_point_less_than_comparator less_than(geom_traits());
    std::sort(wp_list.begin(), wp_list.end(), less_than);
  }


  int n = number_of_vertices();
  if ( n == 1 ) {
    remove_first(v);
  } else if ( n == 2 ) {
    remove_second(v);
  } else if ( n == 3 ) {
    remove_third(v);
  } else {
    int degree = v->degree();
    if ( degree == 2 ) {
      remove_degree_2(v);
    } else if ( degree == 3 ) {
      remove_degree_3(v);
    } else {
      remove_degree_d_vertex(v);
    }
  }

  for (unsigned int i = 0; i < wp_list.size(); i++) {
    insert(wp_list[i]);
  }

  //  CGAL_triangulation_assertion( is_valid(false, 2) );
  return;
}


template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_using_conflict_region(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));

  Weighted_point_list wp_list;
  typename Vertex_base::Weighted_point_list_iterator wpit;

  for (wpit = v->weighted_points_begin();
       wpit != v->weighted_points_end(); ++wpit) {
    wp_list.push_back(*wpit);
  }

  if ( wp_list.size() > 1 ) {
    Weighted_point_less_than_comparator less_than(geom_traits());
    std::sort(wp_list.begin(), wp_list.end(), less_than);
  }


  int n = number_of_vertices();
  if ( n == 1 ) {
    remove_first(v);
  } else if ( n == 2 ) {
    remove_second(v);
  } else if ( n == 3 ) {
    remove_third(v);
  } else {
    int degree = v->degree();
    if ( degree == 2 ) {
      remove_degree_2(v);
    } else if ( degree == 3 ) {
      remove_degree_3(v);
    } else {
      remove_degree_d_vertex_using_conflict_region(v);
    }
  }

  for (unsigned int i = 0; i < wp_list.size(); i++) {
    insert(wp_list[i]);
  }

  //  CGAL_triangulation_assertion( is_valid(false, 2) );
  return;
}


template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_degree_d_vertex(Vertex_handle v)
{
  Additively_weighted_Voronoi_diagram_2< Gt, Tds > awvd_small;

  std::map<Vertex_handle,Vertex_handle> vmap;

  Vertex_circulator vc_start = v->incident_vertices();
  Vertex_circulator vc = v->incident_vertices();
  Vertex_handle vh_large, vh_small;
  do {
    vh_large = Vertex_handle(vc);
    if ( is_infinite(vh_large) ) {
      vh_small = awvd_small.infinite_vertex();
      vmap[vh_small] = vh_large;
    } else { 
      vh_small = awvd_small.insert(vc->point());
      if ( vh_small.ptr() != NULL ) {
	vmap[vh_small] = vh_large;
      }
    }
    ++vc;
  } while ( vc != vc_start );

  std::vector<Vh_triple*> flipped_edges;

  vh_small = awvd_small.insert(v->point(), NULL, &flipped_edges);
  vmap[vh_small] = v;

  std::ofstream f("data/data_small.out");
  assert(f);
  awvd_small.write_non_trivial_weighted_points(f);
  f.close();

  //  awvd_small.is_valid(true, 2);

  if ( v->degree() != vh_small->degree() ) {
#if DO_PRINTOUTS
    std::cout << "degrees (small, large): "
	      << vh_small->degree() << " " << v->degree()
	      << std::endl;
#endif
  }
  //  CGAL_triangulation_assertion( v->degree() == vh_small->degree() );


  Edge_circulator ec;

  unsigned int num_fe = flipped_edges.size();
  for (unsigned int i = 0; i < num_fe; i++) {
    Vh_triple *vhq = flipped_edges[num_fe - i - 1];

    bool found(false);
    ec = v->incident_edges();
    Edge_circulator ec_start = ec;
    do {
      Edge e = *ec;
      if ( (e.first->vertex( ccw(e.second) ) == vmap[(*vhq)[0]] &&
	    e.first->vertex(  cw(e.second) ) == vmap[(*vhq)[1]] &&
	    e.first->vertex(     e.second  ) == vmap[(*vhq)[2]]) ||
	   (e.first->vertex(  cw(e.second) ) == vmap[(*vhq)[0]] &&
	    e.first->vertex( ccw(e.second) ) == vmap[(*vhq)[1]] &&
	    e.first->mirror_vertex(e.second) == vmap[(*vhq)[2]]) ) {
	flip(e);
	found = true;
	break;
      }
      ++ec;
    } while ( ec != ec_start );
    CGAL_assertion( found );
  }
  CGAL_triangulation_precondition( v->degree() == 3 );

  _tds.remove_degree_3(v.ptr(), NULL);

  for (unsigned int i = 0; i < num_fe; i++) {
    delete flipped_edges[i];
  }
}

template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_degree_d_vertex_using_conflict_region(Vertex_handle v)
{
  Additively_weighted_Voronoi_diagram_2< Gt, Tds > awvd_small;

  std::map<Vertex_handle,Vertex_handle> vmap;

  Vertex_circulator vc_start = v->incident_vertices();
  Vertex_circulator vc = v->incident_vertices();
  Vertex_handle vh_large, vh_small;
  do {
    vh_large = Vertex_handle(vc);
    if ( is_infinite(vh_large) ) {
      vh_small = awvd_small.infinite_vertex();
      vmap[vh_small] = vh_large;
    } else { 
      vh_small = awvd_small.insert(vc->point());
      if ( vh_small.ptr() != NULL ) {
	vmap[vh_small] = vh_large;
      }
    }
    ++vc;
  } while ( vc != vc_start );

  Vertex_handle vn =awvd_small.find_nearest_neighbor(v->point());

  List l(*this);
  Face_map fm;
  Vertex_map vm;
  std::vector<Vh_triple*> flipped_edges;  

  awvd_small.find_conflict_region(v->point(), vn, l, fm, vm,
				  &flipped_edges);

  l.clear();
  fm.clear();
  vm.clear();

  Edge_circulator ec;

  unsigned int num_fe = flipped_edges.size();
  for (unsigned int i = 0; i < num_fe; i++) {
    Vh_triple *vhq = flipped_edges[num_fe - i - 1];

    bool found(false);
    ec = v->incident_edges();
    Edge_circulator ec_start = ec;
    do {
      Edge e = *ec;
      if ( (e.first->vertex(  cw(e.second) ) == vmap[(*vhq)[1]] &&
	    e.first->vertex(     e.second  ) == vmap[(*vhq)[2]]) ||
	   (e.first->vertex( ccw(e.second) ) == vmap[(*vhq)[1]] &&
	    e.first->mirror_vertex(e.second) == vmap[(*vhq)[2]]) ) {
	flip(e);
	found = true;
	break;
      }
      ++ec;
    } while ( ec != ec_start );

    if ( ! found ) {
      Face_circulator fc_start = incident_faces(v);
      Face_circulator fc = fc_start;
      std::cout << "degree: " << v->degree() << std::endl;
      do {
	Face_handle f(fc);
	int j = f->index(v);
	Edge e = sym_edge(f, j);
	Sign s;
	if ( !is_infinite(e) ) {
	  s = awvd_test(e.first, v->point());
	} else {
	  if ( is_infinite(f->vertex(ccw(j))) ) {
	    s = awvd_test(e.first, cw(e.second), v->point());
	  } else {
	    s = awvd_test(e.first, ccw(e.second), v->point());
	  }
	}
	std::cout << "sign: " << int(s) << std::endl;
	++fc;
      } while ( fc != fc_start );
      
      std::ofstream f("data/data_small.out");
      assert(f);

      //      awvd_small.write_non_trivial_weighted_points(f);

      Vertex_circulator vc_start = incident_vertices(v);
      Vertex_circulator vc = vc_start;
      f << v->degree() + 1 << std::endl;

      do {
	Weighted_point q = vc->point();
	f << q.x() << " " << q.y() << " " << q.weight() << std::endl;
	vc++;
      } while (vc!=vc_start);

      Weighted_point q = v->point();
      f << q.x() << " " << q.y() << " " << q.weight() << std::endl;
      f.close();
    }



    CGAL_assertion( found );
  }
  CGAL_triangulation_precondition( v->degree() == 3 );

  _tds.remove_degree_3(v.ptr(), NULL);

  for (unsigned int i = 0; i < num_fe; i++) {
    delete flipped_edges[i];
  }
}

template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_degree_d_vertex_using_conflict_region2(Vertex_handle v)
{
  Additively_weighted_Voronoi_diagram_2< Gt, Tds > awvd_small;

  std::map<Vertex_handle,Vertex_handle> vmap;
  std::map<Vertex_triple, bool> vtmap;

  // first create the triangulation without the vertex to be deleted
  // at the same time create a mapping between vertices of the large
  // and the small Voronoi diagrams
  Vertex_circulator vc_start = v->incident_vertices();
  Vertex_circulator vc = v->incident_vertices();
  Vertex_handle vh_large, vh_small;
  do {
    vh_large = Vertex_handle(vc);
    if ( is_infinite(vh_large) ) {
      vh_small = awvd_small.infinite_vertex();
      vmap[vh_small] = vh_large;
    } else { 
      vh_small = awvd_small.insert(vc->point(), false);
      if ( vh_small.ptr() != NULL ) {
	vmap[vh_small] = vh_large;
      }
    }
    ++vc;
  } while ( vc != vc_start );

  All_faces_iterator fit = awvd_small.all_faces_begin();
  for (; fit != awvd_small.all_faces_end(); ++fit) {
    Vertex_handle v[3];
    for (int i = 0; i < 3; i++) { v[i] = vmap[fit->vertex(i)]; }

    Vertex_triple vt(v[0], v[1], v[2]);
    vtmap[vt] = true;
  }

  // now start flipping edges until we get degree 3
  int degree = v->degree();
  while ( degree > 3 ) {
    Face_circulator fc_start = incident_faces(v);
    Face_circulator fc1 = fc_start, fc2 = fc_start, fc3 = fc_start;
    ++fc2;
    ++(++fc3);

    bool found(false);
    do {
      Vertex_handle v0 = fc1->vertex( ccw(fc1->index(v)) );
      Vertex_handle v1 = fc2->vertex( ccw(fc2->index(v)) );
      Vertex_handle v2 = fc3->vertex( ccw(fc3->index(v)) );    

      Vertex_triple vt(v0, v1, v2);

      if ( vtmap.find(vt) != vtmap.end() ) {
	Face_handle ff(fc2);
	int i = cw(ff->index(v));
	flip(ff, i);
	degree--;
	found = true;
	break;
      }

      fc1++;
      fc2++;
      fc3++;
    } while ( fc1 != fc_start );
    if ( ! found ) {
      Face_circulator fc_start = incident_faces(v);
      Face_circulator fc = fc_start;
      std::cout << "degree: " << v->degree() << std::endl;
      do {
	Face_handle f(fc);
	int j = f->index(v);
	Edge e = sym_edge(f, j);
	Sign s;
	if ( !is_infinite(e) ) {
	  s = awvd_test(e.first, v->point());
	} else {
	  if ( is_infinite(f->vertex(ccw(j))) ) {
	    s = awvd_test(e.first, cw(e.second), v->point());
	  } else {
	    s = awvd_test(e.first, ccw(e.second), v->point());
	  }
	}
	std::cout << "sign: " << int(s) << std::endl;
	++fc;
      } while ( fc != fc_start );
      
      std::ofstream f("data/data_small.out");
      assert(f);

      //      awvd_small.write_non_trivial_weighted_points(f);

      Vertex_circulator vc_start = incident_vertices(v);
      Vertex_circulator vc = vc_start;
      f << v->degree() + 1 << std::endl;

      do {
	Weighted_point q = vc->point();
	f << q.x() << " " << q.y() << " " << q.weight() << std::endl;
	vc++;
      } while (vc!=vc_start);

      Weighted_point q = v->point();
      f << q.x() << " " << q.y() << " " << q.weight() << std::endl;
      f.close();
    }
    CGAL_assertion( found );
  }

  

  CGAL_triangulation_precondition( v->degree() == 3 );

  _tds.remove_degree_3(v.ptr(), NULL);

  awvd_small.clear();
}

#if 0
template< class Gt, class Tds >
void
Additively_weighted_Voronoi_diagram_2< Gt, Tds >::
remove_degree_d_vertex_using_conflict_region2(Vertex_handle v)
{
  Additively_weighted_Voronoi_diagram_2< Gt, Tds > awvd_small;

  std::map<Vertex_handle,Vertex_handle> vmap;
  std::map<Vertex_triple, bool> vtmap;

  // first create the triangulation without the vertex to be deleted
  // at the same time create a mapping between vertices of the large
  // and the small Voronoi diagrams
  Vertex_circulator vc_start = v->incident_vertices();
  Vertex_circulator vc = v->incident_vertices();
  Vertex_handle vh_large, vh_small;
  do {
    vh_large = Vertex_handle(vc);
    if ( is_infinite(vh_large) ) {
      vh_small = awvd_small.infinite_vertex();
      vmap[vh_small] = vh_large;
    } else { 
      vh_small = awvd_small.insert(vc->point(), false);
      if ( vh_small.ptr() != NULL ) {
	vmap[vh_small] = vh_large;
      }
    }
    ++vc;
  } while ( vc != vc_start );

  All_faces_iterator fit = awvd_small.all_faces_begin();
  for (; fit != awvd_small.all_faces_end(); ++fit) {
    Vertex_handle v[3];
    for (int i = 0; i < 3; i++) { v[i] = vmap[fit->vertex(i)]; }

    Vertex_triple vt(v[0], v[1], v[2]);
    vtmap[vt] = true;
  }

  // now start flipping edges until we get degree 3
  int degree = v->degree();
  while ( degree > 3 ) {
    Face_circulator fc_start = incident_faces(v);
    Face_circulator fc1 = fc_start, fc2 = fc_start, fc3 = fc_start;
    ++fc2;
    ++(++fc3);

    bool found(false);
    do {
      Vertex_handle v0 = fc1->vertex( ccw(fc1->index(v)) );
      Vertex_handle v1 = fc2->vertex( ccw(fc2->index(v)) );
      Vertex_handle v2 = fc3->vertex( ccw(fc3->index(v)) );    

      Vertex_triple vt(v0, v1, v2);

      if ( vtmap.find(vt) != vtmap.end() ) {
	Face_handle ff(fc2);
	int i = cw(ff->index(v));
	flip(ff, i);
	degree--;
	found = true;
	break;
      }

      fc1++;
      fc2++;
      fc3++;
    } while ( fc1 != fc_start );
    if ( ! found ) {
      Face_circulator fc_start = incident_faces(v);
      Face_circulator fc = fc_start;
      std::cout << "degree: " << v->degree() << std::endl;
      do {
	Face_handle f(fc);
	int j = f->index(v);
	Edge e = sym_edge(f, j);
	Sign s;
	if ( !is_infinite(e) ) {
	  s = awvd_test(e.first, v->point());
	} else {
	  if ( is_infinite(f->vertex(ccw(j))) ) {
	    s = awvd_test(e.first, cw(e.second), v->point());
	  } else {
	    s = awvd_test(e.first, ccw(e.second), v->point());
	  }
	}
	std::cout << "sign: " << int(s) << std::endl;
	++fc;
      } while ( fc != fc_start );
      
      std::ofstream f("data/data_small.out");
      assert(f);

      //      awvd_small.write_non_trivial_weighted_points(f);

      Vertex_circulator vc_start = incident_vertices(v);
      Vertex_circulator vc = vc_start;
      f << v->degree() + 1 << std::endl;

      do {
	Weighted_point q = vc->point();
	f << q.x() << " " << q.y() << " " << q.weight() << std::endl;
	vc++;
      } while (vc!=vc_start);

      Weighted_point q = v->point();
      f << q.x() << " " << q.y() << " " << q.weight() << std::endl;
      f.close();
    }
    CGAL_assertion( found );
  }

  

  CGAL_triangulation_precondition( v->degree() == 3 );

  _tds.remove_degree_3(v.ptr(), NULL);

  awvd_small.clear();
}
#endif

CGAL_END_NAMESPACE

/*
  Methods with predicates:
   1. insert(const Weighted_point&, bool)
      - compare(nt, nt)
   2. insert(const Weighted_point&, std::vector<Vh_triple*>&)
      - awvd_distance2_test(wp, wp)
      - awvd_test(f, int, wp)
      - awvd_test(f, wp)
   3. insert_second(const Weighted_point&)
      - awvd_distance2_test(wp, wp)
   4. insert_third(const Weighted_point&)
      - awvd_distance2_test(wp, wp)
      - awvd_test_denegerated(wp, wp, wp)
   5. insert(Face_handle&, const Weighted_point&)
      - awvd_distance2_test(wp, wp)
      - awvd_test(f, int, p)
      - awvd_test(f, p)
   6. insert_non_conflicting(Vertex_handle&)
      - awvd_circle_existence_test(wp, wp, wp)
      - awvd_test(wp, wp, wp, wp)
      - awvd_ordering_test(wp, wp, wp)
   7. locate(const Point&)
      - is_inside_face(wp, f)
   8. conflict_test(f, int, wp)
      - awvd_circle_existence_test(f)
      - awvd_test(f, wp)
      - awvd_orientation_test(f, wp, wp)
      - awvd_Voronoi_radii_difference_test(f, int, wp)
      - awvd_test(f, int, wp)
      - awvd_ordering_test(wp, wp, wp)
      - awvd_test_denegerated(wp, wp, wp)
      - orientation(p, p, p)
*/


#endif // CGAL_ADDITIVELY_WEIGHTED_VORONOI_DIAGRAM_2_H

#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_2_H

#include <vector>
#include <map>

#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Segment_Voronoi_diagram_site_2.h>
#include <CGAL/Segment_Voronoi_diagram_data_structure_2.h>
#include <CGAL/Segment_Voronoi_diagram_face_base_2.h>
#include <CGAL/Segment_Voronoi_diagram_vertex_base_2.h>

#include <CGAL/Segment_Voronoi_diagram_constructions_C2.h>

#include <CGAL/in_place_edge_list.h>
#include <CGAL/Segment_Voronoi_diagram_traits_wrapper_2.h>

/*
  Conventions:
  ------------
  1. we treat segments as open; the endpoints are separate objects
  2. a segment of length zero is treated as a point
  3. a point is deleted only if it has no segment adjacent to it
  4. when a segment is deleted it's endpoints are not deleted
  5. the user can force the deletion of endpoints; this is only done
     if condition 3 is met.
  6. when objects are written to a stream we distinguish between
     points and segments; points start by a 'p' and segments by an 's'.
*/


CGAL_BEGIN_NAMESPACE


template < class Gt,
  class Svdds = Segment_Voronoi_diagram_data_structure_2 < 
                Segment_Voronoi_diagram_vertex_base_2<Gt>,
                Segment_Voronoi_diagram_face_base_2<Gt> > >
class Segment_Voronoi_diagram_2
  : private Triangulation_2<
          Segment_Voronoi_diagram_traits_wrapper_2<Gt>, Svdds >
{
private:
  static const char point_descriptor;
  static const char segment_descriptor;

private:
  // types and access methods needed for visualization
  //--------------------------------------------------

  // types
  typedef CGAL::Construct_svd_bisector_2<Gt>
  Construct_svd_bisector_2;

  typedef CGAL::Construct_svd_bisector_ray_2<Gt>
  Construct_svd_bisector_ray_2;

  typedef CGAL::Construct_svd_bisector_segment_2<Gt>
  Construct_svd_bisector_segment_2;

  // access
  Construct_svd_bisector_2
  construct_svd_bisector_2_object() const {
    return Construct_svd_bisector_2();
  }

  Construct_svd_bisector_ray_2
  construct_svd_bisector_ray_2_object() const {
    return Construct_svd_bisector_ray_2();
  }

  Construct_svd_bisector_segment_2
  construct_svd_bisector_segment_2_object() const { 
    return Construct_svd_bisector_segment_2(); 
  }

protected:
  // some local types
  typedef Segment_Voronoi_diagram_traits_wrapper_2<Gt>  Modified_traits;
  typedef Triangulation_2<Modified_traits,Svdds>  DG;
  typedef DG                         Delaunay_graph;
  typedef typename DG::Vertex        Vertex;
  typedef typename DG::Face          Face;

public:
  // TYPES
  //------
  typedef Svdds                                  Data_structure;
  typedef Gt                                     Geom_traits;
  typedef typename Gt::Site_2                    Site_2;
  typedef Site_2                                 Site;
  typedef typename Gt::Point_2                   Point;
  typedef typename Gt::Segment_2                 Segment;

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

protected:
  // some more local types
  //  typedef typename Svdds::Vertex_base          Vertex_base;

  typedef std::map<Face_handle,bool>           Face_map;
  typedef std::map<Face_handle, Face_handle>   Face_face_map;
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

public:
  // CREATION
  //---------
  Segment_Voronoi_diagram_2(const Gt& gt=Gt()) :  DG(gt) {}

  template< class Input_iterator >
  Segment_Voronoi_diagram_2(Input_iterator first, Input_iterator beyond,
			    const Gt& gt=Gt())
    : DG(gt)
  {
    insert(first, beyond);
  }

  Segment_Voronoi_diagram_2(const Segment_Voronoi_diagram_2 &svd)
    : DG(svd)
  {
    CGAL_postcondition( is_valid() );
  }

  Segment_Voronoi_diagram_2&
  operator=(const Segment_Voronoi_diagram_2& svd)
  {
    DG::operator=(svd);
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

  unsigned int number_of_incident_segments(Vertex_handle v) const;


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
  // TRAVERSAL OF THE DUAL GRAPH
  //----------------------------
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

    // do it the obvious way: insert them as they come;
    // one might think though that it might be better to first insert
    // all end points and then all segments, or a variation of that.

    for (Input_iterator it = first; it != beyond; ++it) {
      insert(*it);
    }
  }

  inline Vertex_handle  insert(const Point& p) {
    return insert(Site(p), Vertex_handle(NULL));
  }

  inline Vertex_handle  insert(const Segment& s) {
    return insert(Site(s), Vertex_handle(NULL));
  }

  inline Vertex_handle  insert(const Site& t) {
    return insert(t, Vertex_handle(NULL));
  }

  Vertex_handle  insert(const Site& t, Vertex_handle vnear);

public:
  // REMOVAL
  //--------

  // returns the number of sites removed
  // possible answers:
  // 0 : no site was removed; this can only happen if we ask to
  //     remove a point that is the endpoint of a segment
  // 1 : a single site was removed; this can happen if
  //     (1) we ask to remove a point which is not the endpoint of
  //         a segment
  //     (2) we ask to remove an open segment, i.e., we keep its
  //         endpoints
  //     (3) we ask to remove a closed segment, but its two
  //         endpoints are also endpoints of other segments and
  //         thus cannot be removed
  // 2 : two sites were removed; this can happen if we ask to
  //     remove a closed segment, but only one of its endpoints is
  //     removed; the other is also the endpoint of another
  //     segment
  // 2 : three sites where removed; this can happen when we ask to
  //     remove a closed segment; in this case the two endpoints
  //     are not endpoints of other segment, and thus they are
  //     removed as well
  unsigned int remove(Vertex_handle v,
		      bool remove_endpoints = true);


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
  Stream& write_sites(Stream& str) const
  {
    str << number_of_vertices() << std::endl;

    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      if ( vit->is_point() ) {
	str << "p " << vit->point() << std::endl;
      } else {
	str << "s " << vit->segment() << std::endl;
      }
    }
    return str;
  }

  template < class Stream >
  Stream& draw_sites(Stream &str) const
  {
    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      if ( vit->is_point() ) {
	str << vit->point();
      } else {
	str << vit->segment();
      }
    }
    return str;
  }


  template< class Stream >
  inline
  Stream& draw_primal(Stream &str) const
  {
    Finite_edges_iterator eit = finite_edges_begin();
    for (; eit != finite_edges_end(); ++eit) {
      draw_voronoi_edge(*eit, str);
    }
    return str;
  }

  template < class Stream > 
  Stream& draw_skeleton(Stream &str) const
  {
    Finite_edges_iterator eit = finite_edges_begin();
    for (; eit != finite_edges_end(); ++eit) {
      Site p = eit->first->vertex(  cw(eit->second) )->site();
      Site q = eit->first->vertex( ccw(eit->second) )->site();

      bool is_endpoint_of_seg =
	( p.is_segment() && q.is_point() &&
	  is_endpoint_of_segment(q.point(), p.segment()) ) ||
	( p.is_point() && q.is_segment() &&
	  is_endpoint_of_segment(p.point(), q.segment()) );

      if ( !is_endpoint_of_seg ) {
	draw_voronoi_edge(*eit, str);
      }
    }
    return str;
  }

  template < class Stream >
  Stream& draw_Voronoi_circles(Stream& str) const
  {
    Finite_faces_iterator fit = finite_faces_begin();
    for (; fit != finite_faces_end(); ++fit) {
      typename Gt::Circle_2 c = circumcircle(Face_handle(fit));
      str << c;
      str << c.center();
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

  inline void swap(const Segment_Voronoi_diagram_2& svd) {
    DG::swap(svd);
  }

public:
  // MK: THE FOLLOWING ARE NOT IN THE SPEC
  //======================================
  // Primal
  Point  primal(const Face_handle& f) const;
  Object primal(const Edge e) const;
  inline Object primal(const Edge_circulator& ec) const {
    return primal(*ec);
  }
  inline Object primal(const Finite_edges_iterator& ei) const {
    return primal(*ei);
  }


protected:
  // wrappers for the geometric predicates

  bool are_same_points(const Point& p, const Point& q) const;

  inline
  bool is_endpoint_of_segment(const Point& p, const Segment& s) const
  {
    return ( are_same_points(p, s.source()) ||
	     are_same_points(p, s.target()) );
  }

  inline bool is_degenerate_segment(const Segment& s) const {
    return are_same_points(s.source(), s.target());
  }

  // returns:
  //   ON_POSITIVE_SIDE if q is closer to t1
  //   ON_NEGATIVE_SIDE if q is closer to t2
  //   ON_ORIENTED_BOUNDARY if q is on the bisector of t1 and t2
  Oriented_side side_of_bisector(const Site &t1,
				 const Site &t2,
				 const Point &q) const;

  Sign incircle(const Site &t1, const Site &t2,
		const Site &t3, const Site &q) const;

  Sign incircle(const Site &t1, const Site &t2,
		const Site &q) const;


  Sign incircle(const Face_handle& f, const Site& q) const;


  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v) const;

  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Vertex_handle& v) const;


  
  bool finite_edge_interior(const Site& t1, const Site& t2,
			    const Site& t3, const Site& t4,
			    const Site& q,  Sign sgn) const;

  bool finite_edge_interior(const Face_handle& f, int i,
			    const Site& q, Sign sgn) const;

  bool finite_edge_interior(const Vertex_handle& v1,
			    const Vertex_handle& v2,
			    const Vertex_handle& v3,
			    const Vertex_handle& v4,
			    const Vertex_handle& v,
			    Sign sgn) const;

  bool finite_edge_interior_degenerated(const Site& t1,	const Site& t2,
					const Site& t3,	const Site& q,
					Sign sgn) const;


  bool finite_edge_interior_degenerated(const Site& t1,	const Site& t2,
					const Site& q,	Sign sgn) const;

  bool finite_edge_interior_degenerated(const Face_handle& f, int i,
					const Site& p, Sign sgn) const;

  bool finite_edge_interior_degenerated(const Vertex_handle& v1,
					const Vertex_handle& v2,
					const Vertex_handle& v3,
					const Vertex_handle& v4,
					const Vertex_handle& v,
					Sign Sign) const;

  bool infinite_edge_interior(const Site& t2, const Site& t3,
			      const Site& t4, const Site& q,
			      Sign sgn) const;


  bool infinite_edge_interior(const Face_handle& f, int i,
			      const Site& q, Sign sgn) const;

  bool infinite_edge_interior(const Vertex_handle& v1,
			      const Vertex_handle& v2,
			      const Vertex_handle& v3,
			      const Vertex_handle& v4,
			      const Vertex_handle& v,
			      Sign sgn) const;

  Conflict_type
  finite_edge_conflict_type_degenerated(const Site& t1,
					const Site& t2,
					const Site& t) const;

  bool edge_interior(const Face_handle& f, int i,
		     const Site& t, Sign sgn) const;


  inline bool edge_interior(const Edge& e,
			    const Site& t, Sign sgn) const {
    return edge_interior(e.first, e.second, t, sgn);
  }

  bool edge_interior(const Vertex_handle& v1,
		     const Vertex_handle& v2,
		     const Vertex_handle& v3,
		     const Vertex_handle& v4,
		     const Vertex_handle& v,
		     Sign sgn) const;

  inline bool is_degenerate_edge(const Site& t1,
				 const Site& t2,
				 const Site& t3,
				 const Site& t4) const {
    return geom_traits().is_degenerate_edge_2_object()
      (t1, t2, t3, t4);
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


  //protected:
public:
  // wrappers for constructions
  Point circumcenter(const Face_handle& f) const;
  Point circumcenter(const Site& t0, 
		     const Site& t1, 
		     const Site& t2) const;

  typename Gt::Circle_2 circumcircle(const Face_handle& f) const;
  typename Gt::Circle_2 circumcircle(const Site& t0, const Site& t1, 
				     const Site& t2) const;

  typename Gt::Line_2 circumcircle(const Point& p0, const Point& p1) const;

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

  //  Vertex_handle insert_in_face(Face_handle& f, const Weighted_point& p);

  bool          is_degree_2(const Vertex_handle& v) const;

  Vertex_handle insert_degree_2(Edge e);
  Vertex_handle insert_degree_2(Edge e, const Site& t);
  void          remove_degree_2(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v, Face* f);

  // this was defined because the hierarchy needs it
  inline Vertex_handle create_vertex() {
    return _tds.create_vertex();
  }


protected:
  // insertion of the first three sites

  // the first two objects can only be points, since we always
  // add the endpoints first and then the segment.
  Vertex_handle  insert_first(const Point& p);
  Vertex_handle  insert_second(const Point& p);
  Vertex_handle  insert_third(const Site& t);

  // methods for insertion
  void initialize_conflict_region(const Face_handle& f, List& l);

  void expand_conflict_region(const Face_handle& f, const Site& t,
			      List& l, Face_map& fm,
			      std::map<Face_handle,Sign>& sign_map,
			      std::vector<Vh_triple*>* fe);

  Vertex_handle add_bogus_vertex(Edge e, List& l);
  Vertex_list   add_bogus_vertices(List& l);
  void          remove_bogus_vertices(Vertex_list& vl);

  // MK: this is not currently used
  inline  std::vector<Face*> get_faces_for_recycling(Face_map& fm,
					     unsigned int n_wanted);

  Vertex_handle retriangulate_conflict_region(const Site& t, List& l, 
					      Face_map& fm);

protected:
  // methods for removal

  std::pair<Vertex_handle,Vertex_handle >
  endpoint_vertices(Vertex_handle v) const;

  bool is_endpoint_of_segment(Vertex_handle v) const;

  void  remove_first(Vertex_handle v);
  void  remove_second(Vertex_handle v);
  unsigned int remove_third(Vertex_handle v, bool remove_endpoints);
  unsigned int remove_degree_2(Vertex_handle v,
			       bool remove_endpoints);
  unsigned int remove_degree_3(Vertex_handle v,
			       bool remove_endpoints);
  unsigned int remove_degree_d(Vertex_handle v,
			       bool remove_endpoints);

  void  minimize_degree(Vertex_handle v);

  void find_conflict_region_remove(const Vertex_handle& v,
				   const Vertex_handle& vnearest,
				   List& l, Face_map& fm,
				   std::vector<Vh_triple*>* fe);

protected:
  // methods for I/O

  // MK: this has to be rewritten. all the checking must be done in
  // the geometric traits class.

  template< class Stream >
  Stream& draw_voronoi_edge(Edge e, Stream &str) const
  {
    typename Geom_traits::Line_2              l;
    typename Geom_traits::Segment_2           s;
    typename Geom_traits::Ray_2               r;
    CGAL::Parabola_segment_2<Gt>              ps;

    Object o = primal(e);

    if (CGAL::assign(l, o))   str << l;
    if (CGAL::assign(s, o))   str << s; 
    if (CGAL::assign(r, o))   str << r;
    if (CGAL::assign(ps, o))  str << ps;

    return str;
  }


}; // Segment_Voronoi_diagram_2


CGAL_END_NAMESPACE


#include <CGAL/Segment_Voronoi_diagram_2.C>



#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_2_H

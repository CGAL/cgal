#ifndef CGAL_APOLLONIUS_GRAPH_2_TEST_H
#define CGAL_APOLLONIUS_GRAPH_2_TEST_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>

#include <CGAL/Vector_2.h> // this is done in order to avoid error
// when the  Segment_2_Segment_2_intersection.h file is included from
// the Triangulation_euclidean_traits_2.h file.

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_hierarchy_2.h>
#include <CGAL/Apollonius_graph_euclidean_traits_2.h>



CGAL_BEGIN_NAMESPACE


//template<class Kernel, class Method_tag, class InputStream>
template<class Kernel, class Method_tag>
bool test_traits()
  //bool test_traits(InputStream& is)
{
  typedef Apollonius_graph_euclidean_traits_2<Kernel,Method_tag>  Traits;

  // testing typedefs
  //--------------------------------------------------------------------
  typedef typename Traits::Bare_point            Bare_point;
  typedef typename Traits::Weighted_point        Weighted_point;
  typedef typename Traits::Site                  Site;
  typedef typename Traits::Line_2                Line_2;
  typedef typename Traits::Ray_2                 Ray_2;
  typedef typename Traits::Segment_2             Segment_2;
  typedef typename Traits::Parabola_segment_2    Parabola_segment_2;
  typedef typename Traits::Hyperbola_2           Hyperbola_2;
  typedef typename Traits::Hyperbola_ray_2       Hyperbola_ray_2;
  typedef typename Traits::Hyperbola_segment_2   Hyperbola_segment_2;

  typedef typename Traits::Construct_Apollonius_vertex_2
    Construct_Apollonius_vertex_2;
  typedef typename Traits::Construct_Apollonius_weighted_point_2
    Construct_Apollonius_weighted_point_2;
  typedef typename Traits::Construct_Apollonius_bisector_2
    Construct_Apollonius_bisector_2;
  typedef typename Traits::Construct_Apollonius_bisector_ray_2
    Construct_Apollonius_bisector_ray_2;
  typedef typename Traits::Construct_Apollonius_bisector_segment_2
    Construct_Apollonius_bisector_segment_2;
  typedef typename Traits::Construct_Apollonius_primal_ray_2
    Construct_Apollonius_primal_ray_2;
  typedef typename Traits::Construct_Apollonius_primal_segment_2
    Construct_Apollonius_primal_segment_2;

  typedef typename Traits::Compare_x_2       Compare_x_2;
  typedef typename Traits::Compare_y_2       Compare_y_2;
  typedef typename Traits::Compare_weight_2  Compare_weight_2;
  typedef typename Traits::Orientation_2     Orientation_2;
  typedef typename Traits::Is_hidden_2       Is_hidden_2;
  typedef typename Traits::Oriented_side_of_bisector_2
    Oriented_side_of_bisector_2;
  typedef typename Traits::Vertex_conflict_2 Vertex_conflict_2;
  typedef typename Traits::Finite_edge_interior_conflict_2
    Finite_edge_interior_conflict_2;
  typedef typename Traits::Infinite_edge_interior_conflict_2
    Infinite_edge_interior_conflict_2;
  typedef typename Traits::Is_degenerate_edge_2
    Is_degenerate_edge_2;


  // testing constructors
  //--------------------------------------------------------------------
  Traits tr;
  Traits tr1(tr);
  tr1 = tr;
  tr = tr1;


  // testing access to predicates objects
  //--------------------------------------------------------------------
  Compare_x_2 compare_x = tr.compare_x_2_object();
  Compare_y_2 compare_y = tr.compare_y_2_object();
  Compare_weight_2 compare_w = tr.compare_weight_2_object();
  Orientation_2 orienation = tr.orientation_2_object();
  Is_hidden_2 is_hidden = tr.is_hidden_2_object();
  Oriented_side_of_bisector_2 oriented_side_of_bisector =
    tr.oriented_side_of_bisector_2_object();
  Vertex_conflict_2 vertex_conflict = tr.vertex_conflict_2_object();
  Finite_edge_interior_conflict_2 finite_edge_interior_conflict =
    tr.finite_edge_interior_conflict_2_object();
  Infinite_edge_interior_conflict_2 infinite_edge_interior_conflict =
    tr.infinite_edge_interior_conflict_2_object();
  Is_degenerate_edge_2 is_degenerate_edge =
    tr.is_degenerate_edge_2_object();

  // testing access to constructor objects
  //--------------------------------------------------------------------
  Construct_Apollonius_vertex_2 apollonius_vertex =
    tr.construct_Apollonius_vertex_2_object();
  Construct_Apollonius_weighted_point_2 apollonius_weighted_point =
    tr.construct_Apollonius_weighted_point_2_object();
  Construct_Apollonius_bisector_2 apollonius_bisector =
    tr.construct_Apollonius_bisector_2_object();
  Construct_Apollonius_bisector_ray_2 apollonius_bisector_ray =
    tr.construct_Apollonius_bisector_ray_2_object();
  Construct_Apollonius_bisector_segment_2 apollonius_bisector_segment =
    tr.construct_Apollonius_bisector_segment_2_object();
  Construct_Apollonius_primal_ray_2 apollonius_primal_ray =
    tr.construct_Apollonius_primal_ray_2_object();
  Construct_Apollonius_primal_segment_2 apollonius_primal_segment =
    tr.construct_Apollonius_primal_segment_2_object();



  // testing correctness of predicates;
  //--------------------------------------------------------------------
  bool b;
  Weighted_point
    wp1(Bare_point(10.0,-10),20),
    wp2(Bare_point(9,-9),19.),
    wp3(Bare_point(1000000,-1000000),0),
    wp4(Bare_point(10.0,5.0),5),
    wp5;

  // testing coordinate comparison
  //--------------------------------------------------------------------
  CGAL_assertion( compare_x(wp1, wp1) == EQUAL );
  CGAL_assertion( compare_y(wp1, wp1) == EQUAL );
  CGAL_assertion( compare_w(wp1, wp1) == EQUAL );

  CGAL_assertion( compare_x(wp1, wp2) == LARGER );
  CGAL_assertion( compare_y(wp1, wp2) == SMALLER );
  CGAL_assertion( compare_w(wp1, wp2) == LARGER );

  // testing orientation
  //--------------------------------------------------------------------
  CGAL_assertion( orientation(wp1, wp2, wp3) == COLLINEAR );

  // testing is_hidden
  //--------------------------------------------------------------------
  CGAL_assertion( is_hidden(wp1, wp4) == true );
  CGAL_assertion( is_hidden(wp4, wp1) == false );

  // testing oriented_side_of_bisector
  //--------------------------------------------------------------------
  Bare_point p = wp2.point();
  CGAL_assertion( oriented_side_of_bisector(wp1, wp3, p) ==
		  ON_POSITIVE_SIDE );


  // testing vertex_conflict
  //--------------------------------------------------------------------

  // first we consider the case where all vertices are finite
  wp1 = Weighted_point(Bare_point(0,0),0);
  wp2 = Weighted_point(Bare_point(1,1),1);
  wp3 = Weighted_point(Bare_point(2,4),4);
  wp4 = Weighted_point(Bare_point(3,9),9);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp4) == ZERO );
  CGAL_assertion( vertex_conflict(wp2, wp3, wp4, wp1) == ZERO );

  wp4 = Weighted_point(Bare_point(3,9),9.0001);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp4) == NEGATIVE );
  CGAL_assertion( vertex_conflict(wp2, wp3, wp4, wp1) == POSITIVE );


  // then we consider the case where v3 is the vertex at infinity
  wp1 = Weighted_point(Bare_point(1,100),49);
  wp2 = Weighted_point(Bare_point(0,-100),50);
  
  wp3 = Weighted_point(Bare_point(0,0),0);
  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == POSITIVE );

  wp3 = Weighted_point(Bare_point(-100000,0),0);
  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == POSITIVE );

  wp3 = Weighted_point(Bare_point(100000,0),1000);
  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == NEGATIVE );

  wp3 = Weighted_point(Bare_point(-1,0),51);
  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == NEGATIVE );

  wp3 = Weighted_point(Bare_point(-1,200),51);
  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == POSITIVE );

  wp3 = Weighted_point(Bare_point(-1,-200),51);
  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == POSITIVE );


  // testing finite_edge_interior_conflict
  //--------------------------------------------------------------------

  wp1 = Weighted_point(Bare_point(0,-100),50);
  wp2 = Weighted_point(Bare_point(0,100),51);
  wp3 = Weighted_point(Bare_point(-150,0),49);
  wp4 = Weighted_point(Bare_point(150,0),48);


  // first we look at the case where all vertices of the edge are
  // finite...

  // endpoints are not in conflict but the interior is
  wp5 = Weighted_point(Bare_point(0,-49),0);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp5) == POSITIVE );
  CGAL_assertion( vertex_conflict(wp1, wp4, wp2, wp5) == POSITIVE );
  
  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, wp5,  false);

  CGAL_assertion( b );

  // both endpoints and interior are in conflict
  wp5 = Weighted_point(Bare_point(0,-49),10);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp5) == NEGATIVE );
  CGAL_assertion( vertex_conflict(wp1, wp4, wp2, wp5) == NEGATIVE );
  
  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, wp5,  true);

  CGAL_assertion( b );

  // endpoints are in conflict but the interior isn't
  wp5 = Weighted_point(Bare_point(0,-150),99);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp5) == NEGATIVE );
  CGAL_assertion( vertex_conflict(wp1, wp4, wp2, wp5) == NEGATIVE );
  
  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, wp5,  true);

  CGAL_assertion( !b );

  // finally endpoints are not in conflict and interior is not in
  // conflict
  wp5 = Weighted_point(Bare_point(0,-150),70);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp5) == POSITIVE );
  CGAL_assertion( vertex_conflict(wp1, wp4, wp2, wp5) == POSITIVE );
  
  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, wp5,  true);

  CGAL_assertion( !b );


  // then we consider the case where v4 is the vertex at infinity
  wp1 = Weighted_point(Bare_point(0,-100),50);
  wp2 = Weighted_point(Bare_point(1,100),49);
  wp3 = Weighted_point(Bare_point(-150,0),51);

  // endpoints are not in conflict but interior is
  wp4 = Weighted_point(Bare_point(0,-48),0);
  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp4) == POSITIVE );
  CGAL_assertion( vertex_conflict(wp2, wp1, wp4) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, false);
  CGAL_assertion( b );

  // both endpoints and interior are in conflict
  wp4 = Weighted_point(Bare_point(-1,0),51);
  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp4) == NEGATIVE );
  CGAL_assertion( vertex_conflict(wp2, wp1, wp4) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, true);
  CGAL_assertion( b );

  // endpoints are in conflict but interior isn't
  wp4 = Weighted_point(Bare_point(-1,200),149);
  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp4) == NEGATIVE );
  CGAL_assertion( vertex_conflict(wp2, wp1, wp4) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, true);
  CGAL_assertion( !b );

  // neither the endpoints nor the interior are in conflict
  wp4 = Weighted_point(Bare_point(-1,200),51);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3, wp4) == POSITIVE );
  CGAL_assertion( vertex_conflict(wp2, wp1, wp4) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, false);
  CGAL_assertion( !b );


  // finally consider the case where both v3 and v4 are the vertex at
  // infinity.

  wp1 = Weighted_point(Bare_point(10,0),5);
  wp2 = Weighted_point(Bare_point(100,0),50);

  // endpoints are not in conflict but the interior is
  wp3 = Weighted_point(Bare_point(20,0),0);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == POSITIVE );
  CGAL_assertion( vertex_conflict(wp2, wp1, wp3) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, false);
  CGAL_assertion( b );

  // endpoints are in conflict and so is the interior
  wp3 = Weighted_point(Bare_point(20,0),10);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == NEGATIVE );
  CGAL_assertion( vertex_conflict(wp2, wp1, wp3) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, true);
  CGAL_assertion( b );

  // both the endpoints and the interior are not in conflict
  wp3 = Weighted_point(Bare_point(4,0),2);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == POSITIVE );
  CGAL_assertion( vertex_conflict(wp2, wp1, wp3) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, false);
  CGAL_assertion( !b );

  // the endpoints are in conflict but the interior isn't
  wp3 = Weighted_point(Bare_point(4,0),3);

  CGAL_assertion( vertex_conflict(wp1, wp2, wp3) == NEGATIVE );
  CGAL_assertion( vertex_conflict(wp2, wp1, wp3) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, true);
  CGAL_assertion( !b );

  // testing infinite_edge_interior_conflict
  //--------------------------------------------------------------------

  
  wp2 = Weighted_point(Bare_point(0,0),100);
  wp3 = Weighted_point(Bare_point(-100,200),5);
  wp4 = Weighted_point(Bare_point(100,300),4);

  // the endpoints are not in conflict but the interior is
  wp5 = Weighted_point(Bare_point(0,-150),2);
  
  CGAL_assertion( vertex_conflict(wp2, wp3, wp5) == POSITIVE );
  CGAL_assertion( vertex_conflict(wp4, wp2, wp5) == POSITIVE );

  b = infinite_edge_interior_conflict(wp2, wp3, wp4, wp5, false);
  CGAL_assertion( b );


  // the endpoints are in conflict but the interior isn't
  wp5 = Weighted_point(Bare_point(0,150),150);
  
  CGAL_assertion( vertex_conflict(wp2, wp3, wp5) == NEGATIVE );
  CGAL_assertion( vertex_conflict(wp4, wp2, wp5) == NEGATIVE );

  b = infinite_edge_interior_conflict(wp2, wp3, wp4, wp5, true);
  CGAL_assertion( !b );


  // the endpoints are in conflict as well as the interior
  wp5 = Weighted_point(Bare_point(0,-150),150);
  
  CGAL_assertion( vertex_conflict(wp2, wp3, wp5) == NEGATIVE );
  CGAL_assertion( vertex_conflict(wp4, wp2, wp5) == NEGATIVE );

  b = infinite_edge_interior_conflict(wp2, wp3, wp4, wp5, true);
  CGAL_assertion( b );


  // neither the endpoints nor the interior are in conflict
  wp5 = Weighted_point(Bare_point(0,150),50);
  
  CGAL_assertion( vertex_conflict(wp2, wp3, wp5) == POSITIVE );
  CGAL_assertion( vertex_conflict(wp4, wp2, wp5) == POSITIVE );

  b = infinite_edge_interior_conflict(wp2, wp3, wp4, wp5, false);
  CGAL_assertion( !b );



  // testing is_degenerate_edge
  //--------------------------------------------------------------------
  wp1 = Weighted_point(Bare_point(2,4),4);
  wp2 = Weighted_point(Bare_point(10,100),100);
  wp3 = Weighted_point(Bare_point(5,25),25);
  wp4 = Weighted_point(Bare_point(20,400),400);

  b = is_degenerate_edge(wp1, wp2, wp3, wp4);

  CGAL_assertion( b );

  wp1 = Weighted_point(Bare_point(2,4),4);
  wp2 = Weighted_point(Bare_point(10,100),100);
  wp3 = Weighted_point(Bare_point(5,25),25);
  wp4 = Weighted_point(Bare_point(20,400),399);

  b = is_degenerate_edge(wp1, wp2, wp3, wp4);

  CGAL_assertion( !b );

  return true;
}

template<class Kernel, class Method_tag, class InputStream>
bool test_algo(InputStream& is)
{
  typedef Apollonius_graph_euclidean_traits_2<Kernel,Method_tag> Traits;
  typedef Apollonius_graph_2<Traits>  Apollonius_graph;

  // testing typedefs
  //--------------------------------------------------------------------
  typedef typename Apollonius_graph::Data_structure   Data_structure;
  typedef typename Apollonius_graph::Geom_traits      Geom_traits;
  typedef typename Apollonius_graph::Point            Point;
  typedef typename Apollonius_graph::Weighted_point   Weighted_point;
  typedef typename Apollonius_graph::Weight           Weight;

  typedef typename Apollonius_graph::Edge             Edge;
  typedef typename Apollonius_graph::Vertex_handle    Vertex_handle;
  typedef typename Apollonius_graph::Face_handle      Face_handle;

  typedef typename Apollonius_graph::Edge_circulator   Edge_circulator;
  typedef typename Apollonius_graph::Vertex_circulator Vertex_circulator;
  typedef typename Apollonius_graph::Face_circulator   Face_circulator;

  typedef typename Apollonius_graph::All_vertices_iterator
    All_vertices_iterator;
  typedef typename Apollonius_graph::Finite_vertices_iterator
    Finite_vertices_iterator;

  typedef typename Apollonius_graph::All_faces_iterator
    All_faces_iterator;
  typedef typename Apollonius_graph::Finite_faces_iterator
    Finite_faces_iterator;

  typedef typename Apollonius_graph::All_edges_iterator
    All_edges_iterator;
  typedef typename Apollonius_graph::Finite_edges_iterator
    Finite_edges_iterator;

  // testing creation/constructors
  //--------------------------------------------------------------------

  Apollonius_graph ag;
  Apollonius_graph ag1(Traits());
  Apollonius_graph ag2(ag);

  std::vector<Weighted_point> wp_list;
  Weighted_point wp;
  while ( is >> wp ) {
    wp_list.push_back(wp);
  }

  Apollonius_graph ag3(wp_list.begin(), wp_list.end(), Traits());

  ag = ag3;
  ag2 = ag3;

  // testing access functions
  //--------------------------------------------------------------------
  Geom_traits tr = ag.geom_traits();
  int num_vertices = ag.number_of_vertices();
  CGAL_assertion( static_cast<unsigned int>(num_vertices) == wp_list.size() );

  Face_handle inf_f = ag.infinite_face();
  Vertex_handle v1 = ag.infinite_vertex();
  Vertex_handle v2 = ag.finite_vertex();

  // testing traversal
  //--------------------------------------------------------------------

  // finite faces, edges and vertices
  Finite_vertices_iterator fvit = ag.finite_vertices_begin();
  int n_fvertices = 0;
  for (; fvit != ag.finite_vertices_end(); ++fvit) {
    n_fvertices++;
  }
  CGAL_assertion( n_fvertices == num_vertices );

  Finite_edges_iterator feit = ag.finite_edges_begin();
  int n_fedges = 0;
  for (; feit != ag.finite_edges_end(); ++feit) {
    n_fedges++;
  }

  Finite_faces_iterator ffit = ag.finite_faces_begin();
  int n_ffaces = 0;
  for (; ffit != ag.finite_faces_end(); ++ffit) {
    n_ffaces++;
  }

  // all faces, edges and vertices
  All_vertices_iterator avit = ag.all_vertices_begin();
  int n_avertices = 0;
  for (; avit != ag.all_vertices_end(); ++avit) {
    n_avertices++;
  }
  CGAL_assertion( n_avertices == num_vertices + 1 );

  All_edges_iterator aeit = ag.all_edges_begin();
  int n_aedges = 0;
  for (; aeit != ag.all_edges_end(); ++aeit) {
    n_aedges++;
  }

  All_faces_iterator afit = ag.all_faces_begin();
  int n_afaces = 0;
  for (; afit != ag.all_faces_end(); ++afit) {
    n_afaces++;
  }

    
  CGAL_assertion( 2 * n_aedges == 3 * n_afaces );
  CGAL_assertion( n_avertices - n_aedges + n_afaces == 2 );

  // vertex circulators
  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Vertex_circulator vc_start = ag.incident_vertices(v);
    Vertex_circulator vc = vc_start;
    int deg = 0;
    do {
      deg++;
      vc++;
    } while ( vc != vc_start );

    CGAL_assertion( deg == v->degree() );
  }

  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Vertex_circulator vc_start = ag.incident_vertices(v, v->face());
    Vertex_circulator vc = vc_start;
    int deg = 0;
    do {
      deg++;
      vc++;
    } while ( vc != vc_start );

    CGAL_assertion( deg == v->degree() );
  }

  // face circulators
  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Face_circulator fc_start = ag.incident_faces(v);
    Face_circulator fc = fc_start;
    int deg = 0;
    do {
      deg++;
      fc++;
    } while ( fc != fc_start );

    CGAL_assertion( deg == v->degree() );
  }

  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Face_circulator fc_start = ag.incident_faces(v, v->face());
    Face_circulator fc = fc_start;
    int deg = 0;
    do {
      deg++;
      fc++;
    } while ( fc != fc_start );

    CGAL_assertion( deg == v->degree() );
  }

  // edge circulators
  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Edge_circulator ec_start = ag.incident_edges(v);
    Edge_circulator ec = ec_start;
    int deg = 0;
    do {
      deg++;
      ec++;
    } while ( ec != ec_start );

    CGAL_assertion( deg == v->degree() );
  }

  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Edge_circulator ec_start = ag.incident_edges(v, v->face());
    Edge_circulator ec = ec_start;
    int deg = 0;
    do {
      deg++;
      ec++;
    } while ( ec != ec_start );

    CGAL_assertion( deg == v->degree() );
  }


  // testing predicates
  //--------------------------------------------------------------------
  CGAL_assertion( ag.is_infinite(inf_f) );
  CGAL_assertion( ag.is_infinite(v1) );
  CGAL_assertion( !ag.is_infinite(v2) );
  {
    Edge_circulator ec = ag.incident_edges(ag.infinite_vertex());
    CGAL_assertion( ag.is_infinite(*ec) );
    CGAL_assertion( ag.is_infinite(ec) );
    CGAL_assertion( ag.is_infinite(ec->first, ec->second) );
  }

  // testing insertion
  //--------------------------------------------------------------------
  ag.clear();
  ag.insert(wp_list.begin(), wp_list.end());
  CGAL_assertion( ag.is_valid() );

  ag.clear();
  typename std::vector<Weighted_point>::iterator it;
  for (it = wp_list.begin(); it != wp_list.end(); ++it) {
    ag.insert(*it);
  }
  CGAL_assertion( ag.is_valid() );

  {
    ag.clear();
    Vertex_handle v;
    typename std::vector<Weighted_point>::iterator it;
    for (it = wp_list.begin(); it != wp_list.end(); ++it) {
      if ( it == wp_list.begin() ) {
	v = ag.insert(*it);
      } else {
	v = ag.insert(*it, v);
      }
    }
    CGAL_assertion( ag.is_valid() );
  }

  // testing removal
  //--------------------------------------------------------------------
  ag.clear();
  ag.insert(wp_list.begin(), wp_list.end());
  CGAL_assertion( ag.is_valid() );
  {
    while ( ag.number_of_vertices() > 0 ) {
      Vertex_handle v(ag.finite_vertices_begin());
      ag.remove(v);
      CGAL_assertion( ag.is_valid() );
    }
  }

  // testing nearest neighbor location
  //--------------------------------------------------------------------
  ag.clear();
  ag.insert(wp_list.begin(), wp_list.end());
  CGAL_assertion( ag.is_valid() );
  for (fvit = ag.finite_vertices_begin();
       fvit != ag.finite_vertices_end(); ++fvit) {
    Weighted_point wp = fvit->point();
    Vertex_handle nn = ag.nearest_neighbor(wp);
    CGAL_assertion( wp == nn->point() );
  }

  // testing swap
  //--------------------------------------------------------------------
  ag.clear();
  ag.swap(ag2);
  CGAL_assertion( ag.number_of_vertices() ==
		  static_cast<int>(wp_list.size()) );
  CGAL_assertion( ag2.number_of_vertices() == 0 );

  return true;
}


CGAL_END_NAMESPACE



#endif // CGAL_APOLLONIUS_GRAPH_2_TEST_H

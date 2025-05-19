#ifndef CGAL_APOLLONIUS_GRAPH_2_TEST_H
#define CGAL_APOLLONIUS_GRAPH_2_TEST_H

#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_hierarchy_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>

#include <cassert>
#include <CGAL/enum.h>
#include <CGAL/use.h>
#include <CGAL/Random.h>


#include "IO/Null_output_stream.h"


namespace CGAL {

bool dummy(bool b) {
  assert( b );
  return b;
}

#if defined(__INTEL_COMPILER)
template<class Traits>
bool test_traits_base(const Traits& = Traits());

template<class AG, class InputStream>
bool test_algo_generic(InputStream& is, const AG& = AG());

template<class Kernel, class Method_tag, class InputStream>
bool test_algo(InputStream& is, const Kernel& = Kernel(),
               const Method_tag& = Method_tag());

template<class Kernel, class Method_tag, class InputStream>
bool test_hierarchy_algo(InputStream& is, const Kernel& = Kernel(),
                         const Method_tag& = Method_tag());

template<class CK, class CKM, class EK, class EKM, class InputStream>
bool test_filtered_traits_algo(InputStream& is,
                               const CK& = CK(), const CKM& = CKM(),
                               const EK& = EK(),
                               const EKM& = EKM());

template<class CK, class CKM, class EK, class EKM, class InputStream>
bool test_filtered_traits_hierarchy_algo(InputStream& is,
                                         const CK& = CK(),
                                         const CKM& = CKM(),
                                         const EK& = EK(),
                                         const EKM& = EKM());
#endif


template<class Kernel, class Method_tag>
struct Traits_tester
{
  typedef Apollonius_graph_traits_2<Kernel,Method_tag>  Traits;

  bool operator()(int = 0) const {
    return test_traits_base( Traits() );
  }
};

template<class CK, class CK_Method, class EK, class EK_Method>
struct Filtered_traits_tester
{
  typedef
  Apollonius_graph_filtered_traits_2<CK, CK_Method, EK, EK_Method>
  Traits;

  bool operator()(int = 0) const {
    return test_traits_base( Traits() );
  }
};


//template<class Kernel, class Method_tag, class InputStream>
#if defined(__INTEL_COMPILER)
template<class Traits>
bool test_traits_base(const Traits&)
#else
template<class Traits>
bool test_traits_base(const Traits& = Traits())
#endif
  //bool test_traits(InputStream& is)
{
  //  typedef Apollonius_graph_traits_2<Kernel,Method_tag>  Traits;

  // testing typedefs
  //--------------------------------------------------------------------
  typedef typename Traits::FT                    FT;
  typedef typename Traits::RT                    RT;
  typedef typename Traits::Point_2               Point_2;
  typedef typename Traits::Site_2                Site_2;
  typedef typename Traits::Object_2              Object_2;
  typedef typename Traits::Line_2                Line_2;
  typedef typename Traits::Ray_2                 Ray_2;
  typedef typename Traits::Segment_2             Segment_2;
  //  typedef typename Traits::Parabola_segment_2    Parabola_segment_2;
  //  typedef typename Traits::Hyperbola_2           Hyperbola_2;
  //  typedef typename Traits::Hyperbola_ray_2       Hyperbola_ray_2;
  //  typedef typename Traits::Hyperbola_segment_2   Hyperbola_segment_2;

  typedef typename Traits::Construct_object_2    Construct_object_2;
  typedef typename Traits::Assign_2              Assign_2;

  typedef typename Traits::Construct_Apollonius_vertex_2
    Construct_Apollonius_vertex_2;
  typedef typename Traits::Construct_Apollonius_site_2
    Construct_Apollonius_site_2;

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

  CGAL_USE_TYPE(FT);
  CGAL_USE_TYPE(RT);
  CGAL_USE_TYPE(Object_2);
  CGAL_USE_TYPE(Line_2);
  CGAL_USE_TYPE(Ray_2);
  CGAL_USE_TYPE(Segment_2);

  // testing constructors
  //--------------------------------------------------------------------
  Traits tr;
  //  tr = Traits(); // to avoid compiler-warning that variable was not initialized;
  Traits tr1(tr);
  tr1 = tr;
  tr = tr1;


  // testing access to predicates objects
  //--------------------------------------------------------------------
  Compare_x_2 compare_x = tr.compare_x_2_object();
  Compare_y_2 compare_y = tr.compare_y_2_object();
  Compare_weight_2 compare_w = tr.compare_weight_2_object();
  Orientation_2 orientation = tr.orientation_2_object();
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
  Assign_2 assign_2 = tr.assign_2_object();
  CGAL_USE(assign_2);

  Construct_object_2 construct_object_2 =
    tr.construct_object_2_object();
  CGAL_USE(construct_object_2);

  Construct_Apollonius_vertex_2 apollonius_vertex =
    tr.construct_Apollonius_vertex_2_object();
  CGAL_USE(apollonius_vertex);

  Construct_Apollonius_site_2 apollonius_site =
    tr.construct_Apollonius_site_2_object();
  CGAL_USE(apollonius_site);

  // testing correctness of predicates;
  //--------------------------------------------------------------------
  bool b;
  Point_2 p1__icc(10.0,-10);
  Site_2
    wp1(p1__icc,20),
    wp2(Point_2(9,-9),19.),
    wp3(Point_2(1000000,-1000000),0),
    wp4(Point_2(10.0,5.0),5),
    wp5;

  // testing coordinate comparison
  //--------------------------------------------------------------------
  assert( compare_x(wp1, wp1) == EQUAL );
  assert( compare_y(wp1, wp1) == EQUAL );
  assert( compare_w(wp1, wp1) == EQUAL );

  assert( compare_x(wp1, wp2) == LARGER );
  assert( compare_y(wp1, wp2) == SMALLER );
  assert( compare_w(wp1, wp2) == LARGER );

  // testing orientation
  //--------------------------------------------------------------------
  assert( orientation(wp1, wp2, wp3) == COLLINEAR );

  // testing is_hidden
  //--------------------------------------------------------------------
  assert( is_hidden(wp1, wp4) == true );
  assert( is_hidden(wp4, wp1) == false );

  // testing oriented_side_of_bisector
  //--------------------------------------------------------------------
  Point_2 p = wp2.point();
  assert( oriented_side_of_bisector(wp1, wp3, p) ==
                  ON_POSITIVE_SIDE );


  // testing vertex_conflict
  //--------------------------------------------------------------------

  // first we consider the case where all vertices are finite
  wp1 = Site_2(Point_2(0,0),0);
  wp2 = Site_2(Point_2(1,1),1);
  wp3 = Site_2(Point_2(2,4),4);
  wp4 = Site_2(Point_2(3,9),9);

  assert( vertex_conflict(wp1, wp2, wp3, wp4) == ZERO );
  assert( vertex_conflict(wp2, wp3, wp4, wp1) == ZERO );

  wp4 = Site_2(Point_2(3,9),9.0001);

  assert( vertex_conflict(wp1, wp2, wp3, wp4) == NEGATIVE );
  assert( vertex_conflict(wp2, wp3, wp4, wp1) == POSITIVE );


  // then we consider the case where v3 is the vertex at infinity
  wp1 = Site_2(Point_2(1,100),49);
  wp2 = Site_2(Point_2(0,-100),50);

  wp3 = Site_2(Point_2(0,0),0);
  assert( vertex_conflict(wp1, wp2, wp3) == POSITIVE );

  wp3 = Site_2(Point_2(-100000,0),0);
  assert( vertex_conflict(wp1, wp2, wp3) == POSITIVE );

  wp3 = Site_2(Point_2(100000,0),1000);
  assert( vertex_conflict(wp1, wp2, wp3) == NEGATIVE );

  wp3 = Site_2(Point_2(-1,0),51);
  assert( vertex_conflict(wp1, wp2, wp3) == NEGATIVE );

  wp3 = Site_2(Point_2(-1,200),51);
  assert( vertex_conflict(wp1, wp2, wp3) == POSITIVE );

  wp3 = Site_2(Point_2(-1,-200),51);
  assert( vertex_conflict(wp1, wp2, wp3) == POSITIVE );


  // testing finite_edge_interior_conflict
  //--------------------------------------------------------------------

  wp1 = Site_2(Point_2(0,-100),50);
  wp2 = Site_2(Point_2(0,100),51);
  wp3 = Site_2(Point_2(-150,0),49);
  wp4 = Site_2(Point_2(150,0),48);


  // first we look at the case where all vertices of the edge are
  // finite...

  // endpoints are not in conflict but the interior is
  wp5 = Site_2(Point_2(0,-49),0);

  assert( vertex_conflict(wp1, wp2, wp3, wp5) == POSITIVE );
  assert( vertex_conflict(wp1, wp4, wp2, wp5) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, wp5,  false);

  assert( b );

  // both endpoints and interior are in conflict
  wp5 = Site_2(Point_2(0,-49),10);

  assert( vertex_conflict(wp1, wp2, wp3, wp5) == NEGATIVE );
  assert( vertex_conflict(wp1, wp4, wp2, wp5) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, wp5,  true);

  assert( b );

  // endpoints are in conflict but the interior isn't
  wp5 = Site_2(Point_2(0,-150),99);

  assert( vertex_conflict(wp1, wp2, wp3, wp5) == NEGATIVE );
  assert( vertex_conflict(wp1, wp4, wp2, wp5) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, wp5,  true);

  assert( !b );

  // finally endpoints are not in conflict and interior is not in
  // conflict
  wp5 = Site_2(Point_2(0,-150),70);

  assert( vertex_conflict(wp1, wp2, wp3, wp5) == POSITIVE );
  assert( vertex_conflict(wp1, wp4, wp2, wp5) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, wp5,  true);

  assert( !b );


  // then we consider the case where v4 is the vertex at infinity
  wp1 = Site_2(Point_2(0,-100),50);
  wp2 = Site_2(Point_2(1,100),49);
  wp3 = Site_2(Point_2(-150,0),51);

  // endpoints are not in conflict but interior is
  wp4 = Site_2(Point_2(0,-48),0);
  assert( vertex_conflict(wp1, wp2, wp3, wp4) == POSITIVE );
  assert( vertex_conflict(wp2, wp1, wp4) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, false);
  assert( b );

  // both endpoints and interior are in conflict
  wp4 = Site_2(Point_2(-1,0),51);
  assert( vertex_conflict(wp1, wp2, wp3, wp4) == NEGATIVE );
  assert( vertex_conflict(wp2, wp1, wp4) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, true);
  assert( b );

  // endpoints are in conflict but interior isn't
  wp4 = Site_2(Point_2(-1,200),149);
  assert( vertex_conflict(wp1, wp2, wp3, wp4) == NEGATIVE );
  assert( vertex_conflict(wp2, wp1, wp4) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, true);
  assert( !b );

  // neither the endpoints nor the interior are in conflict
  wp4 = Site_2(Point_2(-1,200),51);

  assert( vertex_conflict(wp1, wp2, wp3, wp4) == POSITIVE );
  assert( vertex_conflict(wp2, wp1, wp4) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, wp4, false);
  assert( !b );


  // finally consider the case where both v3 and v4 are the vertex at
  // infinity.

  wp1 = Site_2(Point_2(10,0),5);
  wp2 = Site_2(Point_2(100,0),50);

  // endpoints are not in conflict but the interior is
  wp3 = Site_2(Point_2(20,0),0);

  assert( vertex_conflict(wp1, wp2, wp3) == POSITIVE );
  assert( vertex_conflict(wp2, wp1, wp3) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, false);
  assert( b );

  // endpoints are in conflict and so is the interior
  wp3 = Site_2(Point_2(20,0),10);

  assert( vertex_conflict(wp1, wp2, wp3) == NEGATIVE );
  assert( vertex_conflict(wp2, wp1, wp3) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, true);
  assert( b );

  // both the endpoints and the interior are not in conflict
  wp3 = Site_2(Point_2(4,0),2);

  assert( vertex_conflict(wp1, wp2, wp3) == POSITIVE );
  assert( vertex_conflict(wp2, wp1, wp3) == POSITIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, false);
  assert( !b );

  // the endpoints are in conflict but the interior isn't
  wp3 = Site_2(Point_2(4,0),3);

  assert( vertex_conflict(wp1, wp2, wp3) == NEGATIVE );
  assert( vertex_conflict(wp2, wp1, wp3) == NEGATIVE );

  b = finite_edge_interior_conflict(wp1, wp2, wp3, true);
  assert( !b );

  // testing infinite_edge_interior_conflict
  //--------------------------------------------------------------------


  wp2 = Site_2(Point_2(0,0),100);
  wp3 = Site_2(Point_2(-100,200),5);
  wp4 = Site_2(Point_2(100,300),4);

  // the endpoints are not in conflict but the interior is
  wp5 = Site_2(Point_2(0,-150),2);

  assert( vertex_conflict(wp2, wp3, wp5) == POSITIVE );
  assert( vertex_conflict(wp4, wp2, wp5) == POSITIVE );

  b = infinite_edge_interior_conflict(wp2, wp3, wp4, wp5, false);
  assert( b );


  // the endpoints are in conflict but the interior isn't
  wp5 = Site_2(Point_2(0,150),150);

  assert( vertex_conflict(wp2, wp3, wp5) == NEGATIVE );
  assert( vertex_conflict(wp4, wp2, wp5) == NEGATIVE );

  b = infinite_edge_interior_conflict(wp2, wp3, wp4, wp5, true);
  assert( !b );


  // the endpoints are in conflict as well as the interior
  wp5 = Site_2(Point_2(0,-150),150);

  assert( vertex_conflict(wp2, wp3, wp5) == NEGATIVE );
  assert( vertex_conflict(wp4, wp2, wp5) == NEGATIVE );

  b = infinite_edge_interior_conflict(wp2, wp3, wp4, wp5, true);
  assert( b );


  // neither the endpoints nor the interior are in conflict
  wp5 = Site_2(Point_2(0,150),50);

  assert( vertex_conflict(wp2, wp3, wp5) == POSITIVE );
  assert( vertex_conflict(wp4, wp2, wp5) == POSITIVE );

  b = infinite_edge_interior_conflict(wp2, wp3, wp4, wp5, false);
  assert( !b );



  // testing is_degenerate_edge
  //--------------------------------------------------------------------
  wp1 = Site_2(Point_2(2,4),4);
  wp2 = Site_2(Point_2(10,100),100);
  wp3 = Site_2(Point_2(5,25),25);
  wp4 = Site_2(Point_2(20,400),400);

  b = is_degenerate_edge(wp1, wp2, wp3, wp4);

  assert( b );

  wp1 = Site_2(Point_2(2,4),4);
  wp2 = Site_2(Point_2(10,100),100);
  wp3 = Site_2(Point_2(5,25),25);
  wp4 = Site_2(Point_2(20,400),399);

  b = is_degenerate_edge(wp1, wp2, wp3, wp4);

  assert( !b );

  return true;
}


#if defined(__INTEL_COMPILER)
template<class AG, class InputStream>
bool test_algo_generic(InputStream& is, const AG&)
#else
template<class AG, class InputStream>
bool test_algo_generic(InputStream& is)
#endif
{
  typedef AG                                      Apollonius_graph;
  typedef typename Apollonius_graph::Geom_traits  Traits;

  // testing typedefs
  //--------------------------------------------------------------------
  typedef typename Apollonius_graph::Data_structure    Data_structure;
  typedef typename Apollonius_graph::Geom_traits       Geom_traits;
  typedef typename Apollonius_graph::Point_2           Point_2;
  typedef typename Apollonius_graph::Site_2            Site_2;

  typedef typename Apollonius_graph::Edge              Edge;
  typedef typename Apollonius_graph::Vertex_handle     Vertex_handle;
  typedef typename Apollonius_graph::Face_handle       Face_handle;

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

  typedef typename Apollonius_graph::Sites_iterator Sites_iterator;
  typedef typename Apollonius_graph::Visible_sites_iterator
    Visible_sites_iterator;
  typedef typename Apollonius_graph::Hidden_sites_iterator
    Hidden_sites_iterator;

  typedef typename Apollonius_graph::size_type      size_type;

  CGAL_USE_TYPE(Data_structure);
  CGAL_USE_TYPE(Point_2);
  CGAL_USE_TYPE(Edge);

  Null_output_stream   nos;

  // testing creation/constructors
  //--------------------------------------------------------------------

  Apollonius_graph ag;
  Traits gt = Traits();
  Apollonius_graph ag1(gt);
  Apollonius_graph ag2(ag);

  std::vector<Site_2> wp_list;
  Site_2 wp;
  while ( is >> wp ) {
    wp_list.push_back(wp);
  }

  Apollonius_graph ag3(wp_list.begin(), wp_list.end(), Traits());

  ag = ag3;
  ag2 = ag3;

  // testing access functions
  //--------------------------------------------------------------------
  Geom_traits tr = ag.geom_traits();
  size_type num_vertices = ag.number_of_vertices();
  size_type num_all = num_vertices + ag.number_of_hidden_sites();
  // passing this to a dummy function to avoid warning when
  // CGAL_NO_ASSERTIONS is defined.
  dummy( static_cast<unsigned int>(num_all) == wp_list.size() );

  CGAL_USE(tr);

  Face_handle inf_f = ag.infinite_face();
  Vertex_handle v1 = ag.infinite_vertex();
  Vertex_handle v2 = ag.finite_vertex();

  // testing traversal - iterators
  //--------------------------------------------------------------------

  // finite faces, edges and vertices
  Finite_vertices_iterator fvit = ag.finite_vertices_begin();
  size_type n_fvertices = 0;
  for (; fvit != ag.finite_vertices_end(); ++fvit) {
    n_fvertices++;
  }
  assert( n_fvertices == num_vertices );

  Finite_edges_iterator feit = ag.finite_edges_begin();
  size_type n_fedges = 0;
  for (; feit != ag.finite_edges_end(); ++feit) {
    n_fedges++;
  }

  Finite_faces_iterator ffit = ag.finite_faces_begin();
  size_type n_ffaces = 0;
  for (; ffit != ag.finite_faces_end(); ++ffit) {
    n_ffaces++;
  }

  // all faces, edges and vertices
  All_vertices_iterator avit = ag.all_vertices_begin();
  size_type n_avertices = 0;
  for (; avit != ag.all_vertices_end(); ++avit) {
    n_avertices++;
  }
  assert( n_avertices == num_vertices + 1 );

  All_edges_iterator aeit = ag.all_edges_begin();
  size_type n_aedges = 0;
  for (; aeit != ag.all_edges_end(); ++aeit) {
    n_aedges++;
  }

  All_faces_iterator afit = ag.all_faces_begin();
  size_type n_afaces = 0;
  for (; afit != ag.all_faces_end(); ++afit) {
    n_afaces++;
  }


  assert( 2 * n_aedges == 3 * n_afaces );
  assert( n_avertices - n_aedges + n_afaces == 2 );

  // site iterators
  size_type n_sites(0), n_hidden_sites(0), n_visible_sites(0);

  for (Sites_iterator sit = ag.sites_begin();
       sit != ag.sites_end(); sit++) {
    n_sites++;
    nos << *sit;
    nos << sit->point();
  }

  for (Hidden_sites_iterator sit = ag.hidden_sites_begin();
       sit != ag.hidden_sites_end(); sit++) {
    n_hidden_sites++;
    nos << *sit;
    nos << sit->point();
  }

  for (Visible_sites_iterator sit = ag.visible_sites_begin();
       sit != ag.visible_sites_end(); sit++) {
    n_visible_sites++;
    nos << *sit;
    nos << sit->point();
  }

  assert( n_sites == n_visible_sites + n_hidden_sites );

  // testing traversal - circulators
  //--------------------------------------------------------------------

  // vertex circulators
  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Vertex_circulator vc_start = ag.incident_vertices(v);
    Vertex_circulator vc = vc_start;
    size_type deg = 0;
    do {
      deg++;
      vc++;
    } while ( vc != vc_start );

    assert( deg == ag.data_structure().degree(v) );
  }

  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Vertex_circulator vc_start = ag.incident_vertices(v, v->face());
    Vertex_circulator vc = vc_start;
    size_type deg = 0;
    do {
      deg++;
      vc++;
    } while ( vc != vc_start );

    assert( deg == ag.data_structure().degree(v) );
  }

  // face circulators
  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Face_circulator fc_start = ag.incident_faces(v);
    Face_circulator fc = fc_start;
    size_type deg = 0;
    do {
      deg++;
      fc++;
    } while ( fc != fc_start );

    assert( deg == ag.data_structure().degree(v) );
  }

  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Face_circulator fc_start = ag.incident_faces(v, v->face());
    Face_circulator fc = fc_start;
    size_type deg = 0;
    do {
      deg++;
      fc++;
    } while ( fc != fc_start );

    assert( deg == ag.data_structure().degree(v) );
  }

  // edge circulators
  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Edge_circulator ec_start = ag.incident_edges(v);
    Edge_circulator ec = ec_start;
    size_type deg = 0;
    do {
      deg++;
      ec++;
    } while ( ec != ec_start );

    assert( deg == ag.data_structure().degree(v) );
  }

  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end();
       ++avit) {
    Vertex_handle v(avit);
    Edge_circulator ec_start = ag.incident_edges(v, v->face());
    Edge_circulator ec = ec_start;
    size_type deg = 0;
    do {
      deg++;
      ec++;
    } while ( ec != ec_start );

    assert( deg == ag.data_structure().degree(v) );
  }


  // testing predicates
  //--------------------------------------------------------------------
  assert( ag.is_infinite(inf_f) );
  assert( ag.is_infinite(v1) );
  assert( !ag.is_infinite(v2) );
  {
    Edge_circulator ec = ag.incident_edges(ag.infinite_vertex());
    assert( ag.is_infinite(*ec) );
    assert( ag.is_infinite(ec) );
    assert( ag.is_infinite(ec->first, ec->second) );
  }

  // testing insertion
  //--------------------------------------------------------------------
  ag.clear();
  ag.insert(wp_list.begin(), wp_list.end());
  assert( ag.is_valid() );

  ag.clear();
  typename std::vector<Site_2>::iterator it;
  for (it = wp_list.begin(); it != wp_list.end(); ++it) {
    ag.insert(*it);
  }
  assert( ag.is_valid() );

  {
    ag.clear();
    Vertex_handle v;
    typename std::vector<Site_2>::iterator it;
    for (it = wp_list.begin(); it != wp_list.end(); ++it) {
      if ( it == wp_list.begin() ) {
        v = ag.insert(*it);
      } else {
        v = ag.insert(*it, v);
      }
    }
    assert( ag.is_valid() );
  }

  // testing removal
  //--------------------------------------------------------------------
  ag.clear();
  ag.insert(wp_list.begin(), wp_list.end());
  assert( ag.is_valid() );
  {
    while ( ag.number_of_vertices() > 0 ) {
      Vertex_handle v(ag.finite_vertices_begin());
      ag.remove(v);
      assert( ag.is_valid() );
    }
  }

  // testing nearest neighbor location
  //--------------------------------------------------------------------
  ag.clear();
  ag.insert(wp_list.begin(), wp_list.end());
  assert( ag.is_valid() );
  for (fvit = ag.finite_vertices_begin();
       fvit != ag.finite_vertices_end(); ++fvit) {
    Site_2 wp = fvit->site();
    Vertex_handle nn = ag.nearest_neighbor(wp.point());
    assert( wp == nn->site() );
  }

  // testing swap
  //--------------------------------------------------------------------
  ag.clear();
  ag.swap(ag2);
  assert( ag.number_of_vertices() +
                  ag.number_of_hidden_sites() == wp_list.size() );
  assert( ag2.number_of_vertices() == 0 );


  // drawing methods
  //--------------------------------------------------------------------
  ag.draw_primal(nos);
  ag.draw_dual(nos);

#if 0
  for (fvit = ag.finite_vertices_begin();
       fvit != ag.finite_vertices_end(); ++fvit) {
    ag.draw_primal_vertex(fvit, nos);
  }

  for (ffit = ag.finite_faces_begin();
       ffit != ag.finite_faces_end(); ++ffit) {
    ag.draw_dual_vertex(ffit, nos);
  }
#endif

  for (feit = ag.finite_edges_begin();
       feit != ag.finite_edges_end(); ++feit) {
    ag.draw_primal_edge(feit, nos);
    ag.draw_dual_edge(feit, nos);
  }

  for (aeit = ag.all_edges_begin();
       aeit != ag.all_edges_end(); ++aeit) {
    ag.draw_primal_edge(aeit, nos);
  }

#if 0
  for (afit = ag.all_faces_begin();
       afit != ag.all_faces_end(); ++afit) {
    ag.draw_primal_face(afit, nos);
  }

  for (avit = ag.all_vertices_begin();
       avit != ag.all_vertices_end(); ++avit) {
    ag.draw_dual_face(avit, nos);
  }
#endif

  // file I/O methods
  //--------------------------------------------------------------------
  {
    std::string fname  =  "ag_testsuite_" + std::to_string(CGAL::Random().get_seed())  + ".tmp";
    std::cout << "writing to " << fname << std::endl;

    std::ofstream ofs(fname);
    assert( ofs );
    ag.file_output(ofs);
    ofs.close();

    std::ifstream ifs(fname);
    assert( ifs );
    ag.file_input(ifs);
    ifs.close();
    assert( ag.is_valid() );
  }
  {
    std::string fname = "ag_testsuite_" + std::to_string(CGAL::Random().get_seed()) + ".tmp";
    std::cout << "writing to " << fname << std::endl;
    std::ofstream ofs(fname);
    assert( ofs );
    ofs << ag;
    ofs.close();

    std::ifstream ifs(fname);
    assert( ifs );
    ifs >> ag;
    ifs.close();
    assert( ag.is_valid() );
  }
  return true;
}

#if defined(__INTEL_COMPILER)
template<class Kernel, class Method_tag, class InputStream>
bool test_algo(InputStream& is, const Kernel&,
               const Method_tag&)
#else
template<class Kernel, class Method_tag, class InputStream>
bool test_algo(InputStream& is)
#endif
{
  typedef Apollonius_graph_traits_2<Kernel,Method_tag> Traits;
#if defined( _MSC_VER )
  // Patch for the Microsoft compiler so that it does not produce the
  // nasty warning about decorated name length
  // Basically what I do here is create typedefs for the default
  // template parameters so as to give them shorter names
  typedef Apollonius_graph_vertex_base_2<Traits,true>   Vb;
  typedef Triangulation_face_base_2<Traits>             Fb;
  typedef Triangulation_data_structure_2<Vb,Fb>         Agds;
  typedef Apollonius_graph_2<Traits,Agds>               Apollonius_graph;
#else
  typedef Apollonius_graph_2<Traits>  Apollonius_graph;
#endif

  return test_algo_generic<Apollonius_graph,InputStream>(is);
}

#if defined(__INTEL_COMPILER)
template<class Kernel, class Method_tag, class InputStream>
bool test_hierarchy_algo(InputStream& is, const Kernel&, const Method_tag&)
#else
template<class Kernel, class Method_tag, class InputStream>
bool test_hierarchy_algo(InputStream& is)
#endif
{
  typedef Apollonius_graph_traits_2<Kernel,Method_tag> Traits;
#if defined( _MSC_VER )
  // Patch for the Microsoft compiler so that it does not produce the
  // nasty warning about decorated name length
  // Basically what I do here is create typedefs for the default
  // template parameters so as to give them shorter names
  typedef Apollonius_graph_vertex_base_2<Traits,true>   Vb;
  typedef Apollonius_graph_hierarchy_vertex_base_2<Vb>  HVb;
  typedef Triangulation_face_base_2<Traits>             Fb;
  typedef Triangulation_data_structure_2<HVb,Fb>        Agds;
  typedef Apollonius_graph_hierarchy_2<Traits,Agds>
    Apollonius_graph_hierarchy;
#else
  typedef Apollonius_graph_hierarchy_2<Traits>  Apollonius_graph_hierarchy;
#endif

  return test_algo_generic<Apollonius_graph_hierarchy,InputStream>(is);
}



#if defined(__INTEL_COMPILER)
template<class CK, class CKM, class EK, class EKM, class InputStream>
bool test_filtered_traits_algo(InputStream& is, const CK&,
                               const CKM&, const EK&, const EKM&)
#else
template<class CK, class CKM, class EK, class EKM, class InputStream>
bool test_filtered_traits_algo(InputStream& is)
#endif
{
  typedef Apollonius_graph_filtered_traits_2<CK,CKM,EK,EKM> Traits;
#if defined( _MSC_VER )
  // Patch for the Microsoft compiler so that it does not produce the
  // nasty warning about decorated name length
  // Basically what I do here is create typedefs for the default
  // template parameters so as to give them shorter names
  typedef Apollonius_graph_vertex_base_2<Traits,true>   Vb;
  typedef Triangulation_face_base_2<Traits>             Fb;
  typedef Triangulation_data_structure_2<Vb,Fb>         Agds;
  typedef Apollonius_graph_2<Traits,Agds>               Apollonius_graph;
#else
  typedef Apollonius_graph_2<Traits>  Apollonius_graph;
#endif

  return test_algo_generic<Apollonius_graph,InputStream>(is);
}

#if defined(__INTEL_COMPILER)
template<class CK, class CKM, class EK, class EKM, class InputStream>
bool test_filtered_traits_hierarchy_algo(InputStream& is,
                                         const CK&, const CKM&,
                                         const EK&, const EKM&)
#else
template<class CK, class CKM, class EK, class EKM, class InputStream>
bool test_filtered_traits_hierarchy_algo(InputStream& is)
#endif
{
  typedef Apollonius_graph_filtered_traits_2<CK,CKM,EK,EKM> Traits;
#if defined( _MSC_VER )
  // Patch for the Microsoft compiler so that it does not produce the
  // nasty warning about decorated name length
  // Basically what I do here is create typedefs for the default
  // template parameters so as to give them shorter names
  typedef Apollonius_graph_vertex_base_2<Traits,true>   Vb;
  typedef Apollonius_graph_hierarchy_vertex_base_2<Vb>  HVb;
  typedef Triangulation_face_base_2<Traits>             Fb;
  typedef Triangulation_data_structure_2<HVb,Fb>        Agds;
  typedef Apollonius_graph_hierarchy_2<Traits,Agds>
    Apollonius_graph_hierarchy;
#else
  typedef Apollonius_graph_hierarchy_2<Traits>  Apollonius_graph_hierarchy;
#endif

  return test_algo_generic<Apollonius_graph_hierarchy,InputStream>(is);
}



} //namespace CGAL



#endif // CGAL_APOLLONIUS_GRAPH_2_TEST_H

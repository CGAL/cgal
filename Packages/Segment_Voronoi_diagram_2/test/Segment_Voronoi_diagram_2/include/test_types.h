#ifndef CGAL_SVD_TEST_TYPES_H
#define CGAL_SVD_TEST_TYPES_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <iostream>
#include <algorithm>

#include "IO/Null_output_stream.h"
#include "IO/io_aux.h"

CGAL_BEGIN_NAMESPACE

template<class SVD, class InputStream>
bool test_svd(InputStream& is, const SVD&)
{
  typedef SVD SVD2;

  typedef SVD2 Segment_Voronoi_diagram_2;

  start_testing("type definitions");
  typedef typename SVD2::Geom_traits               Geom_traits;
  typedef typename SVD2::Data_structure            Data_structure;
  typedef typename SVD2::size_type                 size_type;
  typedef typename SVD2::Point_2                   Point_2;
  typedef typename SVD2::Site_2                    Site_2;
  typedef typename SVD2::Point_container           Point_container;
  typedef typename SVD2::Point_handle              Point_handle;

  typedef typename SVD2::Edge                      Edge;
  typedef typename SVD2::Vertex_handle             Vertex_handle;
  typedef typename SVD2::Face_handle               Face_handle;

  typedef typename SVD2::Vertex_circulator         Vertex_circulator;
  typedef typename SVD2::Face_circulator           Face_circulator;
  typedef typename SVD2::Edge_circulator           Edge_circulator;

  typedef typename SVD2::All_vertices_iterator     All_vertices_iterator;
  typedef typename SVD2::Finite_vertices_iterator  Finite_vertices_iterator;

  typedef typename SVD2::All_faces_iterator        All_faces_iterator;
  typedef typename SVD2::Finite_faces_iterator     Finite_faces_iterator;

  typedef typename SVD2::All_edges_iterator        All_edges_iterator;
  typedef typename SVD2::Finite_edges_iterator     Finite_edges_iterator;

  typedef typename SVD2::Input_sites_iterator      Input_sites_iterator;
  typedef typename SVD2::Output_sites_iterator     Output_sites_iterator;
  end_testing("type definitions");

  Point_2 p1(0,0), p2(0,1), p3(1,0), p4(1,1);
  Site_2 t1 = Site_2::construct_site_2(p1);
  Site_2 t2 = Site_2::construct_site_2(p2);
  Site_2 t3 = Site_2::construct_site_2(p3);
  Site_2 t4 = Site_2::construct_site_2(p4);
  Site_2 s1 = Site_2::construct_site_2(p1,p2);
  Site_2 s2 = Site_2::construct_site_2(p2,p3);
  Site_2 s3 = Site_2::construct_site_2(p3,p4);
  Site_2 s4 = Site_2::construct_site_2(p4,p1);

  std::vector<Site_2> site_list;
  {
    site_list.push_back(t1);
    site_list.push_back(t2);
    site_list.push_back(t3);
    site_list.push_back(t4);
    site_list.push_back(s1);
    site_list.push_back(s2);
    site_list.push_back(s3);
    site_list.push_back(s4);
  }

  Null_output_stream null_os;

  Geom_traits gt;

  start_testing("constructors");
  Segment_Voronoi_diagram_2 svd;
  CGAL_assertion( svd.is_valid() );

  Segment_Voronoi_diagram_2 svd2(gt);
  CGAL_assertion( svd2.is_valid() );
  {
    Segment_Voronoi_diagram_2 svd3(svd);

    CGAL_assertion( svd3.is_valid() );

    svd.insert(site_list.begin(), site_list.end());
    Segment_Voronoi_diagram_2 svd4(svd);
    CGAL_assertion( svd4.is_valid() );
    svd.clear();

    Segment_Voronoi_diagram_2 svd5(site_list.begin(), site_list.end());
    CGAL_assertion( svd5.is_valid() );
    Segment_Voronoi_diagram_2 svd6(site_list.begin(), site_list.end(), gt);
    CGAL_assertion( svd6.is_valid() );
  }
  end_testing("constructors");

  start_testing("assignment operator");
  svd.insert(site_list.begin(), site_list.end());

  svd = svd;
  svd2 = svd;

  CGAL_assertion( svd.is_valid() );
  CGAL_assertion( svd2.is_valid() );

  svd.clear();

  CGAL_assertion( svd.is_valid() );
  CGAL_assertion( svd2.is_valid() );

  svd = svd2;

  end_testing("assignment operator");


  start_testing("access methods");
  gt = svd.geom_traits();
  int dim = svd.dimension();
  null_os << dim;

  size_type nv = svd.number_of_vertices();
  null_os << nv;

  size_type nf = svd.number_of_faces();
  null_os << nf;

  size_type nis = svd.number_of_input_sites();
  null_os << nis;

  size_type nos = svd.number_of_output_sites();
  null_os << nos;

  Face_handle inf_f = svd.infinite_face();
  null_os << inf_f;
  Vertex_handle inf_v = svd.infinite_vertex();
  null_os << inf_v;

  if ( nv > 0 ) {
    Vertex_handle fin_v = svd.finite_vertex();
    null_os << fin_v;
  }
  

  Data_structure ds = svd.data_structure();
  Point_container pc = svd.point_container();
  end_testing("access methods");

  start_testing("iterators and circulators");
  {
    Finite_vertices_iterator fit = svd.finite_vertices_begin();
    for (; fit != svd.finite_vertices_end(); ++fit) {
      Vertex_handle v(fit);
      null_os << v << fit;
    }

    All_vertices_iterator ait = svd.all_vertices_begin();
    for (; ait != svd.all_vertices_end(); ++ait) {
      Vertex_handle v(ait);
      null_os << v << ait;      
    }
  }
  {
    Finite_faces_iterator fit = svd.finite_faces_begin();
    for (; fit != svd.finite_faces_end(); ++fit) {
      Face_handle f(fit);
      null_os << f << fit;
    }

    All_faces_iterator ait = svd.all_faces_begin();
    for (; ait != svd.all_faces_end(); ++ait) {
      Face_handle f(ait);
      null_os << f << ait;      
    }
  }
  {
    Finite_edges_iterator fit = svd.finite_edges_begin();
    for (; fit != svd.finite_edges_end(); ++fit) {
      Edge e = *fit;
      null_os << e << fit;
    }

    All_edges_iterator ait = svd.all_edges_begin();
    for (; ait != svd.all_edges_end(); ++ait) {
      Edge e = *ait;
      null_os << e << ait;      
    }
  }
  {
    if ( nv > 0 ) {
      Vertex_circulator vc = svd.incident_vertices(svd.infinite_vertex());
      Vertex_circulator vc_start = vc;
      do {
	null_os << vc;
	vc++;
      } while ( vc != vc_start );

      Face_circulator fc = svd.incident_faces(svd.infinite_vertex());
      Face_circulator fc_start = fc;
      do {
	null_os << fc;
	fc++;
      } while ( fc != fc_start );

      Edge_circulator ec = svd.incident_edges(svd.infinite_vertex());
      Edge_circulator ec_start = ec;
      do {
	null_os << ec;
	ec++;
      } while ( ec != ec_start );
    }
    
  }
  {
    Input_sites_iterator isi = svd.input_sites_begin();
    for (; isi != svd.input_sites_end(); ++isi) {
      null_os << *isi;
    }

    Output_sites_iterator osi = svd.output_sites_begin();
    for (; osi != svd.output_sites_end(); ++osi) {
      null_os << *osi;
    }
  }
  end_testing("iterators and circulators");

  start_testing("predicates");
  {
    bool is_inf;
    is_inf = svd.is_infinite(svd.infinite_vertex());
    null_os << is_inf;

    if ( nv > 0 ) {
      is_inf = svd.is_infinite(svd.finite_vertex());
      null_os << is_inf;

      Face_handle f( svd.finite_faces_begin() );
      is_inf = svd.is_infinite( f );
      null_os << is_inf;

      is_inf = svd.is_infinite( f, 0 );
      null_os << is_inf;      

      Edge e = *svd.finite_edges_begin();
      is_inf = svd.is_infinite( e );
      null_os << is_inf;      

      Edge_circulator ec = svd.incident_edges(svd.infinite_vertex());
      is_inf = svd.is_infinite( ec );
      null_os << is_inf;    
    }
  }
  end_testing("predicates");

  start_testing("insertion methods");
  {
    svd.clear();
    svd.insert(site_list.begin(), site_list.end());
    CGAL_assertion( svd.is_valid() );

    svd.clear();
    svd.insert(site_list.begin(), site_list.end(), Tag_false());
    CGAL_assertion( svd.is_valid() );

    svd.clear();
    svd.insert(site_list.begin(), site_list.end(), Tag_true());
    CGAL_assertion( svd.is_valid() );

    svd.clear();
    Vertex_handle v1 = svd.insert(p1);
    svd.insert(p2, v1);
    svd.insert(p1, p2);
    svd.insert(p2, p3, v1);
    svd.insert(t3);
    svd.insert(t4, v1);

    svd.insert(s3, v1);
    svd.insert(s4, v1);

    CGAL_assertion( svd.is_valid() );
  }
  end_testing("insertion methods");

  start_testing("nearest neighbor methods");
  {
    Vertex_handle vnearest = svd.nearest_neighbor(p1);
    vnearest = svd.nearest_neighbor(p2, vnearest);
  }
  end_testing("nearest neighbor methods");

  start_testing("drawing methods");
  {
    svd.draw_dual(null_os);
    svd.draw_skeleton(null_os);
    svd.draw_dual_edge(*svd.finite_edges_begin(), null_os);
    svd.draw_dual_edge(svd.finite_edges_begin(), null_os);
  }
  end_testing("drawing methods");

  start_testing("swap method");
  {
    svd.swap(svd);

    size_type i1 = svd.number_of_input_sites();
    size_type n1 = svd.number_of_vertices();

    svd2.clear();
    svd.swap(svd2);

    CGAL_assertion( svd2.number_of_input_sites() == i1 );
    CGAL_assertion( svd2.number_of_vertices() == n1 );
    CGAL_assertion( svd2.is_valid() );

    CGAL_assertion( svd.number_of_input_sites() == 0 );
    CGAL_assertion( svd.number_of_vertices() == 0 );
    CGAL_assertion( svd.is_valid() );

    //--------------------------------------------
 
    svd2.swap(svd); // back to original setting (before swap)

    svd2.insert(p4);
    svd2.insert(p3);
    svd2.insert(p3, p4);

    i1 = svd.number_of_input_sites();
    n1 = svd.number_of_vertices();

    size_type i2 = svd2.number_of_input_sites();
    size_type n2 = svd2.number_of_vertices();

    svd.swap(svd2);

    CGAL_assertion( svd2.number_of_input_sites() == i1 );
    CGAL_assertion( svd2.number_of_vertices() == n1 );
    CGAL_assertion( svd2.is_valid() );

    CGAL_assertion( svd.number_of_input_sites() == i2 );
    CGAL_assertion( svd.number_of_vertices() == n2 );
    CGAL_assertion( svd.is_valid() );

    //--------------------------------------------

    i2 = svd2.number_of_input_sites();
    n2 = svd2.number_of_vertices();

    svd.clear();
    svd.swap(svd2);

    CGAL_assertion( svd2.number_of_input_sites() == 0 );
    CGAL_assertion( svd2.number_of_vertices() == 0 );
    CGAL_assertion( svd2.is_valid() );

    CGAL_assertion( svd.number_of_input_sites() == i2 );
    CGAL_assertion( svd.number_of_vertices() == n2 );
    CGAL_assertion( svd.is_valid() );
  }
  end_testing("swap method");

  start_testing("validity check and clear methods");
  std::cout << std::endl;
  {
    std::cout << "  validating: " << std::flush;
    svd.is_valid(true, 1);
    std::cout << std::endl;
    std::cout << "  clearing diagram..." << std::endl;
    svd.clear();
    std::cout << "  validating: " << std::flush;
    svd.is_valid(true, 1);
    std::cout << std::endl;
  }
  start_testing("validity check and clear methods");
  end_testing("validity check and clear method");

  return true;
}

#if 0


#if defined(__INTEL_COMPILER)
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
  int num_vertices = ag.number_of_vertices();
  int num_all = num_vertices + ag.number_of_hidden_sites();
  CGAL_assertion( static_cast<unsigned int>(num_all) == wp_list.size() );

  Face_handle inf_f = ag.infinite_face();
  Vertex_handle v1 = ag.infinite_vertex();
  Vertex_handle v2 = ag.finite_vertex();

  // testing traversal - iterators
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

  // site iterators
  int n_sites(0), n_hidden_sites(0), n_visible_sites(0);

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

  CGAL_assertion( n_sites == n_visible_sites + n_hidden_sites );

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

    CGAL_assertion( deg == v->degree() );
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

    CGAL_assertion( deg == v->degree() );
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

    CGAL_assertion( deg == v->degree() );
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

    CGAL_assertion( deg == v->degree() );
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

    CGAL_assertion( deg == v->degree() );
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
  typename std::vector<Site_2>::iterator it;
  for (it = wp_list.begin(); it != wp_list.end(); ++it) {
    ag.insert(*it);
  }
  CGAL_assertion( ag.is_valid() );

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
    Site_2 wp = fvit->site();
    Vertex_handle nn = ag.nearest_neighbor(wp.point());
    CGAL_assertion( wp == nn->site() );
  }

  // testing swap
  //--------------------------------------------------------------------
  ag.clear();
  ag.swap(ag2);
  CGAL_assertion( ag.number_of_vertices() +
		  ag.number_of_hidden_sites() == wp_list.size() );
  CGAL_assertion( ag2.number_of_vertices() == 0 );


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
  // template paramaters so as to give them shorter names
  typedef Apollonius_graph_vertex_base_2<Traits,true>   Vb;
  typedef Apollonius_graph_face_base_2<Traits>          Fb;
  typedef Apollonius_graph_data_structure_2<Vb,Fb>      Agds;
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
  // template paramaters so as to give them shorter names
  typedef Apollonius_graph_vertex_base_2<Traits,true>   Vb;
  typedef Apollonius_graph_hierarchy_vertex_base_2<Vb>  HVb;
  typedef Apollonius_graph_face_base_2<Traits>          Fb;
  typedef Apollonius_graph_data_structure_2<HVb,Fb>     Agds;
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
  // template paramaters so as to give them shorter names
  typedef Apollonius_graph_vertex_base_2<Traits,true>   Vb;
  typedef Apollonius_graph_face_base_2<Traits>          Fb;
  typedef Apollonius_graph_data_structure_2<Vb,Fb>      Agds;
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
  // template paramaters so as to give them shorter names
  typedef Apollonius_graph_vertex_base_2<Traits,true>   Vb;
  typedef Apollonius_graph_hierarchy_vertex_base_2<Vb>  HVb;
  typedef Apollonius_graph_face_base_2<Traits>          Fb;
  typedef Apollonius_graph_data_structure_2<HVb,Fb>     Agds;
  typedef Apollonius_graph_hierarchy_2<Traits,Agds>
    Apollonius_graph_hierarchy;
#else
  typedef Apollonius_graph_hierarchy_2<Traits>  Apollonius_graph_hierarchy;
#endif

  return test_algo_generic<Apollonius_graph_hierarchy,InputStream>(is);
}

#endif

CGAL_END_NAMESPACE



#endif // CGAL_SVD_TEST_TYPES_H

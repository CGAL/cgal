#ifndef CGAL_SDG_TEST_TYPES_H
#define CGAL_SDG_TEST_TYPES_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/use.h>
#include <iostream>
#include <cctype>
#include <algorithm>
#include <cassert>
#include <string>

#include "IO/Null_output_stream.h"
#include "IO/io_aux.h"

#include <cstring> // for std::strcpy, and std::strcat

//========================================================================

template<class NT>
std::string get_fname(const NT&, std::string ifname)
{
  return std::string("data/") + ifname + ".cin";
}

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpq.h>

template<>
std::string get_fname(const CGAL::Gmpq&, std::string ifname)
{
  return std::string("data/") + ifname + ".Gmpq.cin";
}
#endif

//========================================================================

namespace CGAL {

template<class SDG>
struct Level_finder
{
  typedef typename SDG::size_type      size_type;
  typedef typename SDG::Vertex_handle  Vertex_handle;

  size_type operator()(Vertex_handle v) const {
    assert( v != Vertex_handle() );
    size_type level = 0;
    Vertex_handle vertex = v;
    while ( vertex->up() != Vertex_handle() ) {
      vertex = vertex->up();
      level++;
    }

    return level;
  }
};

template<class Gt, class SDGDS, class LTag>
struct Level_finder< Segment_Delaunay_graph_Linf_2<Gt,SDGDS,LTag> >
{
  typedef Segment_Delaunay_graph_Linf_2<Gt,SDGDS,LTag> SDG;

  typedef typename SDG::size_type      size_type;
  typedef typename SDG::Vertex_handle  Vertex_handle;

  size_type operator()(Vertex_handle ) const { return 0; }
};


//========================================================================

template<class SDG, class InputStream>
bool test_sdg(InputStream&, const SDG&, const char* ifname, const char* ofname,
	      bool test_remove)
{
  std::string ifname_full = get_fname(typename SDG2::Geom_traits::FT(), ifname);

  typedef SDG SDG2;

  typedef SDG2 Segment_Delaunay_graph_2;

  start_testing("type definitions");
  typedef typename SDG2::Geom_traits               Geom_traits;
  typedef typename SDG2::Data_structure            Data_structure;
  typedef typename SDG2::size_type                 size_type;
  typedef typename SDG2::Point_2                   Point_2;
  typedef typename SDG2::Site_2                    Site_2;
  typedef typename SDG2::Point_container           Point_container;
  typedef typename SDG2::Point_handle              Point_handle;
  CGAL_USE_TYPE(Point_handle);

  typedef typename SDG2::Edge                      Edge;
  typedef typename SDG2::Vertex_handle             Vertex_handle;
  typedef typename SDG2::Face_handle               Face_handle;

  typedef typename SDG2::Vertex_circulator         Vertex_circulator;
  typedef typename SDG2::Face_circulator           Face_circulator;
  typedef typename SDG2::Edge_circulator           Edge_circulator;

  typedef typename SDG2::All_vertices_iterator     All_vertices_iterator;
  typedef typename SDG2::Finite_vertices_iterator  Finite_vertices_iterator;

  typedef typename SDG2::All_faces_iterator        All_faces_iterator;
  typedef typename SDG2::Finite_faces_iterator     Finite_faces_iterator;

  typedef typename SDG2::All_edges_iterator        All_edges_iterator;
  typedef typename SDG2::Finite_edges_iterator     Finite_edges_iterator;

  typedef typename SDG2::Input_sites_iterator      Input_sites_iterator;
  typedef typename SDG2::Output_sites_iterator     Output_sites_iterator;
  end_testing("type definitions");

  Point_2 p1(0,0), p2(0,1), p3(1,0), p4(1,1);
  //  Point_2 p5(10,1);
  Site_2 t1 = Site_2::construct_site_2(p1);
  Site_2 t2 = Site_2::construct_site_2(p2);
  Site_2 t3 = Site_2::construct_site_2(p3);
  Site_2 t4 = Site_2::construct_site_2(p4);
  Site_2 s1 = Site_2::construct_site_2(p1,p2);
  Site_2 s2 = Site_2::construct_site_2(p2,p3);
  Site_2 s3 = Site_2::construct_site_2(p3,p4);
  Site_2 s4 = Site_2::construct_site_2(p4,p1);

  //  Site_2 t5 = Site_2::construct_site_2(p5);

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
  Segment_Delaunay_graph_2 sdg;
  assert( sdg.is_valid() );

  Segment_Delaunay_graph_2 sdg2(gt);
  assert( sdg2.is_valid() );
  {
    Segment_Delaunay_graph_2 sdg3(sdg);

    assert( sdg3.is_valid() );

    sdg.insert(site_list.begin(), site_list.end());
    Segment_Delaunay_graph_2 sdg4(sdg);
    assert( sdg4.is_valid() );
    sdg.clear();

    Segment_Delaunay_graph_2 sdg5(site_list.begin(), site_list.end());
    assert( sdg5.is_valid() );
    Segment_Delaunay_graph_2 sdg6(site_list.begin(), site_list.end(), gt);
    assert( sdg6.is_valid() );
  }
  end_testing("constructors");

  start_testing("assignment operator");
  sdg.insert(site_list.begin(), site_list.end());

  sdg = sdg;
  sdg2 = sdg;

  assert( sdg.is_valid() );
  assert( sdg2.is_valid() );

  sdg.clear();

  assert( sdg.is_valid() );
  assert( sdg2.is_valid() );

  sdg = sdg2;

  end_testing("assignment operator");


  start_testing("access methods");
  gt = sdg.geom_traits();
  int dim = sdg.dimension();
  null_os << dim;

  size_type nv = sdg.number_of_vertices();
  null_os << nv;

  size_type nf = sdg.number_of_faces();
  null_os << nf;

  size_type nis = sdg.number_of_input_sites();
  null_os << nis;

  size_type nos = sdg.number_of_output_sites();
  null_os << nos;

  Face_handle inf_f = sdg.infinite_face();
  null_os << inf_f;
  Vertex_handle inf_v = sdg.infinite_vertex();
  null_os << inf_v;

  if ( nv > 0 ) {
    Vertex_handle fin_v = sdg.finite_vertex();
    null_os << fin_v;
  }


  Data_structure ds = sdg.data_structure();
  Point_container pc = sdg.point_container();
  end_testing("access methods");

  start_testing("iterators and circulators");
  {
    Finite_vertices_iterator fit = sdg.finite_vertices_begin();
    for (; fit != sdg.finite_vertices_end(); ++fit) {
      Vertex_handle v(fit);
      null_os << v << fit;
    }

    All_vertices_iterator ait = sdg.all_vertices_begin();
    for (; ait != sdg.all_vertices_end(); ++ait) {
      Vertex_handle v(ait);
      null_os << v << ait;
    }
  }
  {
    Finite_faces_iterator fit = sdg.finite_faces_begin();
    for (; fit != sdg.finite_faces_end(); ++fit) {
      Face_handle f(fit);
      null_os << f << fit;
    }

    All_faces_iterator ait = sdg.all_faces_begin();
    for (; ait != sdg.all_faces_end(); ++ait) {
      Face_handle f(ait);
      null_os << f << ait;
    }
  }
  {
    Finite_edges_iterator fit = sdg.finite_edges_begin();
    for (; fit != sdg.finite_edges_end(); ++fit) {
      Edge e = *fit;
      null_os << e << fit;
    }

    All_edges_iterator ait = sdg.all_edges_begin();
    for (; ait != sdg.all_edges_end(); ++ait) {
      Edge e = *ait;
      null_os << e << ait;
    }
  }
  {
    if ( nv > 0 ) {
      Vertex_circulator vc = sdg.incident_vertices(sdg.infinite_vertex());
      Vertex_circulator vc_start = vc;
      do {
	null_os << vc;
	vc++;
      } while ( vc != vc_start );

      Face_circulator fc = sdg.incident_faces(sdg.infinite_vertex());
      Face_circulator fc_start = fc;
      do {
	null_os << fc;
	fc++;
      } while ( fc != fc_start );

      Edge_circulator ec = sdg.incident_edges(sdg.infinite_vertex());
      Edge_circulator ec_start = ec;
      do {
	null_os << ec;
	ec++;
      } while ( ec != ec_start );
    }

  }
  //  sdg.insert(t5);
  //  sdg2.clear();
  //  std::cout << std::endl;
  {
    Input_sites_iterator isi = sdg.input_sites_begin();
    for (; isi != sdg.input_sites_end(); ++isi) {
      null_os << *isi;
      //      std::cout << *isi << std::endl;
    }
    //    std::cout << std::endl;

    Output_sites_iterator osi = sdg.output_sites_begin();
    for (; osi != sdg.output_sites_end(); ++osi) {
      null_os << *osi;
      //      std::cout << *osi << std::endl;
    }
    //    std::cout << std::endl;
  }
  end_testing("iterators and circulators");

  start_testing("predicates");
  {
    bool is_inf;
    is_inf = sdg.is_infinite(sdg.infinite_vertex());
    null_os << is_inf;

    if ( nv > 0 ) {
      is_inf = sdg.is_infinite(sdg.finite_vertex());
      null_os << is_inf;

      Face_handle f( sdg.finite_faces_begin() );
      is_inf = sdg.is_infinite( f );
      null_os << is_inf;

      is_inf = sdg.is_infinite( f, 0 );
      null_os << is_inf;

      Edge e = *sdg.finite_edges_begin();
      is_inf = sdg.is_infinite( e );
      null_os << is_inf;

      Edge_circulator ec = sdg.incident_edges(sdg.infinite_vertex());
      is_inf = sdg.is_infinite( ec );
      null_os << is_inf;
    }
  }
  end_testing("predicates");

  start_testing("insertion methods");
  {
    sdg.clear();
    sdg.insert(site_list.begin(), site_list.end());
    assert( sdg.is_valid() );

    sdg.clear();
    sdg.insert(site_list.begin(), site_list.end(), Tag_false());
    assert( sdg.is_valid() );

    sdg.clear();
    sdg.insert(site_list.begin(), site_list.end(), Tag_true());
    assert( sdg.is_valid() );

    sdg.clear();
    Vertex_handle v1 = sdg.insert(p1);
    sdg.insert(p2, v1);
    sdg.insert(p1, p2);
    sdg.insert(p2, p3, v1);
    sdg.insert(t3);
    sdg.insert(t4, v1);

    sdg.insert(s3, v1);
    sdg.insert(s4, v1);

    assert( sdg.is_valid() );
  }
  end_testing("insertion methods");

  if ( test_remove ) {
    start_testing("removal methods");
    {
      sdg.clear();

      std::ifstream ifs( ifname_full.c_str() );
      assert( ifs );
      Site_2 t;
      while ( ifs >> t ) {
	sdg.insert(t);
      }
      std::cerr << std::endl;
      assert( sdg.is_valid(true, 1) );

      Finite_vertices_iterator vit;

      std::vector<Vertex_handle> vec;
      for (vit = sdg.finite_vertices_begin();
	   vit != sdg.finite_vertices_end(); ++vit) {
	vec.push_back(vit);
      }
      std::random_shuffle(vec.begin(), vec.end());

      typename std::vector<Vertex_handle>::iterator it = vec.begin();
      std::cerr << std::endl;
      Level_finder<SDG> level;
      do {
	Site_2 tt = (*it)->site();
	std::cerr << "  *** attempting to remove: " << tt << "\t" << std::flush;
	if ( tt.is_point() ) { std::cerr << "\t" << std::flush; }
	std::cerr << "  - LEVEL of v: " << level(*it)
		  << " -  " << std::flush;
	bool success = sdg.remove(*it);
	std::cerr << (success ? " successful" : " UNSUCCESSFUL") << std::endl;
	if ( success ) {
	  vec.erase(it);
	  it = vec.begin();
	} else {
	  ++it;
	}
      } while ( it != vec.end() );

      // second test case
      sdg.clear();

      sdg.insert(site_list.begin(), site_list.end());
      assert( sdg.is_valid() );

      vit = sdg.finite_vertices_begin();
      do {
	bool success = sdg.remove(vit);
	if ( success ) {
	  vit = sdg.finite_vertices_begin();
	} else {
	  ++vit;
	}
      } while ( vit != sdg.finite_vertices_end() );
    }
    end_testing("removal methods");
  }

  {
    // recover state of sdg after testing removals
    Vertex_handle v1 = sdg.insert(p1);
    sdg.insert(p2, v1);
    sdg.insert(p1, p2);
    sdg.insert(p2, p3, v1);
    sdg.insert(t3);
    sdg.insert(t4, v1);

    sdg.insert(s3, v1);
    sdg.insert(s4, v1);

    assert( sdg.is_valid() );
  }

  start_testing("nearest neighbor methods");
  {
    Vertex_handle vnearest = sdg.nearest_neighbor(p1);
    vnearest = sdg.nearest_neighbor(p2, vnearest);
  }
  end_testing("nearest neighbor methods");

  start_testing("drawing methods");
  {
    sdg.draw_dual(null_os);
    sdg.draw_skeleton(null_os);
    sdg.draw_dual_edge(*sdg.finite_edges_begin(), null_os);
    sdg.draw_dual_edge(sdg.finite_edges_begin(), null_os);
  }
  end_testing("drawing methods");

  start_testing("swap method");
  {
    sdg.swap(sdg);

    size_type i1 = sdg.number_of_input_sites();
    size_type n1 = sdg.number_of_vertices();

    sdg2.clear();
    sdg.swap(sdg2);

    assert( sdg2.number_of_input_sites() == i1 );
    assert( sdg2.number_of_vertices() == n1 );
    assert( sdg2.is_valid() );

    assert( sdg.number_of_input_sites() == 0 );
    assert( sdg.number_of_vertices() == 0 );
    assert( sdg.is_valid() );

    //--------------------------------------------

    sdg2.swap(sdg); // back to original setting (before swap)

    sdg2.insert(p4);
    sdg2.insert(p3);
    sdg2.insert(p3, p4);

    i1 = sdg.number_of_input_sites();
    n1 = sdg.number_of_vertices();

    size_type i2 = sdg2.number_of_input_sites();
    size_type n2 = sdg2.number_of_vertices();

    sdg.swap(sdg2);

    assert( sdg2.number_of_input_sites() == i1 );
    assert( sdg2.number_of_vertices() == n1 );
    assert( sdg2.is_valid() );

    assert( sdg.number_of_input_sites() == i2 );
    assert( sdg.number_of_vertices() == n2 );
    assert( sdg.is_valid() );

    //--------------------------------------------

    i2 = sdg2.number_of_input_sites();
    n2 = sdg2.number_of_vertices();

    sdg.clear();
    sdg.swap(sdg2);

    assert( sdg2.number_of_input_sites() == 0 );
    assert( sdg2.number_of_vertices() == 0 );
    assert( sdg2.is_valid() );

    assert( sdg.number_of_input_sites() == i2 );
    assert( sdg.number_of_vertices() == n2 );
    assert( sdg.is_valid() );
  }
  end_testing("swap method");

  Point_2 p5(0.625,0.125), p6(0.125,0.625);
  Site_2 s5 = Site_2::construct_site_2(p5, p6);

  sdg.insert(s5);

  Site_2 s6 = Site_2::construct_site_2(p1, p4);

  sdg.insert(s6);

  start_testing("file I/O methods and I/O operators");
  {
    size_type nisc1 = sdg.number_of_input_sites();
    size_type npc1 = sdg.point_container().size();
    size_type nv1 = sdg.number_of_vertices();

    std::ofstream ofs(ofname);
    assert( ofs );
    sdg.file_output(ofs);
    assert( sdg.is_valid() );
    ofs.close();

    sdg.clear();

    std::ifstream ifs(ofname);
    assert( ifs );
    sdg.file_input(ifs);
    assert( sdg.is_valid() );
    ifs.close();

    size_type nisc2 = sdg.number_of_input_sites();
    size_type npc2 = sdg.point_container().size();
    size_type nv2 = sdg.number_of_vertices();
    if ( nisc1 != nisc2 || npc1 != npc2 || nv1 != nv2 ) { return false; }
  }
  end_testing("file I/O methods and I/O operators");

  start_testing("validity check and clear methods");
  std::cout << std::endl;
  {
    std::cout << "  validating: " << std::flush;
    sdg.is_valid(true, 1);
    std::cout << std::endl;
    std::cout << "  clearing diagram..." << std::endl;
    sdg.clear();
    std::cout << "  validating: " << std::flush;
    sdg.is_valid(true, 1);
    std::cout << std::endl;
  }
  start_testing("validity check and clear methods");
  end_testing("validity check and clear method");


  return true;
}


} //namespace CGAL



#endif // CGAL_SDG_TEST_TYPES_H

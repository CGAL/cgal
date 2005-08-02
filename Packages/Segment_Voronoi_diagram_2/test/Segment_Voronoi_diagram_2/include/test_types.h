#ifndef CGAL_SVD_TEST_TYPES_H
#define CGAL_SVD_TEST_TYPES_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <iostream>
#include <algorithm>
#include <cassert>

#include "IO/Null_output_stream.h"
#include "IO/io_aux.h"

CGAL_BEGIN_NAMESPACE

template<class SVD>
struct Level_finder
{
  typedef typename SVD::size_type      size_type;
  typedef typename SVD::Vertex_handle  Vertex_handle;

  size_type operator()(Vertex_handle v) const {
    CGAL_precondition( v != Vertex_handle() );
    size_type level = 0;
    Vertex_handle vertex = v;
    while ( vertex->up() != Vertex_handle() ) {
      vertex = vertex->up();
      level++;
    }

    return level;
  }
};

template<class Gt, class SVDDS, class LTag>
struct Level_finder< Segment_Voronoi_diagram_2<Gt,SVDDS,LTag> >
{
  typedef Segment_Voronoi_diagram_2<Gt,SVDDS,LTag> SVD;

  typedef typename SVD::size_type      size_type;
  typedef typename SVD::Vertex_handle  Vertex_handle; 

  size_type operator()(Vertex_handle v) const { return 0; }
};


template<class SVD, class InputStream>
bool test_svd(InputStream& is, const SVD&, char* fname)
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
  //  svd.insert(t5);
  //  svd2.clear();
  //  std::cout << std::endl;
  {
    Input_sites_iterator isi = svd.input_sites_begin();
    for (; isi != svd.input_sites_end(); ++isi) {
      null_os << *isi;
      //      std::cout << *isi << std::endl;
    }
    //    std::cout << std::endl;

    Output_sites_iterator osi = svd.output_sites_begin();
    for (; osi != svd.output_sites_end(); ++osi) {
      null_os << *osi;
      //      std::cout << *osi << std::endl;
    }
    //    std::cout << std::endl;
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

  start_testing("removal methods");
  {
    svd.clear();

    std::ifstream ifs("data/bizarre.cin");
    assert( ifs );
    Site_2 t;
    while ( ifs >> t ) {
      svd.insert(t);
    }
    std::cerr << std::endl;
    CGAL_assertion( svd.is_valid(true, 1) );

    Finite_vertices_iterator vit;

    std::vector<Vertex_handle> vec;
    for (vit = svd.finite_vertices_begin();
	 vit != svd.finite_vertices_end(); ++vit) {
      vec.push_back(vit);
    }
    std::random_shuffle(vec.begin(), vec.end());

    typename std::vector<Vertex_handle>::iterator it = vec.begin();
    std::cerr << std::endl;
    Level_finder<SVD> level;
    do {
      Site_2 tt = (*it)->site();
      std::cerr << "  *** attempting to remove: " << tt << "\t" << std::flush;
      if ( tt.is_point() ) { std::cerr << "\t" << std::flush; }
      std::cerr << "  - LEVEL of v: " << level(*it)
		<< " -  " << std::flush;
      bool success = svd.remove(*it);
      std::cerr << (success ? " successful" : " UNSUCCESSFUL") << std::endl;
      if ( success ) {
	vec.erase(it);
	it = vec.begin();
      } else {
	++it;
      }
    } while ( it != vec.end() );

    // second test case
    svd.clear();

    svd.insert(site_list.begin(), site_list.end());
    CGAL_assertion( svd.is_valid() );

    vit = svd.finite_vertices_begin();
    do {
      bool success = svd.remove(vit);
      if ( success ) {
	vit = svd.finite_vertices_begin();
      } else {
	++vit;
      }
    } while ( vit != svd.finite_vertices_end() );
  }
  end_testing("removal methods");

  {
    // recover state of svd after testing removals
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

  Point_2 p5(0.625,0.125), p6(0.125,0.625);
  Site_2 s5 = Site_2::construct_site_2(p5, p6);

  svd.insert(s5);

  Site_2 s6 = Site_2::construct_site_2(p1, p4);

  svd.insert(s6);

  start_testing("file I/O methods and I/O operators");
  {
    size_type nisc1 = svd.number_of_input_sites();
    size_type npc1 = svd.point_container().size();
    size_type nv1 = svd.number_of_vertices();

    std::ofstream ofs(fname);
    assert( ofs );
    svd.file_output(ofs);
    CGAL_assertion( svd.is_valid() );
    ofs.close();

    svd.clear();

    std::ifstream ifs(fname);
    assert( ifs );
    svd.file_input(ifs);
    CGAL_assertion( svd.is_valid() );
    ifs.close();

    size_type nisc2 = svd.number_of_input_sites();
    size_type npc2 = svd.point_container().size();
    size_type nv2 = svd.number_of_vertices();
    if ( nisc1 != nisc2 || npc1 != npc2 || nv1 != nv2 ) { return false; }
  }
  end_testing("file I/O methods and I/O operators");

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


CGAL_END_NAMESPACE



#endif // CGAL_SVD_TEST_TYPES_H

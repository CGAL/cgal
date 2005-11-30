// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef VDA_TEST_LOCATE_H
#define VDA_TEST_LOCATE_H 1

#include <CGAL/basic.h>
#include <iostream>
#include <vector>
#include <boost/variant.hpp>

#include "helper_functions.h"

template<class VDA, class Projector, class QStream, class OStream>
void test_locate(const VDA& vda, const Projector& project,
		 QStream& isq, OStream& os = std::cout,
		 bool print_sites)
{
  std::cout << std::endl;
  std::cout << "is Delaunay graph valid? "
	    << (vda.dual().is_valid() ? "yes" : "no") << std::endl;
  std::cout << std::endl;

  os << "Dimension of Delaunay graph: " << vda.dual().dimension()
     << std::endl << std::endl;

  os << "Vertices of the Delaunay graph:" << std::endl;
  typename VDA::Delaunay_graph::Finite_vertices_iterator vit;
  for (vit = vda.dual().finite_vertices_begin();
       vit != vda.dual().finite_vertices_end(); ++vit) {
    os << project(vit) << std::endl;
  }
  os << std::endl;

  typedef typename VDA::Voronoi_traits::Point_2  Point_2;
  Point_2 p;
  std::vector<Point_2>  vecp;

  while ( isq >> p ) {
    vecp.push_back(p);
  }

  test_locate_dg(vda, project, vecp, os);
  test_locate_vd(vda, vecp, os, print_sites);
}

template<class VDA, class Projector, class Point_vector, class OStream>
void test_locate_dg(const VDA& vda, const Projector& project,
		    const Point_vector& vecp, OStream& os)
{
  typedef typename VDA::Voronoi_traits                Voronoi_traits;
  typedef typename Voronoi_traits::Nearest_site_2     Nearest_site_2;
  typedef typename Voronoi_traits::Vertex_handle      Vertex_handle;
  typedef typename Voronoi_traits::Face_handle        Face_handle;
  typedef typename Voronoi_traits::Edge               Edge;
  typedef typename Nearest_site_2::result_type        result_type;

  Nearest_site_2 nearest_site = vda.voronoi_traits().nearest_site_2_object();

  result_type ns_qr;

  os << "Query sites and location feature in dual graph:" << std::endl;
  for (unsigned int i = 0; i < vecp.size(); ++i) {
    os << vecp[i] << "\t --> \t" << std::flush;
    ns_qr = nearest_site(vda.dual(), vecp[i]);
    if ( Vertex_handle* v = boost::get<Vertex_handle>(&ns_qr) ) {
      os << "FACE";
      kill_warning( v );
    } else if ( Edge* e = boost::get<Edge>(&ns_qr) ) {
      os << "EDGE";
      kill_warning( e );
    } else if ( Face_handle* f = boost::get<Face_handle>(&ns_qr) ) {
      os << "VERTEX";
      kill_warning( f );
    } else {
      os << " *** NOT READY YET *** ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

template<class VDA, class Point_vector, class OStream>
void test_locate_vd(const VDA& vda, const Point_vector& vecp,
		    OStream& os, bool print_sites)
{
  typename VDA::Locate_result lr;
  typename VDA::Voronoi_traits::Access_site_2 get_site =
    vda.voronoi_traits().access_site_2_object();
  os << "Query sites and location feature in dual graph:" << std::endl;
  for (unsigned int i = 0; i < vecp.size(); ++i) {
    os << vecp[i] << "\t --> \t" << std::flush;
    lr = vda.locate(vecp[i]);
    if ( lr.is_edge() ) {
      typename VDA::Halfedge_handle e = lr;
      os << "VORONOI EDGE";
      if ( print_sites ) {
	os << " ---> ";
	os << " up: " << get_site(e->up());
	os << " down: " << get_site(e->down());
	os << " left: ";
	if ( e->has_source() ) {
	  os << get_site(e->left());
	} else {
	  os << "inf";
	}
	os << " right: ";
	if ( e->has_target() ) {
	  os << get_site(e->right());
	} else {
	  os << "inf";
	}
      } // if ( print_sites )
      kill_warning( e );
    } else if ( lr.is_vertex() ) {
      os << "VORONOI VERTEX";
      typename VDA::Vertex_handle v = lr;
      if ( print_sites ) {
	os << " ---> ";
	for (int i = 0; i < 3; ++i) {
	  os << get_site(v->dual()->vertex(i)) << " ";
	}
      }
      kill_warning( v );
    } else if ( lr.is_face() ) {
      typename VDA::Face_handle f = lr;
      kill_warning( f );
      os << "VORONOI FACE";
      if ( print_sites ) {
	os << " ---> " << get_site(f->dual());
      }
    } else {
      os << " *** NOT READY YET *** ";
      CGAL_assertion( false );
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}


#endif // VDA_TEST_LOCATE_H

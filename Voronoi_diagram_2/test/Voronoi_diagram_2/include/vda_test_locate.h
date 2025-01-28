// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef VDA_TEST_LOCATE_H
#define VDA_TEST_LOCATE_H 1

#include <iostream>
#include <vector>
#include <variant>

#include "helper_functions.h"


template<class VDA, class Projector, class Point_vector, class OStream>
void test_locate_dg(const VDA& vda, const Projector& ,
                    const Point_vector& vecp, OStream& os);

template<class VDA, class Point_vector, class OStream>
void test_locate_vd(const VDA& vda, const Point_vector& vecp,
                    OStream& os, bool print_sites);



template<class VDA, class Projector, class QStream, class OStream>
void test_locate(const VDA& vda, const Projector& project,
                 QStream& isq, OStream& os, bool print_sites)
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

  typedef typename VDA::Adaptation_traits::Point_2  Point_2;
  Point_2 p;
  std::vector<Point_2>  vecp;

  while ( isq >> p ) {
    vecp.push_back(p);
  }

  test_locate_dg(vda, project, vecp, os);
  test_locate_vd(vda, vecp, os, print_sites);
}

template<class VDA, class Projector, class Point_vector, class OStream>
void test_locate_dg(const VDA& vda, const Projector& ,
                    const Point_vector& vecp, OStream& os)
{
  typedef typename VDA::Adaptation_traits                     Adaptation_traits;
  typedef typename Adaptation_traits::Nearest_site_2          Nearest_site_2;
  typedef typename Adaptation_traits::Delaunay_vertex_handle  Vertex_handle;
  typedef typename Adaptation_traits::Delaunay_face_handle    Face_handle;
  typedef typename Adaptation_traits::Delaunay_edge           Edge;
  typedef typename Nearest_site_2::result_type                result_type;

  Nearest_site_2 nearest_site = vda.adaptation_traits().nearest_site_2_object();

  result_type ns_qr;

  os << "Query sites and location feature in dual graph:" << std::endl;
  for (unsigned int i = 0; i < vecp.size(); ++i) {
    os << vecp[i] << "\t --> \t" << std::flush;
    ns_qr = nearest_site(vda.dual(), vecp[i]);
    if ( Vertex_handle* v = std::get_if<Vertex_handle>(&ns_qr) ) {
      os << "FACE";
      kill_warning( v );
    } else if ( Edge* e = std::get_if<Edge>(&ns_qr) ) {
      os << "EDGE";
      kill_warning( e );
    } else if ( Face_handle* f = std::get_if<Face_handle>(&ns_qr) ) {
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
  typedef typename VDA::Locate_result      Locate_result;
  typedef typename VDA::Face_handle        Face_handle;
  typedef typename VDA::Vertex_handle      Vertex_handle;
  typedef typename VDA::Halfedge_handle    Halfedge_handle;

  Locate_result lr;
  typename VDA::Adaptation_traits::Access_site_2 get_site =
    vda.adaptation_traits().access_site_2_object();
  os << "Query sites and location feature in dual graph:" << std::endl;
  for (unsigned int i = 0; i < vecp.size(); ++i) {
    os << vecp[i] << "\t --> \t" << std::flush;
    lr = vda.locate(vecp[i]);
    if ( Halfedge_handle* ee = std::get_if<Halfedge_handle>(&lr) ) {
      Halfedge_handle e = *ee;
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
    } else if ( Vertex_handle* vv = std::get_if<Vertex_handle>(&lr) ) {
      os << "VORONOI VERTEX";
      Vertex_handle v = *vv;
      if ( print_sites ) {
        os << " ---> ";
        for (int i = 0; i < 3; ++i) {
          os << get_site(v->dual()->vertex(i)) << " ";
        }
      }
      kill_warning( v );
    } else if ( Face_handle* ff = std::get_if<Face_handle>(&lr) ) {
      typename VDA::Face_handle f = *ff;
      kill_warning( f );
      os << "VORONOI FACE";
      if ( print_sites ) {
        os << " ---> " << get_site(f->dual());
      }
    } else {
      os << " *** NOT READY YET *** ";
      CGAL_error();
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}


#endif // VDA_TEST_LOCATE_H

// Copyright (c) 1997-2000  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Thomas Herrmann, Lutz Kettner

#ifndef CGAL_WIDTH_POLYHEDRON_3_H
#define CGAL_WIDTH_POLYHEDRON_3_H

#include <CGAL/license/Polytope_distance_d.h>


#include <CGAL/basic.h>
#include <CGAL/HalfedgeDS_vertex_base.h>
#include <CGAL/HalfedgeDS_halfedge_base.h>
#include <CGAL/HalfedgeDS_face_base.h>
#include <map>
#include <CGAL/width_assertions.h>

namespace CGAL {



template <class Refs, class Traits>
class Width_vertex_default_base 
    : public CGAL::HalfedgeDS_vertex_base< 
        Refs, Tag_true, typename Traits::Point_3> {
private:
    typedef Traits          WT;
    typedef typename WT::RT RT;
    WT tco;
public:
    typedef typename Traits::Point_3 Point_3;
    typedef CGAL::HalfedgeDS_vertex_base< Refs, Tag_true, Point_3> Vertex_base;
    Width_vertex_default_base() {}
    Width_vertex_default_base( const Point_3& p) : Vertex_base(p) {}
};

struct Width_polyhedron_items_3 {
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef Width_vertex_default_base< Refs, Traits>  Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper {
        typedef CGAL::HalfedgeDS_halfedge_base< Refs>     Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef typename Traits::Plane_3 Plane;
        typedef CGAL::HalfedgeDS_face_base< Refs, CGAL::Tag_true, Plane>  Face;
    };
};


namespace Width_3_internal {

template <class InputPolyhedron, class Width_Traits>
class Data_access {
 public:
  typedef InputPolyhedron Polyhedron;
  typedef typename Polyhedron::Vertex                Vertex;
  typedef typename Polyhedron::Facet                 Facet;
  typedef typename Polyhedron::Halfedge              Halfedge;
  typedef typename Polyhedron::Facet_handle          Facet_handle;
  typedef typename Polyhedron::Vertex_handle         Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle       Halfedge_handle;
  typedef typename Polyhedron::Facet_const_handle    Facet_const_handle;
  typedef typename Polyhedron::Vertex_const_handle   Vertex_const_handle;
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Polyhedron::Point_3               PolyPoint;
  typedef typename Polyhedron::Plane_3               Plane;
  typedef typename Polyhedron::Vertex_iterator       Vertex_iterator;
  typedef typename Polyhedron::Facet_iterator        Facet_iterator;
  typedef typename Polyhedron::Halfedge_iterator     Halfedge_iterator;
  typedef typename Width_Traits::RT                  RT;
 private:
  Width_Traits tco;
 public:
  Data_access() {}
  ~Data_access() {}

 private:
  //Precondition: Plane Equation already computed in a deterministic way
  struct Facet_compare {
    bool operator()(Facet_const_handle f, 
		    Facet_const_handle g) const {
      Width_Traits tco;
      Plane fpp=f->plane();
      Plane gpp=g->plane();
      RT fa,fb,fc,fd,ga,gb,gc,gd;
      tco.get_plane_coefficients(fpp,fa,fb,fc,fd);
      tco.get_plane_coefficients(gpp,ga,gb,gc,gd);
      return fa<ga 
	      || ( fa==ga && fb<gb )
	      || ( fa==ga && fb==gb && fc<gc )
	      || ( fa==ga && fb==gb && fc==gc && fd<gd );
    }
  };
  //Precondition: Plane Equation already computed in a deterministic way
  //and a facet is bounded by exactly 3 edges! 
  struct Halfedge_compare {
    bool operator()(Halfedge_const_handle e, 
		    Halfedge_const_handle h) const {
      Width_Traits tco;
      PolyPoint etail=e->opposite()->vertex()->point();
      PolyPoint ehead=e->vertex()->point();
      PolyPoint htail=h->opposite()->vertex()->point();
      PolyPoint hhead=h->vertex()->point();
      RT etx,ety,etz,eth,ehx,ehy,ehz,ehh;
      tco.get_point_coordinates(etail,etx,ety,etz,eth);
      tco.get_point_coordinates(ehead,ehx,ehy,ehz,ehh);
      RT htx,hty,htz,hth,hhx,hhy,hhz,hhh;
      tco.get_point_coordinates(htail,htx,hty,htz,hth);
      tco.get_point_coordinates(hhead,hhx,hhy,hhz,hhh);
      
      return (etx*hth <htx*eth ||
	      ( etx*hth==htx*eth && ety*hth <hty*eth ) ||
	      ( etx*hth==htx*eth && ety*hth==hty*eth && etz*hth <htz*eth ) ||
	      ( etx*hth==htx*eth && ety*hth==hty*eth && etz*hth==htz*eth
	      && ehx*hhh <hhx*ehh ) ||
	      ( etx*hth==htx*eth && ety*hth==hty*eth && etz*hth==htz*eth 
	      && ehx*hhh==hhx*ehh && ehy*hhh <hhy*ehh ) ||
	      ( etx*hth==htx*eth && ety*hth==hty*eth && etz*hth==htz*eth 
	      && ehx*hhh==hhx*ehh && ehy*hhh==hhy*ehh && ehz*hhh<hhz*ehh )
	      );
    }
  };

  std::map< Facet_const_handle, 
    std::vector<Vertex_handle>,
    Facet_compare> 
    antipodal_vertices;

  std::map< Halfedge_const_handle, 
    bool, 
    Halfedge_compare> 
    visited_halfedges;

  std::map< Halfedge_const_handle, 
    bool, 
    Halfedge_compare> 
    impassable_halfedges;
  
 public:
  bool is_visited(Halfedge_handle& e) const {
    typename std::map < Halfedge_const_handle, bool, 
      Halfedge_compare > ::const_iterator it;
    it=visited_halfedges.find(e);
    CGAL_assertion(it!=visited_halfedges.end());
    DEBUGENDL(VISITED_CHECK,"Visited flag value of edge "
	      <<e->opposite()->vertex()->point()
	      <<" --> ",e->vertex()->point()<<": "<<it->second);
    return it->second;
  }
  void set_visited_flag(Halfedge_handle& e, bool val) {
    DEBUGENDL(VISITED_CHECK,"Set visited flag to: ",val);
    visited_halfedges[e]=val;
  }
  
  bool is_impassable(Halfedge_handle& e) const {
    typename std::map < Halfedge_const_handle, 
      bool, 
      Halfedge_compare > ::const_iterator it;
    it=impassable_halfedges.find(e);
    CGAL_assertion(it!=impassable_halfedges.end());
    DEBUGENDL(IMPASSABLE_CHECK,"Impassable flag value of edge ",
	      e->opposite()->vertex()->point()
	      <<" --> "<<e->vertex()->point()<<": "<<it->second);
    return it->second;
  }
  void set_impassable_flag(Halfedge_handle& e, bool val) {
    DEBUGENDL(IMPASSABLE_CHECK,"Set impassable flag to: ",val);
    impassable_halfedges[e]=val;
  }
  
  void set_antipodal_vertices(Facet_handle& f, 
				 std::vector<Vertex_handle>& V) {
    antipodal_vertices[f]=V;
  }
  
  void get_antipodal_vertices(Facet_handle& f, 
			      std::vector<Vertex_handle>& res) const {
    typename std::map< Facet_const_handle, 
      std::vector<Vertex_handle>, Facet_compare>::const_iterator it;
    it=antipodal_vertices.find(f);
    CGAL_assertion(it!=antipodal_vertices.end());
    res=it->second;
  }

#if !(defined(CGAL_KERNEL_NO_ASSERTIONS) || defined(CGAL_NO_ASSERTIONS) \
      || defined(NDEBUG))
  int size_of_impassable() {
    return(int(impassable_halfedges.size()));
  }
  int size_of_visited() {
    return(int(visited_halfedges.size()));
  }
  int size_of_antipodal_vertices() {
    return (int(antipodal_vertices.size()));
  }
#endif
};

} // namespace Width_3_internal

} //namespace CGAL

#endif //WIDTH_POLYHEDRON_3_H

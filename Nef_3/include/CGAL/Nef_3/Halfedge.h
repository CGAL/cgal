// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
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
// 
//
// Author(s)     : Michael Seel        <seel@mpi-sb.mpg.de>
//                 Miguel Granados     <granados@mpi-sb.mpg.de>
//                 Susan Hert          <hert@mpi-sb.mpg.de>
//                 Lutz Kettner        <kettner@mpi-sb.mpg.de>
//                 Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>
#ifndef CGAL_NEF_HALFEDGE_H
#define CGAL_NEF_HALFEDGE_H

#include <CGAL/license/Nef_3.h>


#include <string>
#include <sstream>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Nef_3/SNC_iteration.h>
#include <CGAL/Origin.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 83
#include <CGAL/Nef_2/debug.h>

#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

namespace CGAL {

template <typename Refs>
class Halfedge_base
{ // == Halfedge
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  typedef void* GenPtr;
  #else
  typedef boost::any GenPtr;
  #endif
  typedef typename Refs::Mark  Mark;
  typedef typename Refs::Vector_3  Vector_3;
  typedef typename Refs::Sphere_point  Sphere_point;
  typedef typename Refs::Vertex_handle    Vertex_handle;
  typedef typename Refs::Halfedge_handle  Halfedge_handle;
  typedef typename Refs::SVertex_handle   SVertex_handle;
  typedef typename Refs::SHalfedge_handle SHalfedge_handle;
  typedef typename Refs::SFace_handle     SFace_handle;
  typedef typename Refs::Vertex_const_handle    Vertex_const_handle;
  typedef typename Refs::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Refs::SVertex_const_handle   SVertex_const_handle;
  typedef typename Refs::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename Refs::SFace_const_handle     SFace_const_handle;

  Vertex_handle      center_vertex_;
  Mark               mark_;
  SVertex_handle     twin_;              
  SHalfedge_handle   out_sedge_;           
  SFace_handle       incident_sface_;    
  GenPtr             info_;                 
  Sphere_point       point_;   

 public:

  Halfedge_base() : center_vertex_(), mark_(), twin_(),
    out_sedge_(), incident_sface_(),
    info_(), point_() {}

    Halfedge_base(Mark m) :  center_vertex_(), mark_(m), twin_(),
      out_sedge_(), incident_sface_(),
      info_(), point_() {}

      ~Halfedge_base() {
	CGAL_NEF_TRACEN("  destroying Halfedge item "<<&*this);
      }

      Halfedge_base(const Halfedge_base<Refs>& e) 
	{ center_vertex_ = e.center_vertex_;
	  point_ = e.point_;
	  mark_ = e.mark_;
	  twin_ = e.twin_;
	  out_sedge_ = e.out_sedge_;
	  incident_sface_ = e.incident_sface_;
	  info_ = 0;
	}

      Halfedge_base<Refs>& operator=(const Halfedge_base<Refs>& e) 
	{ center_vertex_ = e.center_vertex_;
	  point_ = e.point_;
	  mark_ = e.mark_;
	  twin_ = e.twin_;
	  out_sedge_ = e.out_sedge_;
	  incident_sface_ = e.incident_sface_;
	  info_ = 0;
	  return *this;
	}

      Vertex_handle& center_vertex() { return center_vertex_; }
      Vertex_const_handle center_vertex() const { return center_vertex_; }

      Vertex_handle& source() { return center_vertex_; }
      Vertex_const_handle source() const { return center_vertex_; }

      Vertex_handle& target() { return twin()->source(); }
      Vertex_const_handle target() const { return twin()->source(); }

      Mark& mark() { return mark_; }
      const Mark& mark() const { return mark_; }

      Vector_3 vector() const { return (point_ - CGAL::ORIGIN); }
      Sphere_point& point(){ return point_; }
      const Sphere_point& point() const { return point_; }

      SVertex_handle& twin() { return twin_; }
      SVertex_const_handle twin()  const { return twin_; }

      SHalfedge_handle& out_sedge() { return out_sedge_; }
      SHalfedge_const_handle out_sedge() const { return out_sedge_; }

      SFace_handle& incident_sface() { return incident_sface_; } 
      SFace_const_handle incident_sface() const { return incident_sface_; } 

      bool is_isolated() const { return (out_sedge() == SHalfedge_handle()); }

      GenPtr& info() { return info_; }
      const GenPtr& info() const { return info_; }

 public:
      std::string debug() const
	{ std::stringstream os; 
	  set_pretty_mode(os);
	  os<<"sv [ "<<point_
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      <<info_
    #endif
      <<" ] ";
	  return os.str();
	}

      bool is_twin() const { return (&*twin_ < this); }

      bool is_valid( bool verb = false, int level = 0) const {
      
	Verbose_ostream verr(verb);
	verr << "begin CGAL::SNC_items<...>::Halfedge_base::is_valid( verb=true, "
	  "level = " << level << "):" << std::endl;

	bool valid = (center_vertex_ != NULL && center_vertex_ != Vertex_handle());
	valid = valid && (twin_ != NULL && twin_ != SVertex_handle() &&
			  twin_ != SVertex_handle());
	//      valid = valid && (out_sedge_ != NULL);
	//      valid = valid && (incident_sface_ != SFace_handle());
      
	//      valid = valid &&((out_sedge_ != NULL && incident_sface_ == NULL) ||
	//		       (out_sedge_ == NULL && incident_sface_ != NULL));
      
	valid = valid && (out_sedge_ != NULL || incident_sface_ != NULL);

	verr << "end of CGAL::SNC_items<...>::Halfedge_base::is_valid(): structure is "
	     << ( valid ? "valid." : "NOT VALID.") << std::endl;

	return valid;
      }

}; // Halfedge_base

} //namespace CGAL
#endif //CGAL_NEF_HALFEDGE_H

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

#ifndef CGAL_SM_ITEMS_H
#define CGAL_SM_ITEMS_H

#include <CGAL/basic.h>
#include <CGAL/In_place_list.h>
#include <CGAL/Object.h>
#include <CGAL/Nef_S2/Sphere_geometry.h>
#include <string>
#include <sstream>
#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

namespace CGAL {

//template <typename K, typename I,typename C> class Sphere_map;
template <typename SM> class SM_const_decorator;
template <typename SM> class SM_decorator;

struct SM_items {
public:

  template <typename Refs>
  class SVertex
  { 
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    typedef void*  GenPtr;
    #else
    typedef boost::any GenPtr;
    #endif
    typedef typename Refs::Mark                    Mark;
    typedef typename Refs::Sphere_point            Sphere_point;
    typedef typename Refs::SVertex_handle          SVertex_handle;
    typedef typename Refs::SVertex_const_handle    SVertex_const_handle;
    typedef typename Refs::SHalfedge_handle        SHalfedge_handle;    
    typedef typename Refs::SHalfedge_const_handle  SHalfedge_const_handle;
    typedef typename Refs::SFace_handle            SFace_handle;
    typedef typename Refs::SFace_const_handle      SFace_const_handle;

    Sphere_point    point_; 
    Mark            mark_;
    SHalfedge_handle out_sedge_;
    SFace_handle     incident_sface_;
    GenPtr          info_;
    // temporary information:

  public:
    SVertex() : 
      point_(), mark_(), out_sedge_(), incident_sface_(), info_() {}
    SVertex(const Mark& m) : 
      point_(), mark_(m), out_sedge_(), incident_sface_(), info_() {}
    SVertex(const Sphere_point& p) : 
      point_(p), mark_(), out_sedge_(), incident_sface_(), info_() {}

    ~SVertex() {}

    SVertex(const SVertex<Refs>& v)
    { point_ = v.point_;
      mark_ = v.mark_;
      out_sedge_ = v.out_sedge_;
      incident_sface_ = v.incident_sface_;
      info_ = 0;
    }

    SVertex<Refs>& operator=(const SVertex<Refs>& v)
    { point_ = v.point_;
      mark_ = v.mark_;
      out_sedge_ = v.out_sedge_;
      incident_sface_ = v.incident_sface_;
      info_ = 0;
      return *this;
    }

    Mark& mark() { return mark_; }
    const Mark& mark() const { return mark_; }

    Sphere_point& point(){ return point_; }
    const Sphere_point& point() const { return point_; }

    SHalfedge_handle& out_sedge() { return out_sedge_; }
    SHalfedge_const_handle out_sedge() const { return out_sedge_; }

    SFace_handle& incident_sface() { return incident_sface_; } 
    SFace_const_handle incident_sface() const { return incident_sface_; } 

    bool is_isolated() const { return (out_sedge() == SHalfedge_handle()); }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }
                          
    public:
    std::string debug() const
    { std::ostringstream os; set_pretty_mode(os);
      os<<"V"<<point_
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
        <<' '<<info_
      #endif
      <<'\0';
      std::string res(os.str()); return res;
    }
 
  }; // SVertex    


  template <typename Refs>
  class SHalfedge
  { 
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    typedef void*  GenPtr;
    #else
    typedef boost::any GenPtr;
    #endif
    typedef typename Refs::Mark             Mark;
    typedef typename Refs::Sphere_circle    Sphere_circle;
    typedef typename Refs::SVertex_handle          SVertex_handle;
    typedef typename Refs::SVertex_const_handle    SVertex_const_handle;
    typedef typename Refs::SHalfedge_handle        SHalfedge_handle;    
    typedef typename Refs::SHalfedge_const_handle  SHalfedge_const_handle;
    typedef typename Refs::SFace_handle            SFace_handle;
    typedef typename Refs::SFace_const_handle      SFace_const_handle;

    // Role within local graph:
    Sphere_circle      circle_;
    Mark               mark_;
    SHalfedge_handle   twin_, sprev_, snext_;
    SVertex_handle     source_;
    SFace_handle       incident_sface_;
    GenPtr             info_;

  public:
    SHalfedge() : circle_(), mark_(), twin_(), sprev_(), snext_(),
		 source_(), incident_sface_(), info_() {}

    ~SHalfedge() {}

    SHalfedge(const SHalfedge<Refs>& e)
    {
      circle_ = e.circle_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      sprev_ = e.sprev_;
      snext_ = e.snext_;
      source_ = e.source_;
      incident_sface_ = e.incident_sface_;
      info_ = 0;
    }
    SHalfedge<Refs>& operator=(const SHalfedge<Refs>& e)
    {
      circle_ = e.circle_;
      mark_ = e.mark_;
      twin_ = e.twin_;
      sprev_ = e.sprev_;
      snext_ = e.snext_;
      source_ = e.source_;
      incident_sface_ = e.incident_sface_;
      info_ = 0;
      return *this;
    }

    bool is_twin() const { return (&*twin_ < this); }

    Mark& mark() { 
      return mark_;
    }
    const Mark& mark() const {
      return mark_;
    }

    SHalfedge_handle& twin() { return twin_; }
    SHalfedge_const_handle twin() const { return twin_; }

    SVertex_handle& source() { return source_; }
    SVertex_const_handle source() const { return source_; }

    SVertex_handle& target() { return twin()->source(); }
    SVertex_const_handle target() const { return twin()->source(); }

    SHalfedge_handle& sprev() { return sprev_; }
    SHalfedge_const_handle sprev() const { return sprev_; }

    SHalfedge_handle& snext() { return snext_; }
    SHalfedge_const_handle snext() const { return snext_; }

    Sphere_circle& circle() { return circle_; }
    const Sphere_circle& circle() const { return circle_; }
    
    SFace_handle& incident_sface() { return incident_sface_; }
    SFace_const_handle incident_sface() const { return incident_sface_; }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

    std::string debug() const
    { std::ostringstream os; set_pretty_mode(os); 
      os <<"e["<<source_->debug()<<", "
         <<twin_->source_->debug()<<
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
      " "<< info_ <<
      #endif
      "]"<<'\0';
      std::string res(os.str()); return res;
    }

  }; // SHalfedge

  template <typename Refs> 
  class SHalfloop
  {
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    typedef void*  GenPtr;
    #else
    typedef boost::any GenPtr;
    #endif
    typedef typename Refs::Mark             Mark;
    typedef typename Refs::Sphere_circle    Sphere_circle;
    typedef typename Refs::SHalfloop_handle        SHalfloop_handle;    
    typedef typename Refs::SHalfloop_const_handle  SHalfloop_const_handle;
    typedef typename Refs::SFace_handle            SFace_handle;
    typedef typename Refs::SFace_const_handle      SFace_const_handle;


    Sphere_circle   circle_;
    Mark            mark_;
    SHalfloop_handle twin_;
    SFace_handle     incident_sface_;
    GenPtr          info_;
    // temporary needed:

  public:
    SHalfloop() : circle_(), mark_(), twin_(), incident_sface_(), info_() {}
    ~SHalfloop() {}
    SHalfloop(const SHalfloop<Refs>& l)
    { 
      circle_ = l.circle_;
      mark_ = l.mark_;
      twin_ = l.twin_;
      incident_sface_ = l.incident_sface_;
      info_ = 0;
    }
    SHalfloop<Refs>& operator=(const SHalfloop<Refs>& l)
    { 
      circle_ = l.circle_;
      mark_ = l.mark_;
      twin_ = l.twin_;
      incident_sface_ = l.incident_sface_;
      info_ = 0;
      return *this;
    }

    bool is_twin() const { return (&*twin_ < this); }

    Mark& mark() { 
      return mark_;
    }

    const Mark& mark() const { 
      return mark_;
    }

    SHalfloop_handle& twin() { return twin_; }
    SHalfloop_const_handle twin() const { return twin_; }

    Sphere_circle& circle() { return circle_; }
    const Sphere_circle& circle() const { return circle_; }

    SFace_handle& incident_sface() { return incident_sface_; }
    SFace_const_handle incident_sface() const { return incident_sface_; }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

    std::string debug() const
    { std::ostringstream os; set_pretty_mode(os); 
      os<<"l"<<circle_
      #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
        <<' '<<info_
      #endif
        <<'\0';
      std::string res(os.str()); return res;
    }

  }; // SHalfloop

  template <typename Refs>
  class SFace
  { 
    #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
    typedef void*  GenPtr;
    #else
    typedef boost::any GenPtr;
    #endif
    typedef typename Refs::Mark                  Mark;
    typedef typename Refs::Object_handle         Object_handle;
    typedef typename Refs::Object_list           Object_list;
    typedef typename Refs::SFace_cycle_iterator  SFace_cycle_iterator;
    typedef typename Refs::SFace_cycle_const_iterator 
                           SFace_cycle_const_iterator;

    Mark             mark_;
    Object_list      boundary_entry_objects_; 
    // SHalfedges, SHalfloops, Vertices
    GenPtr           info_;
    // temporary needed:

    public:
    SFace() : mark_(), info_() {}
    ~SFace() {}

    SFace(const SFace<Refs>& f)
    { mark_ = f.mark_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
      info_ = 0;
    }
    SFace<Refs>& operator=(const SFace<Refs>& f)
    { if (this == &f) return *this;
      mark_ = f.mark_;
      boundary_entry_objects_ = f.boundary_entry_objects_;
      info_ = 0;
      return *this;
    }

    SFace_cycle_iterator sface_cycles_begin() 
    { return boundary_entry_objects_.begin(); }
    SFace_cycle_iterator sface_cycles_end()
    { return boundary_entry_objects_.end(); }

    SFace_cycle_const_iterator sface_cycles_begin() const
    { return boundary_entry_objects_.begin(); }
    SFace_cycle_const_iterator sface_cycles_end() const
    { return boundary_entry_objects_.end(); }

    Mark& mark() { return mark_; }
    const Mark& mark() const { return mark_; }

    Object_list& boundary_entry_objects() { return boundary_entry_objects_; }
    const Object_list& boundary_entry_objects() const { return boundary_entry_objects_; }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

  }; // SFace

}; // SM_items


} //namespace CGAL
#endif // CGAL_SM_ITEMS_H

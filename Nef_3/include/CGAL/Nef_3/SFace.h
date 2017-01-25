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
#ifndef CGAL_NEF_SFACE_H
#define CGAL_NEF_SFACE_H

#include <CGAL/license/Nef_3.h>


#include <string>
#include <sstream>
#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/Nef_3/SNC_iteration.h>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 83
#include <CGAL/Nef_2/debug.h>

#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

namespace CGAL {

template <typename Refs>
class SFace_base { 
  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  typedef void* GenPtr;
  #else
  typedef boost::any GenPtr;
  #endif
  typedef typename Refs::Mark  Mark;
  typedef typename Refs::Vertex_handle        Vertex_handle;
  typedef typename Refs::Vertex_const_handle  Vertex_const_handle;
  typedef typename Refs::SFace_handle         SFace_handle;
  typedef typename Refs::SFace_const_handle   SFace_const_handle;
  typedef typename Refs::Volume_handle        Volume_handle;
  typedef typename Refs::Volume_const_handle  Volume_const_handle;
  typedef typename Refs::Object_list          Object_list;
  typedef typename Refs::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename Refs::SFace_cycle_const_iterator 
    SFace_cycle_const_iterator;
  Vertex_handle  center_vertex_;
  Volume_handle  volume_;
  //    Object_list   boundary_entry_objects_; // SEdges, SLoops, SVertices
  GenPtr         info_;
  // temporary needed:
  Mark           mark_;

 public:
  Object_list   boundary_entry_objects_; // SEdges, SLoops, SVertices

  SFace_base() : center_vertex_(), volume_(), info_(), mark_() {}

    ~SFace_base() {
      CGAL_NEF_TRACEN("  destroying SFace_base item "<<&*this);
    }

    SFace_base(const SFace_base<Refs>& f)
      { center_vertex_ = f.center_vertex_;
	volume_ = f.volume_;
	boundary_entry_objects_ = f.boundary_entry_objects_;
	info_ = 0;
	mark_ = f.mark_;
      }

    SFace_base<Refs>& operator=(const SFace_base<Refs>& f)
      { if (this == &f) return *this;
	center_vertex_ = f.center_vertex_;
	volume_ = f.volume_;
	boundary_entry_objects_ = f.boundary_entry_objects_;
	info_ = 0;
	mark_ = f.mark_;
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

    Vertex_handle& center_vertex() { return center_vertex_; }
    Vertex_const_handle center_vertex() const { return center_vertex_; }

    Volume_handle& volume() { return volume_; }
    Volume_const_handle volume() const { return volume_; }

    Object_list& boundary_entry_objects() { return boundary_entry_objects_; }
    const Object_list& boundary_entry_objects() const { return boundary_entry_objects_; }

    GenPtr& info() { return info_; }
    const GenPtr& info() const { return info_; }

    bool is_valid( bool verb = false, int level = 0) const {
      
      Verbose_ostream verr(verb);
      verr << "begin CGAL::SNC_items<...>::SFace_base::is_valid( verb=true, "
	"level = " << level << "):" << std::endl;

      bool valid =(center_vertex_ != Vertex_handle() && center_vertex_ != NULL);
      valid = valid && (volume_ != Volume_handle() &&
			volume_ != NULL);

      if(boundary_entry_objects_.empty()) {
	valid = valid && 
	  (center_vertex_->shalfedges_begin() == center_vertex_->shalfedges_end());
      }
      verr << "end of CGAL::SNC_items<...>::SFace_base::is_valid(): structure is "
	   << ( valid ? "valid." : "NOT VALID.") << std::endl;

      return valid;
    }

}; // SFace_base


} //namespace CGAL
#endif //CGAL_NEF_SFACE_H

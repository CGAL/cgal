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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     :     Peter Hachenberger  <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF_S2_ID_SUPPORT_HANDLER
#define CGAL_NEF_S2_ID_SUPPORT_HANDLER

#include <CGAL/license/Nef_S2.h>


#include <CGAL/Unique_hash_map.h>
#include <map>

#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 131
#include <CGAL/Nef_2/debug.h>

namespace CGAL {

template<typename Items, typename Decorator>
class ID_support_handler {

  typedef typename Decorator::SVertex_handle SVertex_handle;
  typedef typename Decorator::SHalfedge_handle SHalfedge_handle;
  
  typedef typename Decorator::SVertex_const_handle SVertex_const_handle;
  typedef typename Decorator::SHalfedge_const_handle SHalfedge_const_handle;
  typedef typename Decorator::SHalfloop_const_handle SHalfloop_const_handle;
  
 public:
  ID_support_handler() {}

  int get_hash(int) { return 0; }
  template<typename Handle> void initialize_hash(Handle /*h*/) {}
  void initialize_hash(int /*i*/) {}
  void handle_support(SVertex_handle , 
		      SHalfedge_const_handle ,
		      SHalfedge_const_handle ) {} 
  
  void handle_support(SVertex_handle ,
		      SHalfloop_const_handle ,
		      SHalfloop_const_handle ) {}
  
  void handle_support(SVertex_handle,
		      SHalfloop_const_handle,
		      SHalfedge_const_handle) {} 
  
  void handle_support(SVertex_handle,
		      SHalfedge_const_handle,
		      SHalfloop_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SHalfedge_const_handle,
		      SVertex_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SVertex_const_handle,
		      SHalfedge_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SVertex_const_handle,
		      SVertex_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SVertex_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SVertex_const_handle,
		      SHalfloop_const_handle) {}
  
  void handle_support(SVertex_handle,
		      SHalfloop_const_handle,
		      SVertex_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfedge_const_handle,
		      SHalfedge_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfedge_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfloop_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfedge_const_handle,
		      SHalfloop_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfloop_const_handle,
		      SHalfedge_const_handle) {}
  
  void handle_support(SHalfedge_handle,
		      SHalfloop_const_handle,
		      SHalfloop_const_handle) {}
};


} //namespace CGAL
#endif // CGAL_NEF_S2_ID_SUPPORT_HANDLER

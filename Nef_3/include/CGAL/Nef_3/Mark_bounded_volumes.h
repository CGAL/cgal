// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Ralf Osbild <osbild@mpi-sb.mpg.de>

#ifndef CGAL_NEF3_MARK_BOUNDED_VOLUMES_H
#define CGAL_NEF3_MARK_BOUNDED_VOLUMES_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/Nef_polyhedron_3.h>

namespace CGAL {

template<typename Decorator, typename Mark>
class Volume_setter {

  typedef typename Decorator::Vertex_handle 
                                    Vertex_handle;
  typedef typename Decorator::Halfedge_handle 
                                    Halfedge_handle;
  typedef typename Decorator::Halffacet_handle 
                                    Halffacet_handle;
  typedef typename Decorator::SHalfedge_handle 
                                    SHalfedge_handle;
  typedef typename Decorator::SHalfloop_handle 
                                    SHalfloop_handle;
  typedef typename Decorator::SFace_handle 
                                    SFace_handle;
  Mark m;

public:
  Volume_setter(Mark m_in = true) : m(m_in) {}
  
  void visit(Vertex_handle ) {}
  void visit(Halfedge_handle ) {}
  void visit(Halffacet_handle ) {}
  void visit(SHalfedge_handle ) {}
  void visit(SHalfloop_handle ) {}
  void visit(SFace_handle sf) {sf->mark() = m;}
};    

template<typename Nef_3>
class Mark_bounded_volumes : public Modifier_base<typename Nef_3::SNC_structure>
{  typedef typename Nef_3::SNC_structure         SNC_structure;
   typedef typename SNC_structure::SNC_decorator SNC_decorator;
   typedef typename SNC_structure::Infi_box      Infi_box;
   typedef typename Nef_3::SFace_handle          SFace_handle;
   typedef typename Nef_3::Volume_iterator       Volume_iterator;
   typedef typename Nef_3::Shell_entry_iterator 
                           Shell_entry_iterator;
   typedef typename Nef_3::Mark                  Mark;

   Mark flag;

public:
   Mark_bounded_volumes (Mark b=true) : flag(b) {}

   void operator()(SNC_structure &snc)
   {  // mark bounded volumes
      Volume_iterator vol_it = snc.volumes_begin();
      CGAL_assertion ( vol_it != snc.volumes_end() );
      if ( Infi_box::extended_kernel() ) ++vol_it; // skip Infi_box
      CGAL_assertion ( vol_it != snc.volumes_end() );
      ++vol_it; // skip unbounded volume
      Volume_setter<SNC_structure,Mark> vs(flag);
      SNC_decorator D(snc);
      for (; vol_it != snc.volumes_end(); ++vol_it)
      {  vol_it->mark() = flag; // mark
	Shell_entry_iterator it;
	CGAL_forall_shells_of(it,vol_it)
	D.visit_shell_objects(SFace_handle(it),vs);
      }
      // simplify(); // in Nef_3.delegate()
   }
};

} //namespace CGAL
#endif // CGAL_NEF3_MARK_BOUNDED_VOLUMES_H

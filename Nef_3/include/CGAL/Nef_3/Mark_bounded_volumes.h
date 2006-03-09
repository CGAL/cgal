// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Ralf Osbild <osbild@mpi-sb.mpg.de>

#ifndef CGAL_NEF3_MARK_BOUNDED_VOLUMES_H
#define CGAL_NEF3_MARK_BOUNDED_VOLUMES_H

#include <CGAL/Nef_polyhedron_3.h>

CGAL_BEGIN_NAMESPACE

template<typename Const_decorator, typename Mark>
class Volume_setter {

  typedef typename Const_decorator::Vertex_const_handle 
                                    Vertex_const_handle;
  typedef typename Const_decorator::Halfedge_const_handle 
                                    Halfedge_const_handle;
  typedef typename Const_decorator::Halffacet_const_handle 
                                    Halffacet_const_handle;
  typedef typename Const_decorator::SHalfedge_const_handle 
                                    SHalfedge_const_handle;
  typedef typename Const_decorator::SHalfloop_const_handle 
                                    SHalfloop_const_handle;
  typedef typename Const_decorator::SFace_const_handle 
                                    SFace_const_handle;
  Mark m;

public:
  Volume_setter(Mark m_in = true) : m(m_in) {}
  
  void visit(Vertex_const_handle v) {}
  void visit(Halfedge_const_handle e) {}
  void visit(Halffacet_const_handle f) {}
  void visit(SHalfedge_const_handle se) {}
  void visit(SHalfloop_const_handle sl) {}
  void visit(SFace_const_handle sf) {sf->mark() = m;}
};    

template<typename Nef_3>
class Mark_bounded_volumes : public Modifier_base<typename Nef_3::SNC_structure>
{  typedef typename Nef_3::SNC_structure        SNC_structure;
   typedef typename SNC_structure::Infi_box     Infi_box;
   typedef typename Nef_3::Volume_iterator      Volume_iterator;
   typedef typename Nef_3::Shell_entry_const_iterator 
                           Shell_entry_const_iterator;
   typedef typename Nef_3::Mark                 Mark;

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
      Volume_setter<Nef_3,Mark> vs(flag);
      Nef_3 cd(snc);
      for (; vol_it != snc.volumes_end(); ++vol_it)
      {  vol_it->mark() = flag; // mark
	Shell_entry_const_iterator it;
	CGAL_forall_shells_of(it,vol_it)
	cd.visit_shell_objects(SFace_const_handle(it),vs);
      }
      // simplify(); // in Nef_3.delegate()
   }
};

CGAL_END_NAMESPACE
#endif // CGAL_NEF3_MARK_BOUNDED_VOLUMES_H

// Copyright (c) 2005  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

#ifndef CGAL_NEF3_RELABEL_VOLUME_H
#define CGAL_NEF3_RELABEL_VOLUME_H

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/Mark_bounded_volumes.h>

namespace CGAL {

/*
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

  void visit(Vertex_handle v) {}
  void visit(Halfedge_handle e) {}
  void visit(Halffacet_handle f) {}
  void visit(SHalfedge_handle se) {}
  void visit(SHalfloop_handle sl) {}
  void visit(SFace_handle sf) {sf->mark() = m;}
};
*/

template<typename Nef_3>
class Relabel_volume : public Modifier_base<typename Nef_3::SNC_structure> {

  typedef typename Nef_3::SNC_structure                SNC_structure;
  typedef typename SNC_structure::SNC_decorator        SNC_decorator;
  typedef typename Nef_3::SFace_handle                 SFace_handle;
  typedef typename Nef_3::Volume_iterator              Volume_iterator;
  typedef typename Nef_3::Shell_entry_iterator         Shell_entry_iterator;
  typedef typename Nef_3::Mark                         Mark;

  Mark flag;

 public:
  Relabel_volume (Mark b=true) : flag(b) {}

    void operator()(SNC_structure &snc) {

      Volume_iterator c = snc.volumes_begin();
      c->mark() = flag;
      SNC_decorator D(snc);
      Volume_setter<SNC_structure,Mark> vs(flag);
      Shell_entry_iterator it;
      CGAL_forall_shells_of(it,c)
        D.visit_shell_objects(SFace_handle(it),vs);
    }
};

} //namespace CGAL
#endif // CGAL_NEF3_RELABEL_VOLUME_H

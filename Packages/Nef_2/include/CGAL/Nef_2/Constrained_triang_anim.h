// Copyright (c) 1997-2000  Max-Planck-Institute Saarbrucken (Germany).
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_PM_CONSTR_TRIANG_ANIM_H
#define CGAL_PM_CONSTR_TRIANG_ANIM_H

#include <CGAL/Nef_2/PM_visualizor.h>

CGAL_BEGIN_NAMESPACE

template <class GT>
class Constrained_triang_anim {

  CGAL::Window_stream _W;
public:
  typedef CGAL::Window_stream   VDEVICE;
  typedef typename GT::GEOMETRY GEOM;
  typedef typename GT::Base     PMDEC;
  typedef typename PMDEC::Point Point;

  Constrained_triang_anim() : _W(400,400) 
  { _W.set_show_coordinates(true); _W.init(-120,120,-120,5); _W.display(); }
  VDEVICE& device() { return _W; } 

void post_init_animation(GT& gpst)
{ 
  PM_visualizor<PMDEC,GEOM> V(_W,gpst);
  V.point(V.target(gpst.e_search)) = Point(-120,0);
  // to draw we have to embed the virtual search vertex
  V.draw_skeleton(CGAL::BLUE);
  _W.read_mouse();
}

void pre_event_animation(GT& gpst)
{ }

void post_event_animation(GT& gpst)
{ PM_visualizor<PMDEC,GEOM> V(_W,gpst);
  V.draw_ending_bundle(gpst.event,CGAL::GREEN);
  _W.read_mouse();
}

void post_completion_animation(GT& gpst)
{ _W.clear();
  PM_visualizor<PMDEC,GEOM> V(_W,gpst);
  V.draw_skeleton(CGAL::BLACK);
  _W.read_mouse(); }

};

CGAL_END_NAMESPACE
#endif // CGAL_PM_CONSTR_TRIANG_ANIM_H


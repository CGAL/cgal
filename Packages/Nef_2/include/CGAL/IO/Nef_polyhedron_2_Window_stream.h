// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/IO/Nef_polyhedron_2_Window_stream.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/Nef_polyhedron_2.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Nef polyhedra in the plane
// ============================================================================

#ifndef NEF_POLYHEDRON_2_WINDOW_STREAM_H
#define NEF_POLYHEDRON_2_WINDOW_STREAM_H

#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/Nef_2/PM_visualizor.h>

CGAL_BEGIN_NAMESPACE

static int frame_default = 100;
static bool show_triangulation = false;

template <typename T>
CGAL::Window_stream& operator<<(CGAL::Window_stream& ws, 
const Nef_polyhedron_2<T>& P)
{
  typedef Nef_polyhedron_2<T> Polyhedron;
  typedef typename T::RT RT;
  typedef typename T::Standard_RT Standard_RT;
  typedef typename Polyhedron::Topological_explorer TExplorer;
  typedef typename Polyhedron::Point            Point;
  typedef typename Polyhedron::Line             Line;
  typedef CGAL::PM_BooleColor<TExplorer> BooleColor;
  typedef CGAL::PM_visualizor<TExplorer,T,BooleColor> Visualizor;

  TExplorer D = P.explorer();
  const T& E = Nef_polyhedron_2<T>::EK;

  Standard_RT frame_radius = frame_default;
  E.determine_frame_radius(D.points_begin(),D.points_end(),frame_radius);
  RT::set_R(frame_radius);
  Visualizor PMV(ws,D); PMV.draw_map();
  if (show_triangulation) {
    P.init_locator();
    Visualizor V(ws,P.locator().triangulation());
    V.draw_skeleton(CGAL::BLUE);
  }  

  return ws;
}


CGAL_END_NAMESPACE

#endif // NEF_POLYHEDRON_2_WINDOW_STREAM_H



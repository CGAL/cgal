// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $$
// release_date  : $$
//
// file          : include/CGAL/IO/Qt_widget_Leda_rat.h
// package       : Arrangement
// maintainer    : Efi Fogel <efif@post.tau.ac.il>
// author(s)     : Efi Fogel <efif@post.tau.ac.il>
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================

#ifndef CGAL_QT_WIDGET_LEDA_RAT_H
#define CGAL_QT_WIDGET_LEDA_RAT_H

#include <CGAL/rat_leda_in_CGAL_2.h>
#include <CGAL/IO/Qt_widget.h>
#include <LEDA/rat_point.h>
#include <LEDA/rat_segment.h>

CGAL_BEGIN_NAMESPACE

Qt_widget & operator<<(Qt_widget & ws, const leda::rat_point & p)
{
  int x = ws.x_pixel(p.xcoordD());
  int y = ws.y_pixel(p.ycoordD());
  // ws.get_painter().drawPoint(x,y);
  ws.get_painter().setBrush(ws.get_painter().pen().color());
  ws.get_painter().drawEllipse(x-2, y-2, 4, 4);

  return ws;
}

Qt_widget & operator<<(Qt_widget & ws, const leda::rat_segment & seg)
{
  ws.get_painter().drawLine(ws.x_pixel(seg.xcoord1D()),
                            ws.y_pixel(seg.ycoord1D()),
                            ws.x_pixel(seg.xcoord2D()),
                            ws.y_pixel(seg.ycoord2D()));
  return ws;
}

CGAL_END_NAMESPACE

#endif

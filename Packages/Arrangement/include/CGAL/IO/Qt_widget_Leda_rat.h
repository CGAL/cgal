// Copyright (c) 2001  Tel-Aviv University (Israel).
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
// Author(s)     : Efi Fogel <efif@post.tau.ac.il>

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

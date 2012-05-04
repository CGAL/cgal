// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef QT_LAYERS_H
#define QT_LAYERS_H

#include <CGAL/IO/Qt_widget_layer.h>
#include "Virtual_Voronoi_diagram_2.h"

template< class T >
class Voronoi_diagram_layer : public CGAL::Qt_widget_layer
{
private:
  T* vd;

public:
  Voronoi_diagram_layer(T* vd) : vd(vd) {}

  void set(T* vvd) { vd = vvd; }

  void draw() {
    *widget << CGAL::BLUE;
#if !defined (__POWERPC__)
    *widget << CGAL::PointSize(3);
    *widget << CGAL::LineWidth(3);
#endif
    vd->draw_diagram(*widget);
  }
};

template< class T >
class Sites_layer : public CGAL::Qt_widget_layer
{
private:
  T* vd;

public:
  Sites_layer(T* vd) : vd(vd) {}

  void set(T* vvd) { vd = vvd; }

  void draw() {
    *widget << CGAL::RED;
#if !defined (__POWERPC__)
    *widget << CGAL::PointSize(6);
    *widget << CGAL::LineWidth(3);
#endif
    vd->draw_sites(*widget);
  }
};


#endif // QT_LAYERS_H

// Copyright (c) 2002  Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Radu Ursu

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Nef_2.h>

template <class Nef_polyhedron>
class Qt_layer_nef_blue : public CGAL::Qt_widget_layer
{
public:
  Qt_layer_nef_blue(Nef_polyhedron &n): Nef(n){}
  void draw()
  {
    *widget << CGAL::FillColor(CGAL::BLUE) << CGAL::GREEN;
    widget->setRasterOp(XorROP);
    *widget << Nef;
  };
private:
  Nef_polyhedron &Nef;
};//end class

template <class Nef_polyhedron>
class Qt_layer_nef_gray : public CGAL::Qt_widget_layer
{
public:

  Qt_layer_nef_gray(Nef_polyhedron &n): Nef(n){}
  void draw()
  {
    *widget << CGAL::FillColor(CGAL::GRAY) << CGAL::WHITE;
    widget->setRasterOp(XorROP);
    *widget << Nef;
  }
private:
  Nef_polyhedron &Nef;
};//end class

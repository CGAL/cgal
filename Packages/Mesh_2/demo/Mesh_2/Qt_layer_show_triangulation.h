// Copyright (c) 2003-2004 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_QT_LAYER_SHOW_TRIANGULATION_H
#define CGAL_QT_LAYER_SHOW_TRIANGULATION_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_Triangulation_2.h>

namespace CGAL {

template <class T>
class Qt_layer_show_triangulation : public Qt_widget_layer
{
public:
	
  Qt_layer_show_triangulation(T *t,
			      CGAL::Color lc = CGAL::BLUE,
			      int linewidth = 1) 
    : tr(t), color(lc), width(linewidth) {};


  void draw()
  {
    QColor old_color = widget->color();
    int old_width = widget->lineWidth();

    widget->setColor(color);
    widget->setLineWidth(width);
      
    *widget << *tr;

    widget->setLineWidth(old_width);
    widget->setColor(old_color);
  };

private:
  T *tr;
  CGAL::Color color;
  int width;
};//end class 

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_TRIANGULATION_H

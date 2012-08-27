// Copyright (c) 2003-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Laurent Rineau

#ifndef CGAL_QT_LAYER_SHOW_TRIANGULATION_H
#define CGAL_QT_LAYER_SHOW_TRIANGULATION_H

#include "Qt_widget_styled_layer.h"
#include <CGAL/IO/Qt_widget_Triangulation_2.h>

namespace CGAL {

template <class T>
class Qt_layer_show_triangulation : public Qt_widget_styled_layer
{
public:

  Qt_layer_show_triangulation(T *t,
			      CGAL::Color lc = CGAL::BLUE,
			      int linewidth = 1,
                              QObject* parent = 0, const char* name = 0)
    : Qt_widget_styled_layer(0, parent, name),
      tr(t)
  {
    color="Color";
    width="Line width";

    setColor(QColor(lc.red(), lc.green(), lc.blue()));
    setLineWidth(linewidth);
  };

  void setColor(QColor c)
  { style()->setColor(color, c); }

  void setLineWidth(int line_width)
  { style()->setInt(width, line_width); }

  void draw()
  {
    QColor old_color = widget->color();
    int old_width = widget->lineWidth();

    widget->setColor(style()->getColor(color));
    widget->setLineWidth(style()->getInt(width));

    *widget << *tr;

    widget->setLineWidth(old_width);
    widget->setColor(old_color);
  };

private:
  T *tr;
  QString color;
  QString width;
};//end class

} // namespace CGAL

#endif // CGAL_QT_LAYER_SHOW_TRIANGULATION_H

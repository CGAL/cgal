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

#include <CGAL/basic.h>


#include "Show_points.h"

namespace CGAL {

  Show_points_base::Show_points_base(Color c,
				     int pointsize,
				     PointStyle pointstyle,
				     QObject * parent,
				     const char * name)
    : Qt_widget_styled_layer(0, parent, name)
  {
    color=tr("Color");
    size=tr("Point size");
    style_name=tr("Point style");

    setColor(QColor(c.red(), c.green(), c.blue()));
    setPointSize(pointsize);
    setPointStyle(pointstyle);
  }

  Show_points_base::Show_points_base(Style* style,
				     QString points_color_name,
				     QString points_size_name,
				     QString points_style_name,
				     QObject * parent,
				     const char * name)
    : Qt_widget_styled_layer(style, parent, name),
      color(points_color_name),
      size(points_size_name),
      style_name(points_style_name)
  {}

  void Show_points_base::setColor(QColor c)
  { style()->setColor(color, c); }

  void Show_points_base::setPointSize(int pointsize)
  { style()->setInt(size, pointsize); }

  void Show_points_base::setPointStyle(PointStyle pointstyle)
  { style()->setInt(style_name, static_cast<int>(pointstyle)); }

} // namespace CGAL

// moc_source_file: Show_points.h
#include "Show_points.moc"


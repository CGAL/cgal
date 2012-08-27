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

#ifndef SHOW_SEGMENTS_BASE_H
#define SHOW_SEGMENTS_BASE_H

#include "Qt_widget_styled_layer.h"

namespace CGAL {

class Show_segments_base: public Qt_widget_styled_layer {
  Q_OBJECT
public:
  typedef Qt_widget_styled_layer::Style Style;

  Show_segments_base(Color c,
                     int linewidth,
                     QObject* parent = 0, const char* name = 0);

  Show_segments_base(Style* style,
                     QString line_color_name,
                     QString line_width_name,
                     QObject* parent = 0, const char* name = 0);

public slots:
  void setColor(QColor);
  void setLineWidth(int);

protected:
  QString color;
  QString width;
}; //end Show_segments_base

} // namespace CGAL

#endif // SHOW_SEGMENTS_BASE_H

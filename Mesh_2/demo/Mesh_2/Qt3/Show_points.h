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

#ifndef SHOW_POINTS_H
#define SHOW_POINTS_H

#include <qmap.h>
#include <qstring.h>
#include <qvariant.h>
#include <qcolor.h>

#include "Qt_widget_styled_layer.h"
#include <CGAL/function_objects.h>

namespace CGAL {

class Show_points_base: public Qt_widget_styled_layer {
  Q_OBJECT
public:
  typedef Qt_widget_styled_layer::Style Style;

  Show_points_base(Color c,
		   int pointsize,
		   PointStyle pointstyle,
		   QObject * parent=0, const char * name=0);

  Show_points_base(Style* style,
		   QString points_color_name,
		   QString points_size_name,
		   QString points_style_name,
		   QObject * parent=0, const char * name=0);

public slots:
  void setColor(QColor);
  void setPointSize(int);
  void setPointStyle(PointStyle);

protected:
  QString color;
  QString size;
  QString style_name;
}; // end Show_points_base

template <class C, class It,
  class Transform = Identity<typename It::value_type> >
class Show_points : public Show_points_base {
public:
  typedef Qt_widget_styled_layer::Style Style;

  typedef It iterator;
  typedef iterator (C::* iterator_function)() const;

  Show_points(C *container,
	      iterator_function begin,
	      iterator_function end,
	      Color c = CGAL::GREEN,
	      int pointsize = 3,
	      PointStyle pointstyle = CGAL::DISC,
	      QObject * parent=0, const char * name=0)
    : Show_points_base(c, pointsize, pointstyle,
		       parent, name),
      cont(container), _begin(begin), _end(end) {};

  Show_points(C *container,
	      iterator_function begin,
	      iterator_function end,
	      Style* style,
	      QString points_color_name,
	      QString points_size_name,
	      QString points_style_name,
	      QObject * parent=0, const char * name=0)
    : Show_points_base(style,
		       points_color_name,
		       points_size_name,
		       points_style_name,
		       parent, name),
      cont(container), _begin(begin), _end(end) {};

  void draw()
  {
    widget->lock();
    {
      QColor old_color = widget->color();
      int old_size = widget->pointSize();
      PointStyle old_point_style = widget->pointStyle();

      widget->setColor(style()->getColor(color));
      widget->setPointSize(style()->getInt(size));
      widget->setPointStyle(static_cast<CGAL::PointStyle>(style()->
							  getInt(style_name)));

      for(iterator it = (cont->*_begin)();
	  it!=(cont->*_end)();
	  ++it)
	*widget << Transform()(*it);

      widget->setColor(old_color);
      widget->setPointSize(old_size);
      widget->setPointStyle(old_point_style);
    }
    widget->unlock();
  };

private:
  C	*cont;
  iterator_function _begin;
  iterator_function _end;
};//end class

} // namespace CGAL

#endif // SHOW_POINTS_H

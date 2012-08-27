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

#ifndef SHOW_SEGMENTS_H
#define SHOW_SEGMENTS_H

#include "Show_segments_base.h"

namespace CGAL {

template <class C, class It,
	  class Transform = Identity<typename It::value_type> >
class Show_segments : public Show_segments_base {
public:
  typedef Qt_widget_styled_layer::Style Style;

  typedef It iterator;
  typedef iterator (C::* iterator_function)() const;

  Show_segments(C *container,
                iterator_function begin,
                iterator_function end,
                Color c=CGAL::GREEN,
                int linewidth=3,
                QObject* parent = 0, const char* name = 0)
    : Show_segments_base(c, linewidth,
		      parent, name),
      cont(container), _begin(begin), _end(end) {};

  Show_segments(C *container,
                iterator_function begin,
                iterator_function end,
                Style* style,
                QString line_color_name,
                QString line_width_name,
                QObject* parent = 0, const char* name = 0)
    : Show_segments_base(style, line_color_name, line_width_name,
		      parent, name),
      cont(container), _begin(begin), _end(end) {};

  void set_container(C* container)
  {
    cont = container;
  }

  void draw()
  {
    if( cont != 0 )
      {
        widget->lock();
        {
          QColor old_color = widget->color();
          int old_width = widget->lineWidth();

          widget->setColor(style()->getColor(color));
          widget->setLineWidth(style()->getInt(width));

          for(iterator it = (cont->*_begin)();
              it!=(cont->*_end)();
              ++it)
            {
              *widget << Transform()(*it);
            }

          widget->setColor(old_color);
          widget->setLineWidth(old_width);
        }
        widget->unlock();
      }
  };

private:
  C	*cont;
  iterator_function _begin;
  iterator_function _end;
};//end class

} // namespace CGAL

#endif // SHOW_SEGMENTS_H

// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL: svn+ssh://guyzucke@scm.gforge.inria.fr/svn/cgal/trunk/Boolean_set_operations_2/demo/Boolean_set_operations_2/Qt_widget_locate_layer.h $
// $Id: Qt_widget_locate_layer.h 37003 2007-03-10 16:55:12Z spion $
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#include "Qt_widget_locate_layer.h"
#include "boolean_operations_2.h"

void Qt_widget_locate_layer::draw()
{
	widget->lock();
	const Polygon_2& outer_boundary = m_pgn.outer_boundary();
	*widget << CGAL::YELLOW;
	if(m_found_pgn)
   {
     if(outer_boundary.is_empty())
     {
       // no boundary -> unbounded polygon
       Iso_rectangle rect(Point_2(widget->x_min(), widget->y_min()),
                         Point_2(widget->x_max(), widget->y_max()));
       *widget << rect;
     }
     else
       *widget << outer_boundary;
     for(Hole_const_iterator hit = m_pgn.holes_begin();
         hit != m_pgn.holes_end();
         ++hit)
     {
       *widget << *hit;
     }
   }
   widget->unlock();
} //end draw 

void Qt_widget_locate_layer::mousePressEvent(QMouseEvent *e)
{
	if(e->button() == Qt::LeftButton)
   {
      Coord_type x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      typedef  Traits::Point_2                   Arc_point_2;
      Arc_point_2 query_pt(x, y);
      if(m_window->red_active)
        m_found_pgn = m_window->red_set.locate(query_pt, m_pgn);
      else
        m_found_pgn = m_window->blue_set.locate(query_pt, m_pgn);
      widget->redraw();
    }
} //end mousePressEvent

void Qt_widget_locate_layer::activating()
{
    m_oldcursor = widget->cursor();
    widget->setCursor(m_cursor);
    m_oldpolicy = widget->focusPolicy();
    widget->setFocusPolicy(QWidget::StrongFocus);
}

void Qt_widget_locate_layer::deactivating()
{
    reset();
    widget->setCursor(m_oldcursor);
    widget->setFocusPolicy(m_oldpolicy);
    widget->redraw();
}

void Qt_widget_locate_layer::reset()
{
    m_found_pgn = false;
    m_pgn.clear();
}

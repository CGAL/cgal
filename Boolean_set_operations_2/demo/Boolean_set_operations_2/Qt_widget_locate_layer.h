// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// $URL$
// $Id$
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_QT_WIDGET_LOCATE_LAYER_H
#define CGAL_QT_WIDGET_LOCATE_LAYER_H

#include <CGAL/IO/Qt_widget_layer.h>
#include <qcursor.h>
#include "Qt_widget_circ_polygon.h"

extern bool                                      red_active;
extern Polygon_set                               red_set;
extern Polygon_set                               blue_set;


class  Qt_widget_locate_layer : public CGAL::Qt_widget_layer
{

    //Data members
    Polygon_with_holes     m_pgn;

    QWidget::FocusPolicy  m_oldpolicy;
    QCursor               m_oldcursor;
    QCursor               m_cursor;

    bool m_found_pgn;


  public:

    Qt_widget_locate_layer(const QCursor c=QCursor(Qt::crossCursor),
                           QObject* parent = 0,
                           const char* name = 0)
      : CGAL::Qt_widget_layer(parent, name),
        m_cursor(c),
        m_found_pgn(false)
    {}

    void draw()
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
    };


  protected:

  void mousePressEvent(QMouseEvent *e)
  {
    if(e->button() == Qt::LeftButton)
    {
      Coord_type x, y;
      widget->x_real(e->x(), x);
      widget->y_real(e->y(), y);
      typedef  Traits::Point_2                   Arc_point_2;
      Arc_point_2 query_pt(x, y);
      if(red_active)
        m_found_pgn = red_set.locate(query_pt, m_pgn);
      else
        m_found_pgn = blue_set.locate(query_pt, m_pgn);
      widget->redraw();
    }

  };//end mousePressEvent

  void activating()
  {
    m_oldcursor = widget->cursor();
    widget->setCursor(m_cursor);
    m_oldpolicy = widget->focusPolicy();
    widget->setFocusPolicy(QWidget::StrongFocus);
  };

  void deactivating()
  {
    reset();
    widget->setCursor(m_oldcursor);
    widget->setFocusPolicy(m_oldpolicy);
    widget->redraw();
  };

  public:
  void reset()
  {
    m_found_pgn = false;
    m_pgn.clear();
  }
};

#endif

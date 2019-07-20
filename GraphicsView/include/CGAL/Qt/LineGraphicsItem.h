// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_LINE_GRAPHICS_ITEM_H
#define CGAL_QT_LINE_GRAPHICS_ITEM_H

#include <CGAL/license/GraphicsView.h>


#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/utility.h>

#include <QGraphicsScene>
#include <QPainter>
#include <QStyleOption>

namespace CGAL {
namespace Qt {

template <typename CK>
class LineGraphicsItem : public GraphicsItem
{
  typedef typename CK::Line_2 Line_2;

public:
  LineGraphicsItem();

  void modelChanged();

public:
  QRectF boundingRect() const;
  
  void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
  

  const QPen& Pen() const
  {
    return this->pen;
  }


  void setPen(const QPen& pen_)
  {
    this->pen = pen_;
  }

  
  void setLine(const Line_2& a);

  Line_2 line() const
  {
    return line_;
  }

protected:

  QPainter* m_painter;
  PainterOstream<CK> painterostream;


  QPen pen;

  Line_2 line_;
};


template <typename CK>
void 
LineGraphicsItem<CK>::setLine(const Line_2& a)
{
  line_ = a;
  update();
}

template <typename CK>
LineGraphicsItem<CK>::LineGraphicsItem()
  : painterostream(0)
{
  this->hide();
  setZValue(13);
}

template <typename CK>
QRectF 
LineGraphicsItem<CK>::boundingRect() const
{
  if(scene()){
    return CGAL::Qt::viewportsBbox(scene());
  }
  return QRectF();
}




template <typename CK>
void 
LineGraphicsItem<CK>::paint(QPainter *painter, 
                                    const QStyleOptionGraphicsItem * /*option*/,
                                    QWidget * /*widget*/)
{
  painter->setPen(this->Pen());
  painterostream = PainterOstream<CK>(painter, boundingRect());
  painterostream << line_;
}




template <typename CK>
void 
LineGraphicsItem<CK>::modelChanged()
{
  update();
}


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_LINE_GRAPHICS_ITEM_H

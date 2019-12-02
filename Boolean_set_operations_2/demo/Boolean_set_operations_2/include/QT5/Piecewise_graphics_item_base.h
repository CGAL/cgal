// Copyright (c) 2012,2018  Tel-Aviv University (Israel).
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
// Author(s) : Apurva Bhatt <response2apurva@gmail.com>
//             Ronnie Gandhi <ronniegandhi19999@gmail.com>
//             Efi Fogel <efifogel@gmain.com>

#ifndef CGAL_QT_PIECEWISE_GRAPHICS_ITEM_BASE_H
#define CGAL_QT_PIECEWISE_GRAPHICS_ITEM_BASE_H

#include <boost/optional.hpp>
#include <boost/utility.hpp>

#include <CGAL/function_objects.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <CGAL/Qt/Converter.h>

#include <QPainter>
#include <QBrush>
#include <QPen>

#include "Typedefs.h"

//This class contains all the necessary drawing tools needed by the demo
namespace CGAL {
namespace Qt {

class Piecewise_graphics_item_base : public GraphicsItem {
protected:
  //constructor
  Piecewise_graphics_item_base() {}

public:
  void updateBoundingBox();

  //updating the box
  void modelChanged()
  {
    updateBoundingBox();
    //updates the widget
    update();
  }

  QRectF boundingRect() const { return m_bounding_rect; }

  void paint(QPainter* aPainter, const QStyleOptionGraphicsItem* aOption,
             QWidget* aWidget);

  const QBrush& brush() const { return m_brush; }

  void setBrush(const QBrush& aBrush) { m_brush = aBrush; }

  const QPen& pen() const{ return m_pen; }

  void setPen(const QPen& aPen) { m_pen = aPen; }

protected:
  //a converter
  typedef Converter<Kernel> ToQtConverter;

  //for adding 2 bbox and initializing if its null
  struct Bbox_builder {
    void add (Bbox_2 const& aBbox)
    {
      if (bbox) bbox = *bbox + aBbox;
      else bbox = aBbox;
    }
    boost::optional<Bbox_2> bbox;
  };

  virtual bool isModelEmpty() const = 0;

  virtual void draw_model (QPainterPath& aPath) = 0;

  virtual void update_bbox(Bbox_builder& aBBoxBuilder) = 0;

protected:
  //qt5 drawing tools
  QRectF m_bounding_rect;
  QBrush m_brush;
  QPen m_pen;
};

//
void Piecewise_graphics_item_base::
paint(QPainter* aPainter,
      const QStyleOptionGraphicsItem* /* aOption */,
      QWidget* /* aWidget */)
{
  //if there is any data to draw
  if (!isModelEmpty()) {
    QPainterPath l_path;
    draw_model(l_path);

    //setting drawing tools
    aPainter->setPen(m_pen);
    aPainter->setBrush(m_brush);
    //drawing l_path
    aPainter->drawPath(l_path);
  }
}

// to let the bounding box only grow, so that when vertices get removed
// the maximal bbox gets refreshed in the GraphicsView
void Piecewise_graphics_item_base::updateBoundingBox()
{
  //if there is any data to draw
  if (!isModelEmpty()) {
    //"Prepares the item for a geometry change
    prepareGeometryChange();
    //update();

    Bbox_builder l_bbox_builder;

    update_bbox(l_bbox_builder);

    //if bbox exits convert it to qt applicable
    if (l_bbox_builder.bbox) {
      ToQtConverter to_Qt;
      m_bounding_rect = to_Qt(*l_bbox_builder.bbox);
    }
  }
}

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PIECEWISE_GRAPHICS_ITEM_BASE_H

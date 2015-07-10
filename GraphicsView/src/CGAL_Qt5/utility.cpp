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
// 
//
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#include <CGAL/Qt/utility.h>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QList>
#include <QPoint>
#include <QPointF>

namespace CGAL {
namespace Qt {

QRectF mapToScene(const QGraphicsView* v, const QRect rect)
{
  QPointF top_left = v->mapToScene(rect.topLeft());
  QPointF size = v->mapToScene(rect.bottomRight());
  size -= top_left;
  return QRectF(top_left.x(),
		top_left.y(),
		size.x(),
		size.y());
}

QRectF viewportsBbox(const QGraphicsScene* scene) {
   QRectF rect;
   Q_FOREACH(QGraphicsView* view, scene->views())
   {
     rect |= mapToScene(view, view->viewport()->rect());
   }
   rect = rect.normalized();
   return rect;
}

} // namespace Qt
} // namespace CGAL

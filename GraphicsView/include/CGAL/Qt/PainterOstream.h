// Copyright (c) 2008  GeometryFactory Sarl (France).
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
// Author(s)     : Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Laurent Rineau <Laurent.Rineau@geometryfactory.com>

#ifndef CGAL_QT_PAINTER_OSTREAM_H
#define CGAL_QT_PAINTER_OSTREAM_H

#include <QPainter>
#include <QPen>
#include <QRectF>
#include <CGAL/Qt/Converter.h>

namespace CGAL {
namespace Qt {

template <typename K>
class PainterOstream {

private:
  QPainter* qp;
  Converter<K> convert;
  
public:
  PainterOstream(QPainter* p, QRectF rect = QRectF())
    : qp(p), convert(rect)
  {}

  PainterOstream& operator<<(const Point_2<K>& p)
  {
    qp->drawPoint(convert(p));
    return *this;
  }
  
  PainterOstream& operator<<(const Segment_2<K>& s)
  {
    qp->drawLine(convert(s.source()), convert(s.target()));
    return *this;
  }
  
  
  PainterOstream& operator<<(const Ray_2<K>& r)
  {
    qp->drawLine(convert(r));
    return *this;
  }

  
  PainterOstream& operator<<(const Line_2<K>& l)
  {
    qp->drawLine(convert(l));
    return *this;
  }


  PainterOstream& operator<<(const Triangle_2<K>& t)
  {
    qp->drawPolygon(convert(t));
    return *this;
  }

  PainterOstream& operator<<(const Iso_rectangle_2<K>& r)
  {
    qp->drawRect(convert(r));
    return *this;
  }

  PainterOstream& operator<<(const Circle_2<K>& c)
  {
    qp->drawRect(convert(c.bbox()));
    return *this;
  }
};

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_PAINTER_OSTREAM_H

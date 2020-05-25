// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef CGAL_QT_CURVE_GRAPHICS_ITEM_H
#define CGAL_QT_CURVE_GRAPHICS_ITEM_H

#include "ArrangementPainterOstream.h"
#include "Utils.h"

#include <CGAL/Qt/Converter.h>
#include <CGAL/Qt/GraphicsItem.h>
#include <QGraphicsScene>

namespace CGAL {
namespace Qt {

/**
   Draws selected curves and vertices of an arrangement.
*/
template < class ArrTraits >
class CurveGraphicsItem : public GraphicsItem, public QGraphicsSceneMixin
{
public:
  // known curve types
  typedef ArrTraits Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
  typedef typename Traits::Curve_2 Curve_2;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Kernel::Point_2 Kernel_point_2;

public: // ctors
  CurveGraphicsItem( );

public: // methods
  void paint(QPainter* painter,
                     const QStyleOptionGraphicsItem* /* option */,
                     QWidget* /* widget */) override;
  QRectF boundingRect( ) const override;
  void insert( const X_monotone_curve_2& curve );
  void insert( const Point_2& point );
  void clear( );

public Q_SLOTS:
  void modelChanged( );
  const QColor& edgeColor( ) const;
  void setEdgeColor( const QColor& color );
  int edgeWidth( ) const;
  void setEdgeWidth( int width );
  const QColor& vertexColor( ) const;
  void setVertexColor( const QColor& color );
  int vertexRadius( ) const;
  void setVertexRadius( int radius );

protected: // methods
  void updateBoundingBox( );

protected: // fields
  CGAL::Qt::Converter< Kernel > convert;
  ArrangementPainterOstream< Traits > painterOstream;
  std::vector< X_monotone_curve_2 > curves;
  std::vector< Point_2 > points;
  CGAL::Bbox_2 boundingBox;

  QColor m_edgeColor;
  int m_edgeWidth;
  QColor m_vertexColor;
  int m_vertexRadius;

}; // class CurveGraphicsItem

} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_CURVE_GRAPHICS_ITEM_H

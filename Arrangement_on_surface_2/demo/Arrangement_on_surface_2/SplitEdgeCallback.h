// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>

#ifndef SPLIT_EDGE_CALLBACK_H
#define SPLIT_EDGE_CALLBACK_H

#include "Callback.h"
#include <QColor>

class QGraphicsSceneMouseEvent;
class QGraphicsScene;
class QGraphicsLineItem;
namespace CGAL
{
template <typename Coefficient_>
class Arr_algebraic_segment_traits_2;
}

class SplitEdgeCallbackBase : public CGAL::Qt::Callback
{
public:
  virtual void setColor( QColor c );
  QColor getColor( ) const;

protected:
  SplitEdgeCallbackBase( QObject* parent );
  QColor color;
}; // SplitEdgeCallbackBase

/**
   Handles splitting of arrangement curves selected from the scene.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template <typename Arr_>
class SplitEdgeCallback : public SplitEdgeCallbackBase
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement::Halfedge_iterator Halfedge_iterator;
  typedef typename Arrangement::Vertex_iterator Vertex_iterator;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Arrangement::Curve_handle Curve_handle;
  typedef typename Arrangement::Originating_curve_iterator
    Originating_curve_iterator;
  typedef typename Arrangement::Induced_edge_iterator Induced_edge_iterator;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Intersect_2 Intersect_2;
  typedef typename Traits::Equal_2 Equal_2;
  typedef typename Traits::Multiplicity Multiplicity;
  typedef typename Traits::Point_2 Point_2;

  SplitEdgeCallback( Arrangement* arr_, QObject* parent );
  void setScene( QGraphicsScene* scene_ );
  virtual void setColor( QColor c );
  void reset( );

  void slotModelChanged( );

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );

  Point_2 snapPoint( QGraphicsSceneMouseEvent *event );

  template < class TTraits >
  void splitEdges( const Point_2& pt, TTraits traits );

  template < class Coefficient_ >
  void splitEdges( const Point_2& pt,
                   CGAL::Arr_algebraic_segment_traits_2<Coefficient_> traits);

  template < class TTraits >
  void updateGuide( const Point_2& pt, TTraits traits );

  Traits traits;
  Arrangement* arr;
  bool hasFirstPoint;
  Point_2 p1;
  Point_2 p2;
  Intersect_2 intersectCurves;
  Equal_2 areEqual;
  QGraphicsLineItem* segmentGuide;
}; // class SplitEdgeCallback

#endif // SPLIT_EDGE_CALLBACK_H

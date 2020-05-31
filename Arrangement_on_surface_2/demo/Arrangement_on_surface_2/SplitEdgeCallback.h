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
#include <QEvent>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Arrangement_with_history_2.h>

#include "CurveGraphicsItem.h"
#include "Utils.h"
#include "ISnappable.h"

class SplitEdgeCallbackBase : public CGAL::Qt::Callback, public ISnappable
{
public:
  void setSnappingEnabled( bool b );
  void setSnapToGridEnabled( bool b );

  virtual void setColor( QColor c );
  QColor getColor( ) const;

protected:
  SplitEdgeCallbackBase( QObject* parent );

  bool snappingEnabled;
  bool snapToGridEnabled;
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
  typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
  typedef typename Traits::Intersect_2 Intersect_2;
  typedef typename Traits::Equal_2 Equal_2;
  typedef typename Traits::Multiplicity Multiplicity;
  typedef typename ArrTraitsAdaptor< Traits >::Point_2 Point_2;
  typedef typename ArrTraitsAdaptor< Traits >::CoordinateType CoordinateType;
  typedef typename Kernel::Segment_2 Segment_2;

  SplitEdgeCallback( Arrangement* arr_, QObject* parent );
  void setScene( QGraphicsScene* scene_ );
  virtual void setColor( QColor c );
  void reset( );

  void slotModelChanged( );

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );

  virtual Point_2 snapPoint( QGraphicsSceneMouseEvent *event );
  template < class TTraits >
  Point_2 snapPoint( QGraphicsSceneMouseEvent *event, TTraits traits );

  template < class TTraits >
  void splitEdges( const Point_2& pt, TTraits traits );
  template < typename RatKernel, typename AlgKernel, typename NtTraits >
  void splitEdges( const Point_2& pt,
                   CGAL::Arr_Bezier_curve_traits_2< RatKernel, AlgKernel, NtTraits> traits);

  template < class Coefficient_ >
  void splitEdges( const Point_2& pt,
                   CGAL::Arr_algebraic_segment_traits_2<Coefficient_> traits);

  template < class TTraits >
  void updateGuide( const Point_2& pt, TTraits traits );

  Traits traits;
  CGAL::Qt::Converter<Kernel> convert;
  Arrangement* arr;
  bool hasFirstPoint;
  Point_2 p1;
  Point_2 p2;
  Intersect_2 intersectCurves;
  Equal_2 areEqual;
  QGraphicsLineItem* segmentGuide;
  SnapToArrangementVertexStrategy<Arrangement> snapToVertexStrategy;
  SnapToGridStrategy< typename Arrangement::Geometry_traits_2 >
    snapToGridStrategy;
}; // class SplitEdgeCallback

#endif // SPLIT_EDGE_CALLBACK_H

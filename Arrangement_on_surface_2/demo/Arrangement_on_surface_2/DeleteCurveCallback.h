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

#ifndef DELETE_CURVE_CALLBACK_H
#define DELETE_CURVE_CALLBACK_H


#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

#include <CGAL/Qt/Converter.h>

#include "Callback.h"
#include "CurveGraphicsItem.h"

/**
   Handles deletion of arrangement curves selected from the scene.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < typename Arr_ >
class DeleteCurveCallback : public CGAL::Qt::Callback
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_const_handle   Halfedge_const_handle;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Halfedge_iterator       Halfedge_iterator;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef typename Arrangement::Curve_handle            Curve_handle;
  typedef typename Arrangement::Originating_curve_iterator
  Originating_curve_iterator;
  typedef typename Arrangement::Induced_edge_iterator   Induced_edge_iterator;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Point;
  typedef typename Kernel::Segment_2                    Segment;

  DeleteCurveCallback( Arrangement* arr_, QObject* parent_ );
  void setScene( QGraphicsScene* scene_ );
  QGraphicsScene* getScene( ) const;
  void reset( );
  virtual void partialReset();
  virtual std::string toString()
  {
    return ( this->deleteOriginatingCurve ) ?  "Delete Curve" : "Delete Edge";
  }

  virtual void changeDeleteMode()
  {
    this->deleteOriginatingCurve = !this->deleteOriginatingCurve;
  }

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
  void highlightNearestCurve( QGraphicsSceneMouseEvent *event );

  Compute_squared_distance_2< Traits > squaredDistance;
  CGAL::Qt::Converter< Kernel > convert;
  QGraphicsScene* scene;
  CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurve;
  Arrangement* arr;
  Halfedge_handle removableHalfedge;
  bool deleteOriginatingCurve;
}; // class DeleteCurveCallback

#endif // DELETE_CURVE_CALLBACK_H

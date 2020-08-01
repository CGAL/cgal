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

#ifndef MERGE_EDGE_CALLBACK_H
#define MERGE_EDGE_CALLBACK_H

#include "Callback.h"
#include <CGAL/Qt/Converter.h>

namespace CGAL
{
namespace Qt
{
template <typename T>
class CurveGraphicsItem;
} // namespace Qt
} // namespace CGAL

class QGraphicsScene;

/**
   Handles merging of arrangement curves selected from the scene.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < typename Arr_ >
class MergeEdgeCallback : public CGAL::Qt::Callback
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_handle         Halfedge_handle;
  typedef typename Arrangement::Halfedge_iterator       Halfedge_iterator;
  typedef typename Arrangement::Vertex_iterator         Vertex_iterator;
  typedef typename Arrangement::Geometry_traits_2       Traits;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

  MergeEdgeCallback( Arrangement* arr_, QObject* parent_ );
  void setScene( QGraphicsScene* scene_ ) override;
  void reset( );

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent* event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent* event );
  Halfedge_handle getNearestMergeableCurve( QGraphicsSceneMouseEvent* event );
  Halfedge_handle getNearestMergeableCurve( Halfedge_handle h,
                                            QGraphicsSceneMouseEvent* event );

  CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurve;
  CGAL::Qt::CurveGraphicsItem< Traits >* highlightedCurve2;
  Arrangement* arr;
  Halfedge_handle mergeableHalfedge;
  bool isFirst;
}; // class MergeEdgeCallback

#endif // MERGE_EDGE_CALLBACK_H

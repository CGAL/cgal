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

#include "GraphicsViewCurveInput.h"
#include "ForwardDeclarations.h"

class QGraphicsSceneMouseEvent;
class QGraphicsScene;
class QGraphicsLineItem;
class PointSnapperBase;
class QColor;

class SplitEdgeCallbackBase :
    public CGAL::Qt::Callback,
    public CGAL::Qt::CurveInputMethodCallback
{
public:
  void setColor( QColor c );
  QColor getColor( ) const;
  void setPointSnapper(PointSnapperBase*);
  bool eventFilter(QObject* object, QEvent* event) override;

protected:
  SplitEdgeCallbackBase( QObject* parent );
  CGAL::Qt::SegmentInputMethod segmentInputMethod;
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
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Intersect_2 Intersect_2;
  typedef typename Traits::Equal_2 Equal_2;
  typedef typename Traits::Multiplicity Multiplicity;
  typedef typename Traits::Point_2 Point_2;
  typedef typename CGAL::Qt::CurveInputMethod::Point_2 Input_point_2;

  SplitEdgeCallback( Arrangement* arr_, QObject* parent );
  void setScene(QGraphicsScene* scene_) override;
  void reset() override;

protected:
  void curveInputDoneEvent(
    const std::vector<Input_point_2>& clickedPoints,
    CGAL::Qt::CurveType type) override;

  template <class TTraits>
  void
  splitEdges(const Input_point_2& p1, const Input_point_2& p2, const TTraits*);

  template <class Coefficient_>
  void splitEdges(
    const Input_point_2& p1, const Input_point_2& p2,
    const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>*);

  Arrangement* arr;
  Intersect_2 intersectCurves;
  Equal_2 areEqual;
}; // class SplitEdgeCallback

#endif // SPLIT_EDGE_CALLBACK_H

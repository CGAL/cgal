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

#ifndef FILL_FACE_CALLBACK_H
#define FILL_FACE_CALLBACK_H

#include "Callback.h"
#include <CGAL/Qt/Converter.h>
#include <QGraphicsSceneMouseEvent>
#include <QColor>

template <typename T>
struct ArrTraitsAdaptor;

namespace CGAL
{
template <typename T>
struct Arr_trapezoid_ric_point_location;
template <typename T>
struct Arr_walk_along_line_point_location;
template <typename T>
struct Arr_simple_point_location;
template <typename T>
struct Supports_landmarks;
} // namespace CGAL

class FillFaceCallbackBase : public CGAL::Qt::Callback
{
public:
  FillFaceCallbackBase( QObject* parent );

  void setColor( QColor c );
  QColor getColor( ) const;

protected:
  QColor fillColor;                       				/*!< Qcolor object to fill a selected space */
};
/**
   Supports visualization of point location on arrangements.

   The template parameter is a CGAL::Arrangement_with_history_2 of some type.
*/
template < class Arr_ >
class FillFaceCallback : public FillFaceCallbackBase
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Halfedge_handle Halfedge_handle;
  typedef typename Arrangement::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arrangement::Halfedge_iterator Halfedge_iterator;
  typedef typename Arrangement::Face_handle Face_handle;
  typedef typename Arrangement::Face_const_handle Face_const_handle;
  typedef typename Arrangement::Vertex_const_handle Vertex_const_handle;
  typedef typename Arrangement::Halfedge_around_vertex_const_circulator
	Halfedge_around_vertex_const_circulator;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Arrangement::Curve_handle Curve_handle;
  typedef typename Arrangement::Originating_curve_iterator
	Originating_curve_iterator;
  typedef typename Arrangement::Induced_edge_iterator Induced_edge_iterator;
  typedef typename Arrangement::Ccb_halfedge_const_circulator
	Ccb_halfedge_const_circulator;
  typedef typename Arrangement::Hole_const_iterator Hole_const_iterator;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel Kernel;
  typedef typename Kernel::Point_2 Kernel_point_2;
  typedef typename Traits::Point_2 Point_2;
  typedef typename Kernel::Segment_2 Segment_2;
  typedef typename CGAL::Arr_trapezoid_ric_point_location< Arrangement >
	TrapezoidPointLocationStrategy;
  typedef typename CGAL::Arr_simple_point_location< Arrangement >
	SimplePointLocationStrategy;
  typedef typename CGAL::Arr_walk_along_line_point_location< Arrangement >
	Walk_pl_strategy;

  FillFaceCallback( Arrangement* arr_, QObject* parent_ );
  void reset( );


protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );

  void fillFace( QGraphicsSceneMouseEvent* event );

  Face_const_handle getFace( const CGAL::Object& o );
  CGAL::Object locate( const Kernel_point_2& point );
  CGAL::Object locate( const Kernel_point_2& point,
					   CGAL::Tag_false/*supportsLandmarks*/ );
  template <typename=void>
  CGAL::Object locate( const Kernel_point_2& point,
					   CGAL::Tag_true /*doesNotSupportLandmarks*/ );

  CGAL::Qt::Converter< Kernel > convert;
  CGAL::Object pointLocationStrategy;
  Arrangement* arr;
}; // class FillFaceCallback

#endif // FILL_FACE_CALLBACK_H

// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef FILL_FACE_CALLBACK_H
#define FILL_FACE_CALLBACK_H

#include "Callback.h"
#include <QColor>

class QGraphicsSceneMouseEvent;

class FillFaceCallbackBase : public CGAL::Qt::Callback
{
public:
  FillFaceCallbackBase( QObject* parent );

  void setColor( QColor c );
  QColor getColor( ) const;

protected:
  QColor fillColor;
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
  typedef typename Arrangement::Face_handle Face_handle;
  typedef typename Arrangement::Face_const_handle Face_const_handle;

  FillFaceCallback( Arrangement* arr_, QObject* parent_ );
  void reset( );

protected:
  void mousePressEvent( QGraphicsSceneMouseEvent *event );
  void mouseMoveEvent( QGraphicsSceneMouseEvent *event );
  void fillFace( QGraphicsSceneMouseEvent* event );

  Arrangement* arr;
}; // class FillFaceCallback

#endif // FILL_FACE_CALLBACK_H

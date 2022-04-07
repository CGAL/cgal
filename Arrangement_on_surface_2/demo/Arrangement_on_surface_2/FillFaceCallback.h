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
#include <CGAL/Object.h>
#include <QColor>

namespace demo_types
{
enum class TraitsType : int;
}

class QGraphicsSceneMouseEvent;

/**
   Supports visualization of point location on arrangements.
*/
class FillFaceCallbackBase : public CGAL::Qt::Callback
{
public:
  static FillFaceCallbackBase*
  create(demo_types::TraitsType, CGAL::Object arr_obj, QObject* parent);

  void setColor( QColor c );
  QColor getColor( ) const;

protected:
  FillFaceCallbackBase( QObject* parent );

  QColor fillColor;
};

#endif // FILL_FACE_CALLBACK_H

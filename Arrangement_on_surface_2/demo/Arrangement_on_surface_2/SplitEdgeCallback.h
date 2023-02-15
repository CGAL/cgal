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

#ifndef SPLIT_EDGE_CALLBACK_H
#define SPLIT_EDGE_CALLBACK_H

#include <CGAL/Object.h>
#include "Callback.h"

class QGraphicsSceneMouseEvent;
class QGraphicsScene;
class QGraphicsLineItem;
class PointSnapperBase;
class QColor;

namespace demo_types
{
enum class TraitsType : int;
}

/**
   Handles splitting of arrangement curves selected from the scene.
*/
class SplitEdgeCallbackBase : public CGAL::Qt::Callback

{
public:
  static SplitEdgeCallbackBase*
  create(demo_types::TraitsType, CGAL::Object arr_obj, QObject* parent);

  virtual void setColor(QColor c) = 0;
  virtual QColor getColor() const = 0;
  virtual void setPointSnapper(PointSnapperBase*) = 0;

protected:
  SplitEdgeCallbackBase(QObject* parent);
}; // SplitEdgeCallbackBase

#endif // SPLIT_EDGE_CALLBACK_H

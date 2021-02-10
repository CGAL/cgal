// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef CGAL_ARRANGEMENTS_DEMO_GRAPHICS_SCENE_MIXIN_H
#define CGAL_ARRANGEMENTS_DEMO_GRAPHICS_SCENE_MIXIN_H

#include <QRectF>
#include <QPoint>
#include <QPointF>

class QGraphicsScene;
class QGraphicsView;

class GraphicsSceneMixin
{
public:
  /*! Costructor */
  GraphicsSceneMixin(QGraphicsScene* scene_ = nullptr);

  /*! Destructor (virtual) */
  virtual ~GraphicsSceneMixin();
  virtual void setScene(QGraphicsScene* scene_);
  QGraphicsScene* getScene() const;
  QRectF viewportRect() const;
  QPoint fromScene(QPointF p);
  QPointF toScene(QPoint p) const;
  QGraphicsView* getView() const;

protected: // fields
  QGraphicsScene* scene;
};

#endif

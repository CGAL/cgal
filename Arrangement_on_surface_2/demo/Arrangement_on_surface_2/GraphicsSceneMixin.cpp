// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#include "GraphicsSceneMixin.h"
#include <QGraphicsScene>
#include <QGraphicsView>

GraphicsSceneMixin::GraphicsSceneMixin(QGraphicsScene* scene_) : scene{scene_}
{
}

GraphicsSceneMixin::~GraphicsSceneMixin() { }

void GraphicsSceneMixin::setScene(QGraphicsScene* scene_)
{
  this->scene = scene_;
}

QGraphicsScene* GraphicsSceneMixin::getScene() const { return this->scene; }

QRectF GraphicsSceneMixin::viewportRect() const
{
  QGraphicsView* view = this->getView();
  return view ? view->mapToScene(QRect{0, 0, view->width(), view->height()})
                  .boundingRect()
              : QRectF{};
}

QPoint GraphicsSceneMixin::fromScene(QPointF p)
{
  QGraphicsView* viewport = this->getView();
  return viewport ? viewport->mapFromScene(p) : QPoint{};
}

QPointF GraphicsSceneMixin::toScene(QPoint p) const
{
  QGraphicsView* viewport = this->getView();
  return viewport ? viewport->mapToScene(p) : QPointF{};
}

QGraphicsView* GraphicsSceneMixin::getView() const
{
  if (!this->scene) return nullptr;

  QList<QGraphicsView*> views = this->scene->views();
  if (views.size() == 0) return nullptr;
  // assumes the first view is the right one
  return views.first();
}

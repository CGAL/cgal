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

#include "FillFaceCallback.h"
#include "ArrangementTypes.h"
#include "PointLocationFunctions.h"
#include "Utils.h"

#include <QGraphicsSceneMouseEvent>


FillFaceCallbackBase::FillFaceCallbackBase(QObject* parent) :
    CGAL::Qt::Callback(parent), fillColor(::Qt::black)
{
}

void FillFaceCallbackBase::setColor(QColor c)
{
  this->fillColor = c;
  Q_EMIT modelChanged();
}

QColor FillFaceCallbackBase::getColor() const { return this->fillColor; }

template <class Arr_>
FillFaceCallback<Arr_>::FillFaceCallback(Arrangement* arr_, QObject* parent_) :
    FillFaceCallbackBase(parent_),
    arr(arr_)
{
}

template <class Arr_>
void FillFaceCallback<Arr_>::reset()
{
}

template <class Arr_>
void FillFaceCallback<Arr_>::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
  this->fillFace(event);
  Q_EMIT modelChanged();
}

template <class Arr_>
void FillFaceCallback<Arr_>::mouseMoveEvent(
  QGraphicsSceneMouseEvent* /* event */)
{
}

template <class Arr_>
void FillFaceCallback<Arr_>::fillFace(QGraphicsSceneMouseEvent* event)
{
  if (!this->fillColor.isValid())
    return;

  Face_const_handle face =
    PointLocationFunctions<Arrangement>{}.getFace(arr, event->scenePos());

  Face_handle f = this->arr->non_const_handle(face);

  if (f->color() != this->fillColor)
    f->set_color(this->fillColor);
  else
    f->set_color(::Qt::white);
}

ARRANGEMENT_DEMO_SPECIALIZE_ARR(FillFaceCallback)

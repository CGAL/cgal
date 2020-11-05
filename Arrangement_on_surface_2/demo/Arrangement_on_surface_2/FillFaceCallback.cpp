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
#include "ArrangementTypesUtils.h"
#include "Utils/PointLocationFunctions.h"
#include "Utils/Utils.h"

#include <QGraphicsSceneMouseEvent>

template <class Arr_>
class FillFaceCallback : public FillFaceCallbackBase
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Face_handle Face_handle;
  typedef typename Arrangement::Face_const_handle Face_const_handle;

  FillFaceCallback(Arrangement* arr_, QObject* parent_);
  void reset();

protected:
  void mousePressEvent(QGraphicsSceneMouseEvent* event);
  void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
  void fillFace(QGraphicsSceneMouseEvent* event);

  Arrangement* arr;
}; // class FillFaceCallback

FillFaceCallbackBase::FillFaceCallbackBase(QObject* parent) :
    CGAL::Qt::Callback(parent), fillColor(::Qt::black)
{
}

// msvc2015 doesn't play well with polymorphic lambdas
namespace
{
struct ExplicitLambda
{
  template <typename Arrangement>
  void operator()(demo_types::TypeHolder<Arrangement>)
  {
    Arrangement* arr = nullptr;
    CGAL::assign(arr, arr_obj);
    res = new FillFaceCallback<Arrangement>(arr, parent);
  }

  FillFaceCallbackBase*& res;
  CGAL::Object& arr_obj;
  QObject* parent;
};
} // anonymous namespace

FillFaceCallbackBase* FillFaceCallbackBase::create(
  demo_types::TraitsType tt, CGAL::Object arr_obj, QObject* parent)
{
  FillFaceCallbackBase* res;
  ExplicitLambda explicit_lambda{res, arr_obj, parent};
  demo_types::visitArrangementType(tt, explicit_lambda);
  return res;
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

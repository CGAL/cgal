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

#include "SplitEdgeCallback.h"
#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "CurveInputMethods.h"
#include "PointSnapper.h"
#include "Utils/ConstructSegment.h"
#include "Utils/IntersectCurves.h"
#include "Utils/SplitAndMerge.h"

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

template <typename Arr_>
class SplitEdgeCallback :
    public SplitEdgeCallbackBase,
    public CGAL::Qt::CurveInputMethodCallback
{
public:
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2 Traits;
  typedef typename Traits::X_monotone_curve_2 X_monotone_curve_2;
  typedef typename Traits::Multiplicity Multiplicity;
  typedef typename Traits::Point_2 Point_2;
  typedef typename CGAL::Qt::CurveInputMethod::Point_2 Input_point_2;

  SplitEdgeCallback(Arrangement* arr_, QObject* parent);
  void setScene(QGraphicsScene* scene_) override;
  void reset() override;

protected:
  void curveInputDoneEvent(
    const std::vector<Input_point_2>& clickedPoints,
    CGAL::Qt::CurveType type) override;

  void splitEdges(const Input_point_2& p1, const Input_point_2& p2);

  void setColor(QColor c) override;
  QColor getColor() const override;
  void setPointSnapper(PointSnapperBase*) override;
  bool eventFilter(QObject* object, QEvent* event) override;

  Arrangement* arr;
  CGAL::Qt::SegmentInputMethod segmentInputMethod;
  Intersect_curves<Traits> intersectCurves;
  Split_edge<Arrangement> edgeSplitter;
}; // class SplitEdgeCallback

SplitEdgeCallbackBase::SplitEdgeCallbackBase(QObject* parent) :
    CGAL::Qt::Callback(parent)
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
    res = new SplitEdgeCallback<Arrangement>(arr, parent);
  }

  SplitEdgeCallbackBase*& res;
  CGAL::Object& arr_obj;
  QObject* parent;
};
} // anonymous namespace

SplitEdgeCallbackBase* SplitEdgeCallbackBase::create(
  demo_types::TraitsType tt, CGAL::Object arr_obj, QObject* parent)
{
  SplitEdgeCallbackBase* res;
  ExplicitLambda explicit_lambda{res, arr_obj, parent};
  demo_types::visitArrangementType(tt, explicit_lambda);
  return res;
}

template <typename Arr_>
SplitEdgeCallback<Arr_>::SplitEdgeCallback(Arrangement* arr_, QObject* parent) :
    SplitEdgeCallbackBase(parent), arr(arr_),
    intersectCurves(this->arr->traits())
{
  segmentInputMethod.setCallback(this);
  this->setColor(::Qt::darkGray);
}

template <typename Arr_>
bool SplitEdgeCallback<Arr_>::eventFilter(QObject* object, QEvent* event)
{
  return segmentInputMethod.eventFilter(object, event);
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::setColor(QColor c)
{
  this->segmentInputMethod.setColor(c);
}

template <typename Arr_>
QColor SplitEdgeCallback<Arr_>::getColor() const
{
  return this->segmentInputMethod.getColor();
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::setPointSnapper(PointSnapperBase* snapper_)
{
  this->segmentInputMethod.setPointSnapper(snapper_);
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::setScene(QGraphicsScene* scene_)
{
  CGAL::Qt::Callback::setScene(scene_);
  this->segmentInputMethod.setScene(scene_);
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::reset()
{
  this->segmentInputMethod.reset();
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::curveInputDoneEvent(
  const std::vector<Input_point_2>& clickedPoints, CGAL::Qt::CurveType)
{
  try
  {
    this->splitEdges(clickedPoints[0], clickedPoints[1]);
  }
  catch (const std::exception& ex)
  {
    std::cerr << ex.what() << '\n';
    std::cerr << __FILE__ << ':' << __LINE__ << '\n';
    return;
  }
}

template <typename Arr_>
void SplitEdgeCallback<Arr_>::splitEdges(
  const Input_point_2& p1, const Input_point_2& p2)
{
  Construct_segment<Traits> construct_segment;

  X_monotone_curve_2 splitCurve =
    construct_segment(this->arr->traits(), p1, p2);

  for (auto hei = this->arr->halfedges_begin();
       hei != this->arr->halfedges_end(); ++hei)
  {
    X_monotone_curve_2 curve = hei->curve();
    std::vector<CGAL::Object> intersection_points;
    this->intersectCurves(splitCurve, curve, intersection_points);

    for (auto& pt_obj : intersection_points)
    {
      std::pair<Point_2, Multiplicity> pt_pair;
      if (CGAL::assign(pt_pair, pt_obj))
      {
        Point_2 splitPoint = pt_pair.first;
        this->edgeSplitter(this->arr, hei, splitPoint);
      }
    }
  }
  this->reset();
  Q_EMIT modelChanged();
}

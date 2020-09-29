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

#include "ArrangementTypes.h"
#include "ArrangementTypesUtils.h"
#include "ArrangementDemoTab.h"
#include "ArrangementGraphicsItem.h"
#include "ArrangementDemoGraphicsView.h"
#include "ArrangementCurveInputCallback.h"
#include "DeleteCurveCallback.h"
#include "PointLocationCallback.h"
#include "VerticalRayShootCallback.h"
#include "MergeEdgeCallback.h"
#include "SplitEdgeCallback.h"
#include "EnvelopeCallback.h"
#include "FillFaceCallback.h"
#include "GridGraphicsItem.h"
#include "PointSnapper.h"
#include "ConstructBoundingBox.h"

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include <QGridLayout>
#include <limits>


ArrangementDemoTabBase::ArrangementDemoTabBase( QWidget* parent ) :
  QWidget( parent ),
  GraphicsSceneMixin( new QGraphicsScene() ),
  graphicsView( new ArrangementDemoGraphicsView( this ) ),
  layout( new QGridLayout( this ) ),
  navigation( std::make_unique<CGAL::Qt::GraphicsViewNavigation>() ),
  activeCallback( nullptr ),
  arrangementGraphicsItem( nullptr ),
  gridGraphicsItem( nullptr )
{
  this->setupUi( );
}

ArrangementDemoTabBase::~ArrangementDemoTabBase( )
{
  this->unhookCallbacks();
}

void ArrangementDemoTabBase::setupUi( )
{
  auto scene = this->getScene();

  this->layout->addWidget( this->graphicsView, 0, 0 );
  this->graphicsView->setScene( scene );

  double MAX_WIDTH = (std::numeric_limits<double>::max)() / 1048576;

  double xymin = -MAX_WIDTH / 2;
  double wh = MAX_WIDTH;
  scene->setSceneRect(xymin, xymin, wh, wh);

  this->getView()->installEventFilter(this->navigation.get());
  this->getView()->viewport()->installEventFilter(this->navigation.get());
}

QGraphicsView* ArrangementDemoTabBase::getView() const
{
  return this->graphicsView;
}

CGAL::Qt::ArrangementGraphicsItemBase*
ArrangementDemoTabBase::getArrangementGraphicsItem( ) const
{
  return this->arrangementGraphicsItem;
}

GridGraphicsItem* ArrangementDemoTabBase::getGridGraphicsItem() const
{
  return this->gridGraphicsItem;
}

void ArrangementDemoTabBase::showGrid(bool val)
{
  this->gridGraphicsItem->setVisible(val);
}

void ArrangementDemoTabBase::setSnapToGrid(bool val)
{
  this->snapper->setSnapToGrid(val);
}

void ArrangementDemoTabBase::setSnapToArrangement(bool val)
{
  this->snapper->setSnapToArrangement(val);
}

bool ArrangementDemoTabBase::isSnapToGridEnabled()
{
  return this->snapper->isSnapToGridEnabled();
}

bool ArrangementDemoTabBase::isSnapToArrangementEnabled()
{
  return this->snapper->isSnapToArrangementEnabled();
}

bool ArrangementDemoTabBase::isGridVisible()
{
  return this->gridGraphicsItem->isVisible();
}

CGAL::Qt::GraphicsViewCurveInputBase*
ArrangementDemoTabBase::getCurveInputCallback( ) const
{
  return this->curveInputCallback.get();
}

CGAL::Qt::Callback* ArrangementDemoTabBase::getDeleteCurveCallback( ) const
{
  return this->deleteCurveCallback.get();
}

CGAL::Qt::Callback* ArrangementDemoTabBase::getPointLocationCallback( ) const
{
  return this->pointLocationCallback.get();
}

VerticalRayShootCallbackBase*
ArrangementDemoTabBase::getVerticalRayShootCallback( ) const
{
  return this->verticalRayShootCallback.get();
}

CGAL::Qt::Callback* ArrangementDemoTabBase::getMergeEdgeCallback( ) const
{
  return this->mergeEdgeCallback.get();
}

SplitEdgeCallbackBase* ArrangementDemoTabBase::getSplitEdgeCallback( ) const
{
  return this->splitEdgeCallback.get();
}

EnvelopeCallbackBase* ArrangementDemoTabBase::getEnvelopeCallback( ) const
{
  return this->envelopeCallback.get();
}

FillFaceCallbackBase* ArrangementDemoTabBase::getFillFaceCallback( ) const
{
  return this->fillFaceCallback.get();
}

auto ArrangementDemoTabBase::getFillFaceColor() const -> QColor
{
  return this->fillFaceCallback->getColor();
}

void ArrangementDemoTabBase::setFillFaceColor(QColor color)
{
  this->fillFaceCallback->setColor(color);
}

CGAL::Qt::GraphicsViewNavigation*
ArrangementDemoTabBase::getGraphicsViewNavigation() const
{
  return this->navigation.get();
}

void ArrangementDemoTabBase::activateCurveInputCallback(CGAL::Qt::CurveType type)
{
  this->unhookCallbacks();

  this->curveInputCallback->setCurveType(type);
  this->getScene()->installEventFilter(this->curveInputCallback.get());
  this->activeCallback = this->curveInputCallback.get();
}

void ArrangementDemoTabBase::showLowerEnvelope(bool show)
{
  this->envelopeCallback->showLowerEnvelope(show);
  this->update();
}

void ArrangementDemoTabBase::showUpperEnvelope(bool show)
{
  this->envelopeCallback->showUpperEnvelope(show);
  this->update();
}

bool ArrangementDemoTabBase::isUpperEnvelopeShown()
{
  return this->envelopeCallback->isUpperEnvelopeShown();
}

bool ArrangementDemoTabBase::isLowerEnvelopeShown()
{
  return this->envelopeCallback->isLowerEnvelopeShown();
}

void ArrangementDemoTabBase::unhookCallbacks()
{
  if (this->activeCallback)
  {
    this->getScene()->removeEventFilter(this->activeCallback);

    activeCallback->reset();
    this->activeCallback = nullptr;
  }
}

void ArrangementDemoTabBase::unhookAndInstallEventFilter(
  CGAL::Qt::Callback* obj)
{
  this->unhookCallbacks();
  this->getScene()->installEventFilter(obj);
  this->activeCallback = obj;
}

void ArrangementDemoTabBase::activateDeleteCurveCallback()
{
  // TODO: Create different button for modes of delete
  if (
    this->activeCallback ==
    static_cast<CGAL::Qt::Callback*>(this->deleteCurveCallback.get()))
  {
    auto deleteMode = this->deleteCurveCallback->getDeleteMode();
    if (deleteMode == DeleteMode::DeleteOriginatingCuve)
      this->deleteCurveCallback->setDeleteMode(DeleteMode::DeleteEdge);
    else
      this->deleteCurveCallback->setDeleteMode(
        DeleteMode::DeleteOriginatingCuve);
  }
  else
  {
    this->unhookAndInstallEventFilter(this->deleteCurveCallback.get());
  }
}

void ArrangementDemoTabBase::activatePointLocationCallback()
{
  this->unhookAndInstallEventFilter(this->pointLocationCallback.get());
}

void ArrangementDemoTabBase::activateVerticalRayShootCallback(bool shootingUp)
{
  this->verticalRayShootCallback->setShootingUp(shootingUp);
  this->unhookAndInstallEventFilter(this->verticalRayShootCallback.get());
}

void ArrangementDemoTabBase::activateMergeEdgeCallback()
{
  this->unhookAndInstallEventFilter(this->mergeEdgeCallback.get());
}

void ArrangementDemoTabBase::activateSplitEdgeCallback()
{
  this->unhookAndInstallEventFilter(this->splitEdgeCallback.get());
}

void ArrangementDemoTabBase::activateFillFaceCallback()
{
  this->unhookAndInstallEventFilter(this->fillFaceCallback.get());
}

void ArrangementDemoTabBase::updatePreferences(const Preferences& pref)
{
  auto agi = this->getArrangementGraphicsItem();
  auto envelopeCallback = this->getEnvelopeCallback();
  auto verticalRayShootCallback = this->getVerticalRayShootCallback();
  auto splitEdgeCallback = this->getSplitEdgeCallback();
  auto gridGraphicsItem = this->getGridGraphicsItem();

  QPen edgesPen(QBrush(pref.edgeColor), pref.edgeWidth);
  edgesPen.setCosmetic(true);
  QPen verticesPen(QBrush(pref.vertexColor), pref.vertexRadius);
  verticesPen.setCosmetic(true);
  agi->setEdgesPen(edgesPen);
  agi->setVerticesPen(verticesPen);
  gridGraphicsItem->setAxesColor(pref.axesColor);
  gridGraphicsItem->setGridColor(pref.gridColor);
  envelopeCallback->setEnvelopeEdgeColor(pref.envelopeEdgeColor);
  envelopeCallback->setEnvelopeEdgeWidth(pref.envelopeEdgeWidth);
  envelopeCallback->setEnvelopeVertexColor(pref.envelopeVertexColor);
  envelopeCallback->setEnvelopeVertexRadius(pref.envelopeVertexRadius);
  verticalRayShootCallback->setEdgeColor(pref.verticalRayEdgeColor);
  verticalRayShootCallback->setEdgeWidth(pref.verticalRayEdgeWidth);
  splitEdgeCallback->setColor(pref.edgeColor);

  Q_EMIT modelChanged();
}

template <class Arr_>
ArrangementDemoTab<Arr_>::ArrangementDemoTab(
  QWidget* parent, std::unique_ptr<Arrangement> arrangement_) :
    Superclass(parent),
    arrangement(std::move(arrangement_))
{
  if (!this->arrangement) this->initArrangement();
  this->initComponents();
  this->setupCallbacks();
}

template <class Arr_>
demo_types::TraitsType ArrangementDemoTab<Arr_>::traitsType() const
{
  return demo_types::enumFromArrType<Arrangement>();
}

template <class Arr_>
ArrangementDemoTab<Arr_>::~ArrangementDemoTab()
{
}

template <class Arr_>
void ArrangementDemoTab<Arr_>::initArrangement()
{
  this->arrangement = std::make_unique<Arrangement>();
}

template <class Arr_>
void ArrangementDemoTab<Arr_>::initComponents()
{
  auto scene = this->getScene();

  this->gridGraphicsItem = new GridGraphicsItem();
  this->snapper = std::make_unique<PointSnapper<Arrangement>>(
    scene, this->gridGraphicsItem, this->arrangement.get());

  this->curveInputCallback =
    std::make_unique<ArrangementCurveInputCallback<Arrangement>>(
      this->arrangement.get(), this, scene);
  this->deleteCurveCallback =
    std::make_unique<DeleteCurveCallback<Arrangement>>(
      this->arrangement.get(), this);
  this->pointLocationCallback =
    std::make_unique<PointLocationCallback<Arrangement>>(
      this->arrangement.get(), this);
  this->verticalRayShootCallback =
    std::make_unique<VerticalRayShootCallback<Arrangement>>(
      this->arrangement.get(), this);
  this->mergeEdgeCallback = std::make_unique<MergeEdgeCallback<Arrangement>>(
    this->arrangement.get(), this);
  this->splitEdgeCallback = std::make_unique<SplitEdgeCallback<Arrangement>>(
    this->arrangement.get(), this);
  this->envelopeCallback = std::make_unique<EnvelopeCallback<Arrangement>>(
    this->arrangement.get(), this);
  this->fillFaceCallback = std::make_unique<FillFaceCallback<Arrangement>>(
    this->arrangement.get(), this);
  this->arrangementGraphicsItem =
    new CGAL::Qt::ArrangementGraphicsItem<Arrangement>(this->arrangement.get());

  this->curveInputCallback->setPointSnapper(snapper.get());
  this->splitEdgeCallback->setPointSnapper(snapper.get());

  scene->addItem(this->arrangementGraphicsItem);
  scene->addItem(this->gridGraphicsItem);

  this->arrangementGraphicsItem->setScene(scene);
  this->curveInputCallback->setScene(scene);
  this->deleteCurveCallback->setScene(scene);
  this->pointLocationCallback->setScene(scene);
  this->verticalRayShootCallback->setScene(scene);
  this->mergeEdgeCallback->setScene(scene);
  this->splitEdgeCallback->setScene(scene);
  this->envelopeCallback->setScene(scene);
  this->fillFaceCallback->setScene(scene);
}

template <class Arr_>
void ArrangementDemoTab<Arr_>::setupCallbacks()
{
  // set up callbacks
  QObject::connect(
    this->curveInputCallback.get(), SIGNAL(modelChanged()), this,
    SIGNAL(modelChanged()));
  QObject::connect(
    this->deleteCurveCallback.get(), SIGNAL(modelChanged()), this,
    SIGNAL(modelChanged()));
  QObject::connect(
    this->fillFaceCallback.get(), SIGNAL(modelChanged()), this,
    SIGNAL(modelChanged()));
  QObject::connect(
    this, SIGNAL(modelChanged()), this->arrangementGraphicsItem,
    SLOT(modelChanged()));
  QObject::connect(
    this, SIGNAL(modelChanged()), this->envelopeCallback.get(),
    SLOT(slotModelChanged()));
  QObject::connect(
    this->splitEdgeCallback.get(), SIGNAL(modelChanged()), this,
    SIGNAL(modelChanged()));
  QObject::connect(
    this->mergeEdgeCallback.get(), SIGNAL(modelChanged()), this,
    SIGNAL(modelChanged()));
  QObject::connect(
    this, SIGNAL(modelChanged()), this, SLOT(slotModelChanged()));
}

template <class Arr_>
CGAL::Object ArrangementDemoTab<Arr_>::getArrangement() const
{
  return CGAL::make_object(this->arrangement.get());
}

template <class Arr_>
void ArrangementDemoTab<Arr_>::slotModelChanged()
{
}

static CGAL::Bbox_2 reject_not_in_allowable_range(
  const CGAL::Bbox_2& box, const CGAL::Bbox_2& allowable_range)
{
  double xmin = std::numeric_limits<double>::infinity();
  double ymin = std::numeric_limits<double>::infinity();
  double xmax = -std::numeric_limits<double>::infinity();
  double ymax = -std::numeric_limits<double>::infinity();

  if (box.xmin() > allowable_range.xmin()) xmin = box.xmin();
  if (box.ymin() > allowable_range.ymin()) ymin = box.ymin();
  if (box.xmax() < allowable_range.xmax()) xmax = box.xmax();
  if (box.ymax() < allowable_range.ymax()) ymax = box.ymax();
  return {xmin, ymin, xmax, ymax};
}

static bool isFinite(const CGAL::Bbox_2& box)
{
  return !std::isinf(box.xmin()) && !std::isinf(box.xmax()) &&
         !std::isinf(box.ymin()) && !std::isinf(box.ymax());
}

static CGAL::Bbox_2 addMargins(const CGAL::Bbox_2& box)
{
  // add margin to bounding box
  double x_margin;
  double y_margin;
  if (box.xmin() == box.xmax() || box.ymin() == box.ymax())
  {
    static constexpr float const_margin = 50;
    x_margin = const_margin;
    y_margin = const_margin;
  }
  else
  {
    static constexpr double prop_margin = 0.10;
    x_margin = (box.xmax() - box.xmin()) * prop_margin;
    y_margin = (box.ymax() - box.ymin()) * prop_margin;
  }
  return {
    box.xmin() - x_margin, box.ymin() - y_margin,
    box.xmax() + x_margin, box.ymax() + y_margin};
}

// TODO: clean this up, it's ugly!
template <class Arr_>
CGAL::Bbox_2
findOtherInterestingPoints(const std::unique_ptr<Arr_>&, const CGAL::Bbox_2&)
{
  return {};
}

#ifdef CGAL_USE_CORE
template <typename Coefficient_>
static const auto&
getXyCurves(const CGAL::Arr_algebraic_segment_traits_2<Coefficient_>* traits)
{
  // the traits object is only needed the first time
  // this assumes that X_monotone_curves created from the first traits object
  // will work with arrangements with a different object
  using Traits = CGAL::Arr_algebraic_segment_traits_2<Coefficient_>;
  static std::vector<typename Traits::X_monotone_curve_2> xy_curves;
  if (xy_curves.empty())
  {
    typedef typename Traits::Polynomial_2 Polynomial_2;
    auto construct_curve = traits->construct_curve_2_object();
    auto make_x_monotone = traits->make_x_monotone_2_object();

    Polynomial_2 x = CGAL::shift(Polynomial_2(1), 1, 0);
    Polynomial_2 y = CGAL::shift(Polynomial_2(1), 1, 1);
    auto x_cv = construct_curve(x);
    auto y_cv = construct_curve(y);

    std::vector<CGAL::Object> arcs;
    make_x_monotone(x_cv, std::back_inserter(arcs));
    make_x_monotone(y_cv, std::back_inserter(arcs));
    for (auto& arc_obj : arcs)
    {
      typename Traits::X_monotone_curve_2 arc;
      CGAL::assign(arc, arc_obj);
      xy_curves.push_back(arc);
    }
  }
  return xy_curves;
}

template <>
CGAL::Bbox_2 findOtherInterestingPoints<demo_types::Alg_seg_arr>(
  const std::unique_ptr<demo_types::Alg_seg_arr>& arr,
  const CGAL::Bbox_2& allowable_range)
{
  using Traits = demo_types::Alg_seg_traits;
  CGAL::Bbox_2 bb = {};
  std::vector<CGAL::Object> intersections;
  for (auto it = arr->edges_begin(); it != arr->edges_end(); ++it)
  {
    for (auto& arc : getXyCurves(arr->traits()))
      if (arc.is_vertical() != it->curve().is_vertical())
        it->curve().intersections(arc, std::back_inserter(intersections));
  }
  for (auto it = intersections.begin(); it != intersections.end(); it++)
  {
    std::pair<typename Traits::Point_2, unsigned int> point_multiplicity;
    CGAL::assign(point_multiplicity, *it);
    auto& point = point_multiplicity.first;
    if (point.location() == CGAL::ARR_INTERIOR)
    {
      auto xy = point.to_double();
      bb += reject_not_in_allowable_range(
        {xy.first, xy.second, xy.first, xy.second}, allowable_range);
    }
  }
  return bb;
}
#endif // CGAL_USE_CORE

template <class Arr_>
void ArrangementDemoTab<Arr_>::adjustViewport()
{
  QRectF scene_rect = this->getScene()->sceneRect();
  CGAL::Bbox_2 scene_bbox = {
    scene_rect.left(), scene_rect.top(), scene_rect.right(),
    scene_rect.bottom()};

  ConstructBoundingBox<Traits> construct_bounding_box;
  CGAL::Bbox_2 bb = {};
  for (auto it = this->arrangement->edges_begin();
       it != this->arrangement->edges_end(); ++it)
  {
    bb += reject_not_in_allowable_range(
      construct_bounding_box(it->curve()), scene_bbox);
  }
  for (auto it = this->arrangement->vertices_begin();
       it != this->arrangement->vertices_end(); ++it)
  {
    bb += reject_not_in_allowable_range(
      construct_bounding_box(it->point()), scene_bbox);
  }

  if (!isFinite(bb))
    bb += findOtherInterestingPoints(this->arrangement, scene_bbox);

  // ideally this should happen only if the arrangement is empty
  // (if findOtherInterestingPoints worked for all arrangement types)
  if (!isFinite(bb))
    bb += {0, 0, 0, 0};

  bb = addMargins(bb);

  auto viewportRect =
    QRectF(bb.xmin(), bb.ymin(), bb.xmax() - bb.xmin(), bb.ymax() - bb.ymin());

  this->graphicsView->resetTransform();
  this->graphicsView->fitInView(viewportRect, Qt::KeepAspectRatio);

  Q_EMIT modelChanged();
}

ARRANGEMENT_DEMO_SPECIALIZE_ARR(ArrangementDemoTab)

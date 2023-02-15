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

#include "ArrangementDemoTab.h"
#include "ArrangementDemoGraphicsView.h"
#include "ArrangementGraphicsItem.h"
#include "DeleteCurveCallback.h"
#include "EnvelopeCallback.h"
#include "FillFaceCallback.h"
#include "GraphicsViewCurveInput.h"
#include "GridGraphicsItem.h"
#include "MergeEdgeCallback.h"
#include "PointLocationCallback.h"
#include "PointSnapper.h"
#include "SplitEdgeCallback.h"
#include "Utils/Utils.h"
#include "VerticalRayShootCallback.h"

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include <QGridLayout>
#include <limits>


ArrangementDemoTab::ArrangementDemoTab(
  QWidget* parent, demo_types::TraitsType tt, CGAL::Object arrangement_) :
    QWidget(parent),
    GraphicsSceneMixin(new QGraphicsScene()),
    graphicsView(new ArrangementDemoGraphicsView(this)),
    layout(new QGridLayout(this)),
    navigation(std::make_unique<CGAL::Qt::GraphicsViewNavigation>()),
    activeCallback(nullptr), arrangementGraphicsItem(nullptr),
    gridGraphicsItem(nullptr), ttype(tt), arrangement(arrangement_)
{
  this->setupUi();
  if (!this->arrangement) this->initArrangement();
  this->initComponents();
  this->setupCallbacks();
}

ArrangementDemoTab::~ArrangementDemoTab()
{
  this->unhookCallbacks();
  deleteArrangement(this->traitsType(), this->getArrangement());
}

void ArrangementDemoTab::setupUi()
{
  auto scene = this->getScene();

  this->layout->addWidget(this->graphicsView, 0, 0);
  this->graphicsView->setScene(scene);

  double MAX_WIDTH = (std::numeric_limits<double>::max)() / 1048576;

  double xymin = -MAX_WIDTH / 2;
  double wh = MAX_WIDTH;
  scene->setSceneRect(xymin, xymin, wh, wh);

  this->getView()->installEventFilter(this->navigation.get());
  this->getView()->viewport()->installEventFilter(this->navigation.get());
}

QGraphicsView* ArrangementDemoTab::getView() const
{
  return this->graphicsView;
}

CGAL::Qt::ArrangementGraphicsItemBase*
ArrangementDemoTab::getArrangementGraphicsItem() const
{
  return this->arrangementGraphicsItem;
}

GridGraphicsItem* ArrangementDemoTab::getGridGraphicsItem() const
{
  return this->gridGraphicsItem;
}

void ArrangementDemoTab::showGrid(bool val)
{
  this->gridGraphicsItem->setVisible(val);
}

void ArrangementDemoTab::setSnapToGrid(bool val)
{
  this->snapper->setSnapToGrid(val);
}

void ArrangementDemoTab::setSnapToArrangement(bool val)
{
  this->snapper->setSnapToArrangement(val);
}

bool ArrangementDemoTab::isSnapToGridEnabled()
{
  return this->snapper->isSnapToGridEnabled();
}

bool ArrangementDemoTab::isSnapToArrangementEnabled()
{
  return this->snapper->isSnapToArrangementEnabled();
}

bool ArrangementDemoTab::isGridVisible()
{
  return this->gridGraphicsItem->isVisible();
}

CGAL::Qt::GraphicsViewCurveInputBase*
ArrangementDemoTab::getCurveInputCallback() const
{
  return this->curveInputCallback.get();
}

CGAL::Qt::Callback* ArrangementDemoTab::getDeleteCurveCallback() const
{
  return this->deleteCurveCallback.get();
}

PointLocationCallbackBase*
ArrangementDemoTab::getPointLocationCallback() const
{
  return this->pointLocationCallback.get();
}

VerticalRayShootCallbackBase*
ArrangementDemoTab::getVerticalRayShootCallback() const
{
  return this->verticalRayShootCallback.get();
}

MergeEdgeCallbackBase* ArrangementDemoTab::getMergeEdgeCallback() const
{
  return this->mergeEdgeCallback.get();
}

SplitEdgeCallbackBase* ArrangementDemoTab::getSplitEdgeCallback() const
{
  return this->splitEdgeCallback.get();
}

EnvelopeCallbackBase* ArrangementDemoTab::getEnvelopeCallback() const
{
  return this->envelopeCallback.get();
}

FillFaceCallbackBase* ArrangementDemoTab::getFillFaceCallback() const
{
  return this->fillFaceCallback.get();
}

auto ArrangementDemoTab::getFillFaceColor() const -> QColor
{
  return this->fillFaceCallback->getColor();
}

void ArrangementDemoTab::setFillFaceColor(QColor color)
{
  this->fillFaceCallback->setColor(color);
}

CGAL::Qt::GraphicsViewNavigation*
ArrangementDemoTab::getGraphicsViewNavigation() const
{
  return this->navigation.get();
}

void ArrangementDemoTab::activateCurveInputCallback(
  CGAL::Qt::CurveType type)
{
  this->unhookCallbacks();

  this->curveInputCallback->setCurveType(type);
  this->getScene()->installEventFilter(this->curveInputCallback.get());
  this->activeCallback = this->curveInputCallback.get();
}

void ArrangementDemoTab::showLowerEnvelope(bool show)
{
  this->envelopeCallback->showLowerEnvelope(show);
  this->update();
}

void ArrangementDemoTab::showUpperEnvelope(bool show)
{
  this->envelopeCallback->showUpperEnvelope(show);
  this->update();
}

bool ArrangementDemoTab::isUpperEnvelopeShown()
{
  return this->envelopeCallback->isUpperEnvelopeShown();
}

bool ArrangementDemoTab::isLowerEnvelopeShown()
{
  return this->envelopeCallback->isLowerEnvelopeShown();
}

void ArrangementDemoTab::unhookCallbacks()
{
  if (this->activeCallback)
  {
    this->getScene()->removeEventFilter(this->activeCallback);

    activeCallback->reset();
    this->activeCallback = nullptr;
  }
}

void ArrangementDemoTab::unhookAndInstallEventFilter(
  CGAL::Qt::Callback* obj)
{
  this->unhookCallbacks();
  this->getScene()->installEventFilter(obj);
  this->activeCallback = obj;
}

void ArrangementDemoTab::activateDeleteCurveCallback()
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

void ArrangementDemoTab::activatePointLocationCallback()
{
  this->unhookAndInstallEventFilter(this->pointLocationCallback.get());
}

void ArrangementDemoTab::activateVerticalRayShootCallback(bool shootingUp)
{
  this->verticalRayShootCallback->setShootingUp(shootingUp);
  this->unhookAndInstallEventFilter(this->verticalRayShootCallback.get());
}

void ArrangementDemoTab::activateMergeEdgeCallback()
{
  this->unhookAndInstallEventFilter(this->mergeEdgeCallback.get());
}

void ArrangementDemoTab::activateSplitEdgeCallback()
{
  this->unhookAndInstallEventFilter(this->splitEdgeCallback.get());
}

void ArrangementDemoTab::activateFillFaceCallback()
{
  this->unhookAndInstallEventFilter(this->fillFaceCallback.get());
}

void ArrangementDemoTab::updatePreferences(const Preferences& pref)
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

demo_types::TraitsType ArrangementDemoTab::traitsType() const
{
  return ttype;
}

void ArrangementDemoTab::initArrangement()
{
  this->arrangement = createArrangement(this->traitsType());
}

void ArrangementDemoTab::initComponents()
{
  auto scene = this->getScene();

  this->gridGraphicsItem = new GridGraphicsItem();

  this->snapper = std::unique_ptr<PointSnapperBase>(PointSnapperBase::create(
    this->traitsType(), scene, this->gridGraphicsItem, this->getArrangement()));

  this->curveInputCallback =
    std::unique_ptr<CGAL::Qt::GraphicsViewCurveInputBase>(
      CGAL::Qt::GraphicsViewCurveInputBase::create(
        this->traitsType(), this->getArrangement(), this, scene));
  this->deleteCurveCallback =
    std::unique_ptr<DeleteCurveCallbackBase>(DeleteCurveCallbackBase::create(
      this->traitsType(), this->getArrangement(), this));
  this->pointLocationCallback = std::unique_ptr<PointLocationCallbackBase>(
    PointLocationCallbackBase::create(
      this->traitsType(), this->getArrangement(), this));
  this->verticalRayShootCallback =
    std::unique_ptr<VerticalRayShootCallbackBase>(
      VerticalRayShootCallbackBase::create(
        this->traitsType(), this->getArrangement(), this));
  this->mergeEdgeCallback =
    std::unique_ptr<MergeEdgeCallbackBase>(MergeEdgeCallbackBase::create(
      this->traitsType(), this->getArrangement(), this));
  this->splitEdgeCallback =
    std::unique_ptr<SplitEdgeCallbackBase>(SplitEdgeCallbackBase::create(
      this->traitsType(), this->getArrangement(), this));
  this->envelopeCallback =
    std::unique_ptr<EnvelopeCallbackBase>(EnvelopeCallbackBase::create(
      this->traitsType(), this->getArrangement(), this));
  this->fillFaceCallback =
    std::unique_ptr<FillFaceCallbackBase>(FillFaceCallbackBase::create(
      this->traitsType(), this->getArrangement(), this));
  this->arrangementGraphicsItem = CGAL::Qt::ArrangementGraphicsItemBase::create(
    this->traitsType(), this->getArrangement());

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

void ArrangementDemoTab::setupCallbacks()
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
}

CGAL::Object ArrangementDemoTab::getArrangement() const
{
  return this->arrangement;
}

void ArrangementDemoTab::adjustViewport()
{
  this->graphicsView->resetTransform();
  this->graphicsView->fitInView(
    this->arrangementGraphicsItem->getInterestingViewport(),
    Qt::KeepAspectRatio);

  Q_EMIT modelChanged();
}

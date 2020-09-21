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

#ifndef ARRANGEMENT_DEMO_TAB_H
#define ARRANGEMENT_DEMO_TAB_H

#include "GraphicsSceneMixin.h"
#include <CGAL/Object.h>
#include <QWidget>
#include <QColor>
#include <memory>

class QGraphicsScene;
class QGridLayout;
class ArrangementDemoGraphicsView;
class VerticalRayShootCallbackBase;
class EnvelopeCallbackBase;
class FillFaceCallbackBase;
class SplitEdgeCallbackBase;
class DeleteCurveCallbackBase;
class GridGraphicsItem;
class PointSnapperBase;

namespace CGAL
{
namespace Qt
{
class ArrangementGraphicsItemBase;
class Callback;
class ArrangementGraphicsItemBase;
class GraphicsViewCurveInputBase;
class GraphicsViewNavigation;
enum class CurveType;
} // namespace Qt
} // namespace CGAL

namespace demo_types
{
enum class TraitsType;
}

class ArrangementDemoTabBase : public QWidget, public GraphicsSceneMixin
{
  Q_OBJECT

Q_SIGNALS:
  void modelChanged( );

public:
  ArrangementDemoTabBase( QWidget* parent );
  virtual ~ArrangementDemoTabBase( );

  virtual CGAL::Object getArrangement() const = 0;
  virtual void adjustViewport() = 0;
  virtual demo_types::TraitsType traitsType() const = 0;

  QGraphicsView* getView() const;
  void showGrid(bool);
  bool isGridVisible();
  void setSnapToGrid(bool);
  void setSnapToArrangement(bool);
  bool isSnapToGridEnabled();
  bool isSnapToArrangementEnabled();

  auto getArrangementGraphicsItem() const
    -> CGAL::Qt::ArrangementGraphicsItemBase*;
  auto getGridGraphicsItem() const -> GridGraphicsItem*;
  auto getCurveInputCallback() const -> CGAL::Qt::GraphicsViewCurveInputBase*;
  auto getDeleteCurveCallback() const -> CGAL::Qt::Callback*;
  auto getPointLocationCallback() const -> CGAL::Qt::Callback*;
  auto getVerticalRayShootCallback() const -> VerticalRayShootCallbackBase*;
  auto getMergeEdgeCallback() const -> CGAL::Qt::Callback*;
  auto getSplitEdgeCallback() const -> SplitEdgeCallbackBase*;
  auto getEnvelopeCallback() const -> EnvelopeCallbackBase*;
  auto getFillFaceCallback() const -> FillFaceCallbackBase*;
  auto getGraphicsViewNavigation() const -> CGAL::Qt::GraphicsViewNavigation*;
  auto getFillFaceColor() const -> QColor;
  void setFillFaceColor(QColor);
  void activateCurveInputCallback(CGAL::Qt::CurveType);
  void activateDeleteCurveCallback();
  void activatePointLocationCallback();
  void activateVerticalRayShootCallback(bool);
  void activateMergeEdgeCallback();
  void activateSplitEdgeCallback();
  void activateEnvelopeCallback();
  void activateFillFaceCallback();
  void activatePanCallback();
  void showLowerEnvelope(bool);
  void showUpperEnvelope(bool);
  bool isUpperEnvelopeShown();
  bool isLowerEnvelopeShown();
  void unhookCallbacks();

  struct Preferences
  {
    QColor edgeColor = {};
    QColor vertexColor = {};
    QColor envelopeEdgeColor = {};
    QColor envelopeVertexColor = {};
    QColor verticalRayEdgeColor = {};
    QColor axesColor = {};
    QColor gridColor = {};
    uint32_t edgeWidth = 0;
    uint32_t vertexRadius = 0;
    uint32_t envelopeEdgeWidth = 0;
    uint32_t envelopeVertexRadius = 0;
    uint32_t verticalRayEdgeWidth = 0;
  };
  void updatePreferences(const Preferences&);

protected Q_SLOTS:
  virtual void slotModelChanged() = 0;

protected:
  virtual void setupUi( );
  void unhookAndInstallEventFilter(CGAL::Qt::Callback*);

  ArrangementDemoGraphicsView* graphicsView;
  QGridLayout* layout;

  std::unique_ptr<CGAL::Qt::GraphicsViewCurveInputBase> curveInputCallback;
  std::unique_ptr<DeleteCurveCallbackBase> deleteCurveCallback;
  std::unique_ptr<CGAL::Qt::Callback> pointLocationCallback;
  std::unique_ptr<VerticalRayShootCallbackBase> verticalRayShootCallback;
  std::unique_ptr<CGAL::Qt::Callback> mergeEdgeCallback;
  std::unique_ptr<SplitEdgeCallbackBase> splitEdgeCallback;
  std::unique_ptr<EnvelopeCallbackBase> envelopeCallback;
  std::unique_ptr<FillFaceCallbackBase> fillFaceCallback;
  std::unique_ptr<PointSnapperBase> snapper;
  std::unique_ptr<CGAL::Qt::GraphicsViewNavigation> navigation;

  CGAL::Qt::Callback* activeCallback;
  CGAL::Qt::ArrangementGraphicsItemBase* arrangementGraphicsItem;
  GridGraphicsItem* gridGraphicsItem;
}; // class ArrangementDemoTabBase

template < class Arr_ >
class ArrangementDemoTab : public ArrangementDemoTabBase
{
public:
  typedef ArrangementDemoTabBase Superclass;
  typedef Arr_ Arrangement;
  typedef typename Arrangement::Geometry_traits_2 Traits;

  ArrangementDemoTab(
    QWidget* parent, std::unique_ptr<Arrangement> arrangement_ = nullptr);
  ~ArrangementDemoTab();
  void adjustViewport() override;
  CGAL::Object getArrangement() const override;
  demo_types::TraitsType traitsType() const override;

private:
  void initArrangement();
  void initComponents();
  void setupCallbacks();

protected:
  void slotModelChanged() override;

protected:
  std::unique_ptr<Arrangement> arrangement;

}; // class ArrangementDemoTab

#endif // ARRANGEMENT_DEMO_TAB_H

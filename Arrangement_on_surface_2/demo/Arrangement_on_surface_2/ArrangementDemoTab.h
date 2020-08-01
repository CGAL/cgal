// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef ARRANGEMENT_DEMO_TAB_H
#define ARRANGEMENT_DEMO_TAB_H

#include "GraphicsSceneMixin.h"
#include <CGAL/Object.h>
#include <QWidget>
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
enum class CurveType;
} // namespace Qt
} // namespace CGAL

class ArrangementDemoTabBase : public QWidget, public QGraphicsSceneMixin
{
  Q_OBJECT

Q_SIGNALS:
  void modelChanged( );

public:
  ArrangementDemoTabBase( QWidget* parent );
  virtual ~ArrangementDemoTabBase( );

  QGraphicsView* getView() const;
  virtual CGAL::Object getArrangement() const = 0;
  virtual void adjustViewport() = 0;
  void showGrid(bool);
  bool isGridVisible();
  void setSnapToGrid(bool);
  void setSnapToArrangement(bool);
  bool isSnapToGridEnabled();
  bool isSnapToArrangementEnabled();

  CGAL::Qt::ArrangementGraphicsItemBase* getArrangementGraphicsItem() const;
  GridGraphicsItem* getGridGraphicsItem() const;
  CGAL::Qt::GraphicsViewCurveInputBase* getCurveInputCallback() const;
  CGAL::Qt::Callback* getDeleteCurveCallback() const;
  CGAL::Qt::Callback* getPointLocationCallback() const;
  VerticalRayShootCallbackBase* getVerticalRayShootCallback() const;
  CGAL::Qt::Callback* getMergeEdgeCallback() const;
  SplitEdgeCallbackBase* getSplitEdgeCallback() const;
  EnvelopeCallbackBase* getEnvelopeCallback() const;
  FillFaceCallbackBase* getFillFaceCallback() const;

  void activateCurveInputCallback(CGAL::Qt::CurveType);
  void activateDeleteCurveCallback();
  void activatePointLocationCallback();
  void activateVerticalRayShootCallback(bool);
  void activateMergeEdgeCallback();
  void activateSplitEdgeCallback();
  void activateEnvelopeCallback();
  void activateFillFaceCallback();
  void activatePanCallback();
  void unhookCallbacks();

protected Q_SLOTS:
  virtual void slotModelChanged() = 0;

protected:
  virtual void setupUi( );
  void unhookAndInstallEventFilter(QObject*);

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

  QObject* activeCallback;

  CGAL::Qt::ArrangementGraphicsItemBase* arrangementGraphicsItem;
  GridGraphicsItem* gridGraphicsItem;
}; // class ArrangementDemoTabBase

template < class Arr_ >
class ArrangementDemoTab : public ArrangementDemoTabBase
{
public:
  typedef ArrangementDemoTabBase Superclass;
  typedef Arr_ Arrangement;

  ArrangementDemoTab(
    QWidget* parent, std::unique_ptr<Arrangement> arrangement_ = nullptr);
  ~ArrangementDemoTab();
  void adjustViewport();

private:
  void initArrangement();
  void initComponents();
  void setupCallbacks();
  CGAL::Object getArrangement() const override;

protected:
  void slotModelChanged() override;

protected:
  std::unique_ptr<Arrangement> arrangement;

}; // class ArrangementDemoTab

#endif // ARRANGEMENT_DEMO_TAB_H

// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef ARRANGEMENT_DEMO_WINDOW_H
#define ARRANGEMENT_DEMO_WINDOW_H

#include <CGAL/Qt/DemosMainWindow.h>

#include <Qt>
#include <utility>
#include <vector>

namespace Ui
{
class ArrangementDemoWindow;
}

namespace CGAL
{
class Object;
}

class ArrangementDemoTabBase;
class QActionGroup;

class ArrangementDemoWindow : public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT

public:
  typedef enum TraitsType
  {
    SEGMENT_TRAITS,
    POLYLINE_TRAITS,
    CONIC_TRAITS,
    LINEAR_TRAITS,
    ALGEBRAIC_TRAITS,
    BEZIER_TRAITS,
    NONE,
  } TraitsType;

  ArrangementDemoWindow(QWidget* parent = nullptr);
  ~ArrangementDemoWindow();

  std::vector<CGAL::Object> getArrangements() const;
  std::vector<QString> getTabLabels() const;
  std::pair<ArrangementDemoTabBase*, TraitsType> getCurrentTab();

public Q_SLOTS:
  void updateEnvelope(QAction*);
  void updateInputType(QAction*);
  void on_actionNewTab_triggered();
  void on_actionSaveAs_triggered();
  void on_actionOpen_triggered();
  void on_actionQuit_triggered();
  void on_tabWidget_currentChanged(int);
  void on_actionOverlay_triggered();
  void on_actionCloseTab_triggered();
  void on_actionZoomIn_triggered();
  void on_actionZoomOut_triggered();
  void on_actionZoomReset_triggered();
  void on_actionPreferences_triggered();
  void on_actionFillColor_triggered();
  void on_actionInsert_toggled(bool);
  void on_actionDrag_toggled(bool);
  void on_actionDelete_toggled(bool);
  void on_actionDelete_triggered();
  void on_actionPointLocation_toggled(bool);
  void on_actionRayShootingUp_toggled(bool);
  void on_actionRayShootingDown_toggled(bool);
  void on_actionShowGrid_toggled(bool);
  void on_actionGridSnapMode_toggled(bool);
  void on_actionArrangementSnapMode_toggled(bool);
  void on_actionMerge_toggled(bool);
  void on_actionSplit_toggled(bool);
  void on_actionFill_toggled(bool);
  void on_actionAddAlgebraicCurve_triggered();

Q_SIGNALS:
  void modelChanged();

protected:
  void setupUi();
  ArrangementDemoTabBase* makeTab(TraitsType);
  void addTab(ArrangementDemoTabBase*, QString, TraitsType);
  void resetCallbackState(ArrangementDemoTabBase*);
  void resetActionGroups(ArrangementDemoTabBase*, TraitsType);
  void hideInsertMethods();
  void showInsertMethods();
  void updateFillColorSwatch(ArrangementDemoTabBase*);
  QString makeTabLabel(TraitsType);
  ArrangementDemoTabBase* openArrFile(QString filename);

  template <class ArrType>
  ArrangementDemoTabBase*
  makeTab(std::unique_ptr<ArrType> arr, QString, TraitsType);

  template <class ArrType>
  void makeOverlayTab(ArrType* arr1, ArrType* arr2, TraitsType);

private:
  Ui::ArrangementDemoWindow* ui;
  std::vector<std::pair<ArrangementDemoTabBase*, TraitsType>> tabs;
  QActionGroup* modeGroup;
  QActionGroup* envelopeGroup;
  QActionGroup* inputTypeGroup;
  int tabLabelCounter;
};

#endif // ARRANGEMENT_DEMO_WINDOW_H

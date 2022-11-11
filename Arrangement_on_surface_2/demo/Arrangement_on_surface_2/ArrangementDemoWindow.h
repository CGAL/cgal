// Copyright (c) 2012, 2020 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#ifndef ARRANGEMENT_DEMO_WINDOW_H
#define ARRANGEMENT_DEMO_WINDOW_H

#include <CGAL/Object.h>
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
namespace Qt
{
class GraphicsViewNavigation;
}
} // namespace CGAL

namespace demo_types
{
enum class TraitsType;
}

class ArrangementDemoTab;
class QActionGroup;

class ArrangementDemoWindow : public CGAL::Qt::DemosMainWindow
{
  Q_OBJECT

public:
  ArrangementDemoWindow(QWidget* parent = nullptr);
  ~ArrangementDemoWindow();

  ArrangementDemoTab* getCurrentTab();

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
  void on_actionAddRationalCurve_triggered();

Q_SIGNALS:
  void modelChanged();

protected:
  void setupUi();
  ArrangementDemoTab* makeTab(
    demo_types::TraitsType, QString label = {}, CGAL::Object arr_obj = {});
  void addTab(ArrangementDemoTab*, QString);
  void resetCallbackState(ArrangementDemoTab*);
  void resetActionGroups(ArrangementDemoTab*, demo_types::TraitsType);
  void hideInsertMethods();
  void showInsertMethods(demo_types::TraitsType);
  void updateFillColorSwatch(ArrangementDemoTab*);
  QString makeTabLabel(demo_types::TraitsType);
  ArrangementDemoTab* openArrFile(QString filename);

private:
  Ui::ArrangementDemoWindow* ui;
  std::vector<ArrangementDemoTab*> tabs;
  QActionGroup* modeGroup;
  QActionGroup* envelopeGroup;
  QActionGroup* inputTypeGroup;
  int tabLabelCounter;
};

#endif // ARRANGEMENT_DEMO_WINDOW_H

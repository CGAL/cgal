// Copyright (c) 2012, 2020  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>
//            Ahmed Essam <theartful.ae@gmail.com>

#include "ArrangementDemoWindow.h"
#include "AlgebraicCurveInputDialog.h"
#include "RationalCurveInputDialog.h"
#include "ArrangementDemoPropertiesDialog.h"
#include "ArrangementDemoTab.h"
#include "NewTabDialog.h"
#include "OverlayDialog.h"
#include "GraphicsViewCurveInput.h"
#include "CurveInputMethods.h"
#include "ArrangementTypesUtils.h"
#include "ArrangementIO.h"
#include "Utils/Utils.h"

#include <QActionGroup>
#include <QColorDialog>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QString>
#include <QGraphicsView>
#include <QLabel>

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "ui_ArrangementDemoWindow.h"

using TraitsType = demo_types::TraitsType;

ArrangementDemoWindow::ArrangementDemoWindow(QWidget* parent) :
    CGAL::Qt::DemosMainWindow(parent), ui(new Ui::ArrangementDemoWindow),
    tabLabelCounter(1)
{
  this->setupUi();

  QObject::connect(
    this->inputTypeGroup, SIGNAL(triggered(QAction*)), this,
    SLOT(updateInputType(QAction*)));

  QObject::connect(
    this->envelopeGroup, SIGNAL(triggered(QAction*)), this,
    SLOT(updateEnvelope(QAction*)));

  this->makeTab(TraitsType::SEGMENT_TRAITS);

  // Call inherited functions
  this->setupStatusBar();
  // this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/Arrangement_on_surface_2/about.html");
  this->addAboutCGAL();
}

ArrangementDemoWindow::~ArrangementDemoWindow() { }

void ArrangementDemoWindow::setupUi()
{
  this->ui->setupUi(this);

  this->modeGroup = new QActionGroup(this);
  this->modeGroup->addAction(this->ui->actionDrag);
  this->modeGroup->addAction(this->ui->actionInsert);
  this->modeGroup->addAction(this->ui->actionDelete);
  this->modeGroup->addAction(this->ui->actionPointLocation);
  this->modeGroup->addAction(this->ui->actionRayShootingUp);
  this->modeGroup->addAction(this->ui->actionRayShootingDown);
  this->modeGroup->addAction(this->ui->actionMerge);
  this->modeGroup->addAction(this->ui->actionSplit);
  this->modeGroup->addAction(this->ui->actionFill);

  this->envelopeGroup = new QActionGroup(this);
  this->envelopeGroup->addAction(this->ui->actionLowerEnvelope);
  this->envelopeGroup->addAction(this->ui->actionUpperEnvelope);
  this->envelopeGroup->setExclusive(false);

  this->inputTypeGroup = new QActionGroup(this);
  this->inputTypeGroup->addAction(this->ui->actionSegment);
  this->inputTypeGroup->addAction(this->ui->actionPolyline);
  this->inputTypeGroup->addAction(this->ui->actionCircle);
  this->inputTypeGroup->addAction(this->ui->actionEllipse);
  this->inputTypeGroup->addAction(this->ui->actionConicThreePoint);
  this->inputTypeGroup->addAction(this->ui->actionConicFivePoint);
  this->inputTypeGroup->addAction(this->ui->actionRay);
  this->inputTypeGroup->addAction(this->ui->actionLine);
  this->inputTypeGroup->addAction(this->ui->actionBezier);
}

QString ArrangementDemoWindow::makeTabLabel(TraitsType tt)
{
  static const char* typeNames[] = {
    "Segment", "Polyline", "Linear", "Conic", "Algebraic",
    "Bezier",  "Rational Functions"
  };
  return QString("%1 - %2")
    .arg(this->tabLabelCounter++)
    .arg(typeNames[static_cast<int>(tt)]);
}

ArrangementDemoTab* ArrangementDemoWindow::makeTab(
  TraitsType tt, QString tabLabel, CGAL::Object arr_obj)
{
  ArrangementDemoTab* demoTab = new ArrangementDemoTab(this, tt, arr_obj);
  if (tabLabel.isNull())
    tabLabel = makeTabLabel(tt);
  this->addTab(demoTab, tabLabel);
  return demoTab;
}

void ArrangementDemoWindow::addTab(
  ArrangementDemoTab* demoTab, QString tabLabel)
{
  this->tabs.push_back(demoTab);

  this->ui->tabWidget->addTab(demoTab, tabLabel);
  this->ui->tabWidget->setCurrentWidget(demoTab);

  // avoid using addNavigation since it will cause memory leaks when there are
  // multiple tabs (no way to delete the navigation item!)
  QObject::connect(
    demoTab->getGraphicsViewNavigation(), SIGNAL(mouseCoordinates(QString)),
    xycoord, SLOT(setText(QString)));
}

void ArrangementDemoWindow::resetActionGroups(
  ArrangementDemoTab* tab, TraitsType tt)
{
  this->hideInsertMethods();
  this->ui->actionInsert->setChecked(false);
  this->ui->actionDrag->setChecked(false);
  this->ui->actionLowerEnvelope->setChecked(false);
  this->ui->actionUpperEnvelope->setChecked(false);
  this->ui->actionShowGrid->setChecked(tab->isGridVisible());
  this->ui->actionGridSnapMode->setChecked(tab->isSnapToGridEnabled());
  this->ui->actionArrangementSnapMode->setChecked(
    tab->isSnapToArrangementEnabled());
  this->ui->actionLowerEnvelope->setChecked(tab->isLowerEnvelopeShown());
  this->ui->actionUpperEnvelope->setChecked(tab->isUpperEnvelopeShown());

  if (
    tt != TraitsType::SEGMENT_TRAITS && tt != TraitsType::POLYLINE_TRAITS &&
    tt != TraitsType::LINEAR_TRAITS)
    this->ui->actionArrangementSnapMode->setVisible(false);
  else
    this->ui->actionArrangementSnapMode->setVisible(true);

  // default action group is scrolling
  this->ui->actionDrag->activate(QAction::Trigger);
}

void ArrangementDemoWindow::resetCallbackState(ArrangementDemoTab* tab)
{
  if (tab)
  {
    tab->getView()->setDragMode(QGraphicsView::NoDrag);
    tab->unhookCallbacks();
  }
}

void ArrangementDemoWindow::updateEnvelope(QAction* newMode)
{
  auto currentTab = this->getCurrentTab();
  if (!currentTab) return;

  bool show = newMode->isChecked();
  if (newMode == this->ui->actionLowerEnvelope)
    currentTab->showLowerEnvelope(show);
  else if (newMode == this->ui->actionUpperEnvelope)
    currentTab->showUpperEnvelope(show);
}

void ArrangementDemoWindow::hideInsertMethods()
{
  this->ui->actionAddAlgebraicCurve->setVisible(false);
  this->ui->actionAddAlgebraicCurve->setChecked(false);
  this->ui->actionAddRationalCurve->setVisible(false);
  this->ui->actionAddRationalCurve->setChecked(false);
  this->ui->actionSegment->setVisible(false);
  this->ui->actionSegment->setChecked(false);
  this->ui->actionCircle->setVisible(false);
  this->ui->actionCircle->setChecked(false);
  this->ui->actionEllipse->setVisible(false);
  this->ui->actionEllipse->setChecked(false);
  this->ui->actionConicThreePoint->setVisible(false);
  this->ui->actionConicThreePoint->setChecked(false);
  this->ui->actionConicFivePoint->setVisible(false);
  this->ui->actionConicFivePoint->setChecked(false);
  this->ui->actionRay->setVisible(false);
  this->ui->actionRay->setChecked(false);
  this->ui->actionLine->setVisible(false);
  this->ui->actionLine->setChecked(false);
  this->ui->actionPolyline->setVisible(false);
  this->ui->actionPolyline->setChecked(false);
  this->ui->actionBezier->setVisible(false);
  this->ui->actionBezier->setChecked(false);
}

void ArrangementDemoWindow::showInsertMethods(demo_types::TraitsType tabType)
{
  switch(tabType)
  {
  case TraitsType::NONE:
    return;
  case TraitsType::SEGMENT_TRAITS:
    this->ui->actionSegment->setVisible(true);
    this->ui->actionSegment->activate(QAction::Trigger);
    break;
  case TraitsType::POLYLINE_TRAITS:
    this->ui->actionPolyline->setVisible(true);
    this->ui->actionPolyline->activate(QAction::Trigger);
    break;
  case TraitsType::LINEAR_TRAITS:
    this->ui->actionSegment->setVisible(true);
    this->ui->actionSegment->activate(QAction::Trigger);
    this->ui->actionRay->setVisible(true);
    this->ui->actionLine->setVisible(true);
    break;
#ifdef CGAL_USE_CORE
  case TraitsType::CONIC_TRAITS:
    this->ui->actionSegment->setVisible(true);
    this->ui->actionCircle->setVisible(true);
    this->ui->actionEllipse->setVisible(true);
    this->ui->actionConicThreePoint->setVisible(true);
    this->ui->actionConicFivePoint->setVisible(true);
    break;
  case TraitsType::ALGEBRAIC_TRAITS:
    this->ui->actionCircle->setVisible(true);
    this->ui->actionCircle->activate(QAction::Trigger);
    this->ui->actionEllipse->setVisible(true);
    this->ui->actionLine->setVisible(true);
    this->ui->actionAddAlgebraicCurve->setVisible(true);
    break;
  case TraitsType::BEZIER_TRAITS:
    this->ui->actionBezier->setVisible(true);
    this->ui->actionBezier->activate(QAction::Trigger);
    break;
  case TraitsType::RATIONAL_FUNCTION_TRAITS:
    this->ui->actionAddRationalCurve->setVisible(true);
    break;
#endif
  default:
    CGAL_error();
  }
}

void ArrangementDemoWindow::updateInputType(QAction* a)
{
  auto tab = this->getCurrentTab();
  if (!tab) return;

  using namespace CGAL::Qt;
  CurveType curveType = CurveType::None;

  if (a == this->ui->actionSegment)
    curveType = CurveType::Segment;
  else if (a == this->ui->actionPolyline)
    curveType = CurveType::Polyline;
  else if (a == this->ui->actionLine)
    curveType = CurveType::Line;
  else if (a == this->ui->actionRay)
    curveType = CurveType::Ray;
  else if (a == this->ui->actionCircle)
    curveType = CurveType::Circle;
  else if (a == this->ui->actionEllipse)
    curveType = CurveType::Ellipse;
  else if (a == this->ui->actionConicThreePoint)
    curveType = CurveType::ThreePointCircularArc;
  else if (a == this->ui->actionConicFivePoint)
    curveType = CurveType::FivePointConicArc;
  else if (a == this->ui->actionBezier)
    curveType = CurveType::Bezier;

  tab->activateCurveInputCallback(curveType);
}

void ArrangementDemoWindow::on_actionAddAlgebraicCurve_triggered()
{
#ifdef CGAL_USE_CORE
  AlgebraicCurveInputDialog newDialog;

  if (newDialog.exec() == QDialog::Accepted)
  {
    auto currentTab = this->getCurrentTab();

    std::string expression = newDialog.getLineEditText();
    bool is_first_curve = false;

    auto cv = algebraicCurveFromExpression(
      currentTab->getArrangement(), expression,
      is_first_curve);

    if (!cv)
    {
      QMessageBox msgBox;
      msgBox.setWindowTitle("Invalid Expression");
      msgBox.setIcon(QMessageBox::Critical);
      msgBox.setText(QString::fromStdString("Unable to parse expression!"));
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.exec();
      return;
    }

    // adding curve to the arrangement
    auto algCurveInputCallback = currentTab->getCurveInputCallback();
    algCurveInputCallback->generate(cv);
    if (is_first_curve) currentTab->adjustViewport();
  }
#endif
}

void ArrangementDemoWindow::on_actionAddRationalCurve_triggered()
{
#ifdef CGAL_USE_CORE
  RationalCurveInputDialog newDialog;

  if (newDialog.exec() == QDialog::Accepted)
  {
    auto currentTab = this->getCurrentTab();

    std::string num = newDialog.getNumeratorText();
    std::string den = newDialog.getDenominatorText();

    bool is_first_curve = false;

    auto cv = rationalCurveFromExpression(
      currentTab->getArrangement(), num, den, is_first_curve);

    if (!cv)
    {
      QMessageBox msgBox;
      msgBox.setWindowTitle("Invalid Expression");
      msgBox.setIcon(QMessageBox::Critical);
      msgBox.setText(QString::fromStdString("Unable to parse expression!"));
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.exec();
      return;
    }

    // adding curve to the arrangement
    auto algCurveInputCallback = currentTab->getCurveInputCallback();
    algCurveInputCallback->generate(cv);
    if (is_first_curve) currentTab->adjustViewport();
  }
#endif
}

void ArrangementDemoWindow::on_actionQuit_triggered() { qApp->exit(); }

void ArrangementDemoWindow::on_actionNewTab_triggered()
{
  NewTabDialog newTabDialog;
  if (newTabDialog.exec() == QDialog::Accepted)
  {
    int id = newTabDialog.checkedId();
    this->makeTab(static_cast<TraitsType>(id));
  }
}

void ArrangementDemoWindow::on_tabWidget_currentChanged(int)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    auto tt = currentTab->traitsType();
    this->resetCallbackState(currentTab);
    this->resetActionGroups(currentTab, tt);
    this->updateFillColorSwatch(currentTab);
  }
}

void ArrangementDemoWindow::on_actionInsert_toggled(bool checked)
{
  this->hideInsertMethods();

  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    if (checked)
      this->showInsertMethods(currentTab->traitsType());
    else
      this->resetCallbackState(currentTab);
  }
}

void ArrangementDemoWindow::on_actionDrag_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    // TODO: Move this to DemoTab
    QGraphicsView* activeView = currentTab->getView();
    if (!checked)
      activeView->setDragMode(QGraphicsView::NoDrag);
    else
      activeView->setDragMode(QGraphicsView::ScrollHandDrag);
  }
}

void ArrangementDemoWindow::on_actionDelete_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab && !checked)
    this->resetCallbackState(currentTab);
}

void ArrangementDemoWindow::on_actionDelete_triggered()
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
    currentTab->activateDeleteCurveCallback();
}

void ArrangementDemoWindow::on_actionPointLocation_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    if (!checked)
      this->resetCallbackState(currentTab);
    else
      currentTab->activatePointLocationCallback();
  }
}

void ArrangementDemoWindow::on_actionRayShootingUp_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    if (!checked)
      this->resetCallbackState(currentTab);
    else
      currentTab->activateVerticalRayShootCallback(true);
  }
}

void ArrangementDemoWindow::on_actionRayShootingDown_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    if (!checked)
      this->resetCallbackState(currentTab);
    else
      currentTab->activateVerticalRayShootCallback(false);
  }
}

void ArrangementDemoWindow::on_actionMerge_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    if (!checked)
      this->resetCallbackState(currentTab);
    else
      currentTab->activateMergeEdgeCallback();
  }
}

void ArrangementDemoWindow::on_actionSplit_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    if (!checked)
      this->resetCallbackState(currentTab);
    else
      currentTab->activateSplitEdgeCallback();
  }
}

void ArrangementDemoWindow::on_actionFill_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    if (!checked)
      this->resetCallbackState(currentTab);
    else
      currentTab->activateFillFaceCallback();
  }
}


void ArrangementDemoWindow::on_actionShowGrid_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    currentTab->showGrid(checked);
    if (!checked && this->ui->actionGridSnapMode->isChecked())
      this->ui->actionGridSnapMode->activate(QAction::Trigger);
  }
}

void ArrangementDemoWindow::on_actionGridSnapMode_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    currentTab->setSnapToGrid(checked);
    if (checked && !this->ui->actionShowGrid->isChecked())
      this->ui->actionShowGrid->activate(QAction::Trigger);
  }
}

void ArrangementDemoWindow::on_actionArrangementSnapMode_toggled(bool checked)
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    currentTab->setSnapToArrangement(checked);
  }
}

void ArrangementDemoWindow::on_actionCloseTab_triggered()
{
  auto currentTab = this->getCurrentTab();
  if (!currentTab) return;

  // delete the tab
  auto currentTabIndex = this->ui->tabWidget->currentIndex();
  this->ui->tabWidget->removeTab(currentTabIndex);

  // release memory
  delete currentTab;

  // remove the tab
  this->tabs.erase(this->tabs.begin() + currentTabIndex);
}

void ArrangementDemoWindow::on_actionZoomIn_triggered()
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    QGraphicsView* view = currentTab->getView();
    view->scale(2.0, 2.0);
  }
}

void ArrangementDemoWindow::on_actionZoomOut_triggered()
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
  {
    QGraphicsView* view = currentTab->getView();
    view->scale(0.5, 0.5);
  }
}

void ArrangementDemoWindow::on_actionZoomReset_triggered()
{
  auto currentTab = this->getCurrentTab();
  if (currentTab)
    currentTab->adjustViewport();
}

void ArrangementDemoWindow::on_actionFillColor_triggered()
{
  auto currentTab = this->getCurrentTab();
  if (!currentTab) return;

  QColor fillColor = currentTab->getFillFaceColor();

  QColor selectedColor = QColorDialog::getColor(fillColor);
  if (selectedColor.isValid())
  {
    currentTab->setFillFaceColor(selectedColor);
    this->updateFillColorSwatch(currentTab);
  }
}

void ArrangementDemoWindow::updateFillColorSwatch(ArrangementDemoTab* tab)
{
  if (!tab) return;

  QColor fillColor = tab->getFillFaceColor();
  if (!fillColor.isValid()) { fillColor = ::Qt::black; }

  QPixmap fillColorPixmap(16, 16);
  fillColorPixmap.fill(fillColor);
  QIcon fillColorIcon(fillColorPixmap);
  this->ui->actionFillColor->setIcon(fillColorIcon);
}

ArrangementDemoTab* ArrangementDemoWindow::getCurrentTab()
{
  int tabIndex = this->ui->tabWidget->currentIndex();
  if (tabIndex == -1) return nullptr;

  return this->tabs[tabIndex];
}

void ArrangementDemoWindow::on_actionOverlay_triggered()
{
  std::vector<OverlayDialog::ArrangementInfo> arr_infos;
  for (int i = 0; i < this->ui->tabWidget->count(); ++i)
  {
    // this assumes that all tabs are of type ArrangementDemoTab
    ArrangementDemoTab* tab =
      static_cast<ArrangementDemoTab*>(this->ui->tabWidget->widget(i));

    arr_infos.push_back({
        i, tab->traitsType(), this->ui->tabWidget->tabText(i)});
  }

  OverlayDialog overlayDialog{this, arr_infos};
  if (overlayDialog.exec() == QDialog::Accepted)
  {
    std::vector<CGAL::Object> arrs;
    auto selected_arrs = overlayDialog.selectedArrangements();
    if (selected_arrs.empty()) return;

    demo_types::TraitsType overlay_tt =
      static_cast<ArrangementDemoTab*>(
        this->ui->tabWidget->widget(selected_arrs[0]))
        ->traitsType();

    for (auto idx : selected_arrs)
    {
      ArrangementDemoTab* tab =
        static_cast<ArrangementDemoTab*>(this->ui->tabWidget->widget(idx));

      arrs.push_back(tab->getArrangement());
    }

    auto arr_obj = makeOverlayArrangement(arrs);
    if (arr_obj)
    {
      QString tab_label = QString("%1 (Overlay)").arg(makeTabLabel(overlay_tt));
      auto tab = makeTab(overlay_tt, tab_label, arr_obj);

      tab->adjustViewport();
    }
  }
}

void ArrangementDemoWindow::on_actionSaveAs_triggered()
{
  auto currentTab = this->getCurrentTab();
  if (!currentTab)
  {
    QMessageBox::information(this, "Oops", "Create a new tab first");
    return;
  }

  std::pair<CGAL::Object, demo_types::TraitsType> arr_pair = {
    currentTab->getArrangement(), currentTab->traitsType()
  };

  QString filename = QFileDialog::getSaveFileName(
    this, tr("Save file"), "", "Arrangement (*.arr)");
  if (filename.isNull()) return;

  QByteArray ba = filename.toLocal8Bit();
  std::ofstream ofs(ba.data());
  if (!ArrangementIO{}.write(arr_pair, ofs))
    QMessageBox::information(this, "Oops", "Error saving file!");
}

void ArrangementDemoWindow::on_actionOpen_triggered()
{
  QString filename = QFileDialog::getOpenFileName(
    this, tr("Open file"), "", "Arrangement files (*.arr)");
  if (filename.isNull()) return;

  if (filename.endsWith(".arr"))
  {
    auto tab = this->openArrFile(filename);
    if (tab)
    {
      tab->setParent(this);
      this->addTab(tab, this->makeTabLabel(tab->traitsType()));
      tab->adjustViewport();
    }
  }
  else
  {
    QMessageBox::information(this, "Oops", "Unsupported file format");
  }
}

ArrangementDemoTab* ArrangementDemoWindow::openArrFile(QString filename)
{
  if (filename.isNull()) return nullptr;

  QByteArray filename_ba = filename.toLocal8Bit();
  std::ifstream ifs(filename_ba.data());

  ArrangementDemoTab* tab = nullptr;

  std::pair<CGAL::Object, demo_types::TraitsType> arr_pair =
    ArrangementIO{}.read(ifs);
  if (arr_pair.first)
    tab =
      makeTab(arr_pair.second, makeTabLabel(arr_pair.second), arr_pair.first);

  return tab;
}

void ArrangementDemoWindow::on_actionPreferences_triggered()
{
  auto currentTab = this->getCurrentTab();
  if (!currentTab) return;

  ArrangementDemoPropertiesDialog dialog{this};

  if (dialog.exec() == QDialog::Accepted)
  {
    typedef ArrangementDemoPropertiesDialog Dialog;
    ArrangementDemoTab::Preferences pref;

    pref.edgeColor = dialog.property(Dialog::EDGE_COLOR_KEY).value<QColor>();
    pref.edgeWidth =
      dialog.property(Dialog::EDGE_WIDTH_KEY).value<unsigned int>();
    pref.vertexColor =
      dialog.property(Dialog::VERTEX_COLOR_KEY).value<QColor>();
    pref.vertexRadius =
      dialog.property(Dialog::VERTEX_RADIUS_KEY).value<unsigned int>();
    pref.envelopeEdgeColor =
      dialog.property(Dialog::ENVELOPE_EDGE_COLOR_KEY).value<QColor>();
    pref.envelopeEdgeWidth =
      dialog.property(Dialog::ENVELOPE_EDGE_WIDTH_KEY).value<unsigned int>();
    pref.envelopeVertexColor =
      dialog.property(Dialog::ENVELOPE_VERTEX_COLOR_KEY).value<QColor>();
    pref.envelopeVertexRadius =
      dialog.property(Dialog::ENVELOPE_VERTEX_RADIUS_KEY).value<unsigned int>();
    pref.verticalRayEdgeColor =
      dialog.property(Dialog::VERTICAL_RAY_EDGE_COLOR_KEY).value<QColor>();
    pref.verticalRayEdgeWidth =
      dialog.property(Dialog::VERTICAL_RAY_EDGE_WIDTH_KEY)
        .value<unsigned int>();
    pref.axesColor = dialog.property(Dialog::GRID_COLOR_KEY).value<QColor>();
    pref.gridColor = pref.axesColor;
    pref.gridColor.setAlphaF(0.5);
    currentTab->updatePreferences(pref);
  }
}

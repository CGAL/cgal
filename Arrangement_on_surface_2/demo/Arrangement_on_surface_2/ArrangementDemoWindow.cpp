// Copyright (c) 2012  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Alex Tsui <alextsui05@gmail.com>

#include "ArrangementDemoWindow.h"
#include "AlgebraicCurveInputDialog.h"
#include "AlgebraicCurveParser.h"
#include "ArrangementDemoPropertiesDialog.h"
#include "ArrangementDemoTab.h"
#include "Conic_reader.h"
#include "NewTabDialog.h"
#include "OverlayDialog.h"
#include "DeleteCurveCallback.h"
#include "EnvelopeCallback.h"
#include "MergeEdgeCallback.h"
#include "PointLocationCallback.h"
#include "SplitEdgeCallback.h"
#include "VerticalRayShootCallback.h"
#include "ArrangementGraphicsItem.h"
#include "DeleteCurveMode.h"
#include "GraphicsViewCurveInput.h"
#include "FillFaceCallback.h"

#include <QActionGroup>
#include <QColorDialog>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QString>

#include <CGAL/IO/Arr_text_formatter.h>
#include <CGAL/IO/Arr_with_history_iostream.h>
#include <CGAL/IO/Arr_with_history_text_formatter.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Arr_overlay_2.h>

#include "ui_AlgebraicCurveInputDialog.h"

using namespace CGAL::Qt;

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

  this->makeTab(SEGMENT_TRAITS);

  // Call inherited functions
  this->setupStatusBar();
  // this->setupOptionsMenu();
  this->addAboutDemo(":/help/about.html");
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
}

template <class T>
struct TypeHolder { using type = T; };

template <class Lambda>
void ArrangementDemoWindow::visitTraitsType(TraitsType tt, Lambda lambda)
{
  switch (tt)
  {
  default:
  case SEGMENT_TRAITS:
    lambda(TypeHolder<Seg_arr>{});
    break;
  case POLYLINE_TRAITS:
    lambda(TypeHolder<Pol_arr>{});
    break;
#ifdef CGAL_USE_CORE
  case CONIC_TRAITS:
    lambda(TypeHolder<Conic_arr>{});
    break;
#endif
  case LINEAR_TRAITS:
    lambda(TypeHolder<Lin_arr>{});
    break;
  case ALGEBRAIC_TRAITS:
    lambda(TypeHolder<Alg_seg_arr>{});
    break;
  }
}

QString ArrangementDemoWindow::makeTabLabel(TraitsType tt)
{
  static const QString typeNames[] = {
    "Segment", "Polyline", "Conic", "Linear", "Algebraic"};
  return QString("%1 - %2").arg(this->tabLabelCounter++).arg(typeNames[tt]);
}

ArrangementDemoTabBase* ArrangementDemoWindow::makeTab(TraitsType tt)
{
  ArrangementDemoTabBase* demoTab;
  QString tabLabel = makeTabLabel(tt);
  visitTraitsType(tt, [&](auto type_holder) {
    demoTab =
      new ArrangementDemoTab<typename decltype(type_holder)::type>(this);
  });
  this->addTab(demoTab, tabLabel, tt);
  return demoTab;
}

void ArrangementDemoWindow::addTab(
  ArrangementDemoTabBase* demoTab, QString tabLabel, TraitsType tt)
{
  this->tabs.push_back({demoTab, tt});

  this->ui->tabWidget->addTab(demoTab, tabLabel);
  this->ui->tabWidget->setCurrentWidget(demoTab);
  // TODO: This causes memory leak because of the existence of multiple tabs
  // fix it
  this->addNavigation(demoTab->getView());
}

void ArrangementDemoWindow::resetActionGroups()
{
  if (!this->ui->actionInsert->isChecked())
  { this->ui->actionInsert->setChecked(true); }
  else
  {
    this->hideInsertMethods();
    this->showInsertMethods();
  }
}

void ArrangementDemoWindow::resetCallbackState()
{
  auto activeTab = this->getActiveTab().first;
  if (activeTab) activeTab->unhookCallbacks();
}

void ArrangementDemoWindow::updateEnvelope(QAction* newMode)
{
  auto activeTab = this->getActiveTab().first;
  if (!activeTab) return;

  bool show = newMode->isChecked();
  if (newMode == this->ui->actionLowerEnvelope)
  {
    activeTab->getEnvelopeCallback()->showLowerEnvelope(show);
    activeTab->update();
  }
  else if (newMode == this->ui->actionUpperEnvelope)
  {
    activeTab->getEnvelopeCallback()->showUpperEnvelope(show);
    activeTab->update();
  }
}

void ArrangementDemoWindow::hideInsertMethods()
{
  this->ui->actionAddAlgebraicCurve->setVisible(false);
  this->ui->actionAddAlgebraicCurve->setChecked(false);
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
}

void ArrangementDemoWindow::showInsertMethods()
{
  auto tabType = this->getActiveTab().second;
  switch(tabType)
  {
  case NONE:
    return;
  case SEGMENT_TRAITS:
    this->ui->actionSegment->setVisible(true);
    this->ui->actionSegment->activate(QAction::Trigger);
    break;
  case POLYLINE_TRAITS:
    this->ui->actionPolyline->setVisible(true);
    this->ui->actionPolyline->activate(QAction::Trigger);
    break;
  case CONIC_TRAITS:
    this->ui->actionSegment->setVisible(true);
    this->ui->actionCircle->setVisible(true);
    this->ui->actionEllipse->setVisible(true);
    this->ui->actionConicThreePoint->setVisible(true);
    this->ui->actionConicFivePoint->setVisible(true);
    break;
  case LINEAR_TRAITS:
    this->ui->actionSegment->setVisible(true);
    this->ui->actionSegment->activate(QAction::Trigger);
    this->ui->actionRay->setVisible(true);
    this->ui->actionLine->setVisible(true);
    break;
  case ALGEBRAIC_TRAITS:
    this->ui->actionCircle->setVisible(true);
    this->ui->actionCircle->activate(QAction::Trigger);
    this->ui->actionEllipse->setVisible(true);
    this->ui->actionLine->setVisible(true);
    this->ui->actionAddAlgebraicCurve->setVisible(true);
    break;
  case BEZIER_TRAITS:
    break;
  }
}

void ArrangementDemoWindow::updateInputType(QAction* a)
{
  auto tab = this->getActiveTab().first;
  if (!tab) return;

  using namespace CGAL::Qt;
  CurveType curveType;

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

  tab->activateCurveInputCallback(curveType);
}

// should refactor this to GraphicsViewCurveInput if any other use case arises
void ArrangementDemoWindow::on_actionAddAlgebraicCurve_triggered()
{
  AlgebraicCurveInputDialog newDialog;
  newDialog.getUi()->lineEdit->setFocus();

  if (newDialog.exec() == QDialog::Accepted)
  {
    typedef Alg_seg_arr::Geometry_traits_2 Alg_seg_geom_traits;
    typedef typename Alg_seg_geom_traits::Polynomial_2 Polynomial_2;

    std::string algebraicExpression = newDialog.getLineEditText();
    AlgebraicCurveParser<Polynomial_2> parser(algebraicExpression);
    auto poly = parser.parse();

    if (!poly)
    {
      QMessageBox msgBox;
      msgBox.setWindowTitle("Wrong Expression");
      msgBox.setIcon(QMessageBox::Critical);
      msgBox.setText(QString::fromStdString("Invalid Expression"));
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.exec();
      return;
    }

    // To create a curve
    Alg_seg_geom_traits alg_traits;
    auto construct_curve = alg_traits.construct_curve_2_object();

    // adding curve to the arrangement
    auto activeTab = this->getActiveTab().first;
    auto algCurveInputCallback = activeTab->getCurveInputCallback();
    auto cv = construct_curve(poly.value());
    Q_EMIT algCurveInputCallback->generate(CGAL::make_object(cv));
    activeTab->adjustViewport();
  }
}

void ArrangementDemoWindow::on_actionQuit_triggered() { qApp->exit(); }

void ArrangementDemoWindow::on_actionNewTab_triggered()
{
  NewTabDialog newTabDialog;
  if (newTabDialog.exec() == QDialog::Accepted)
  {
    int id = newTabDialog.checkedId();
    // this assumes that options have the same order as TraitsType
    this->makeTab(static_cast<TraitsType>(id));
  }
}

void ArrangementDemoWindow::on_tabWidget_currentChanged(int index)
{
  if (index == -1) return;
  this->resetCallbackState();
  this->resetActionGroups();
  this->updateFillColorSwatch();
}

void ArrangementDemoWindow::on_actionInsert_toggled(bool checked)
{
  this->hideInsertMethods();
  if (checked)
    this->showInsertMethods();
  else
    this->resetCallbackState();
}

void ArrangementDemoWindow::on_actionDrag_toggled(bool checked)
{
  auto activeTab = this->getActiveTab().first;
  if (activeTab)
  {
    // TODO: Move this to DemoTab
    QGraphicsView* activeView = activeTab->getView();
    if (!checked)
      activeView->setDragMode(QGraphicsView::NoDrag);
    else
      activeView->setDragMode(QGraphicsView::ScrollHandDrag);
  }
}

void ArrangementDemoWindow::on_actionDelete_toggled(bool checked)
{
  if (!checked)
    this->resetCallbackState();
}

void ArrangementDemoWindow::on_actionDelete_triggered()
{
  auto activeTab = this->getActiveTab().first;
  if (activeTab)
    activeTab->activateDeleteCurveCallback();
}

void ArrangementDemoWindow::on_actionPointLocation_toggled(bool checked)
{
  if (!checked)
  {
    this->resetCallbackState();
    return;
  }
  auto activeTab = this->getActiveTab().first;
  if (activeTab)
    activeTab->activatePointLocationCallback();
}

void ArrangementDemoWindow::on_actionRayShootingUp_toggled(bool checked)
{
  if (!checked)
  {
    this->resetCallbackState();
    return;
  }
  auto activeTab = this->getActiveTab().first;
  if (activeTab)
    activeTab->activateVerticalRayShootCallback(true);
}

void ArrangementDemoWindow::on_actionRayShootingDown_toggled(bool checked)
{
  if (!checked)
  {
    this->resetCallbackState();
    return;
  }
  auto activeTab = this->getActiveTab().first;
  if (activeTab)
    activeTab->activateVerticalRayShootCallback(false);
}

void ArrangementDemoWindow::on_actionMerge_toggled(bool checked)
{
  if (!checked)
  {
    this->resetCallbackState();
    return;
  }
  auto activeTab = this->getActiveTab().first;
  if (activeTab)
    activeTab->activateMergeEdgeCallback();
}

void ArrangementDemoWindow::on_actionSplit_toggled(bool checked)
{
  if (!checked)
  {
    this->resetCallbackState();
    return;
  }
  auto activeTab = this->getActiveTab().first;
  if (activeTab)
    activeTab->activateSplitEdgeCallback();
}

void ArrangementDemoWindow::on_actionFill_toggled(bool checked)
{
  if (!checked)
  {
    this->resetCallbackState();
    return;
  }
  auto activeTab = this->getActiveTab().first;
  if (activeTab)
    activeTab->activateFillFaceCallback();
}

void ArrangementDemoWindow::on_actionCloseTab_triggered()
{
  auto curTab = this->getActiveTab().first;
  if (!curTab) return;

  // delete the tab
  auto currentTabIndex = this->ui->tabWidget->currentIndex();
  this->ui->tabWidget->removeTab(currentTabIndex);

  // release memory
  delete curTab;

  // remove the tab
  this->tabs.erase(this->tabs.begin() + currentTabIndex);
}

void ArrangementDemoWindow::on_actionZoomIn_triggered()
{
  auto tab = this->getActiveTab().first;
  if (tab)
  {
    QGraphicsView* view = tab->getView();
    view->scale(2.0, 2.0);
  }
}

void ArrangementDemoWindow::on_actionZoomOut_triggered()
{
  auto tab = this->getActiveTab().first;
  if (tab)
  {
    QGraphicsView* view = tab->getView();
    view->scale(0.5, 0.5);
  }
}

void ArrangementDemoWindow::on_actionZoomReset_triggered()
{
  auto tab = this->getActiveTab().first;
  if (tab)
  {
    tab->adjustViewport();
  }
}

void ArrangementDemoWindow::on_actionFillColor_triggered()
{
  auto currentTab = this->getActiveTab().first;
  if (!currentTab) return;

  FillFaceCallbackBase* fillFaceCallback = currentTab->getFillFaceCallback();
  QColor fillColor = fillFaceCallback->getColor();

  QColor selectedColor = QColorDialog::getColor(fillColor);
  if (selectedColor.isValid())
  {
    fillFaceCallback->setColor(selectedColor);
    this->updateFillColorSwatch();
  }
}

void ArrangementDemoWindow::updateFillColorSwatch()
{
  auto currentTab = this->getActiveTab().first;
  if (!currentTab) return;

  QColor fillColor = currentTab->getFillFaceCallback()->getColor();
  if (!fillColor.isValid()) { fillColor = ::Qt::black; }

  QPixmap fillColorPixmap(16, 16);
  fillColorPixmap.fill(fillColor);
  QIcon fillColorIcon(fillColorPixmap);
  this->ui->actionFillColor->setIcon(fillColorIcon);
}

auto ArrangementDemoWindow::getActiveTab()
  -> std::pair<ArrangementDemoTabBase*, TraitsType>
{
  int tabIndex = this->ui->tabWidget->currentIndex();
  if (tabIndex == -1) return {nullptr, NONE};

  return this->tabs[tabIndex];
}

std::vector<QString> ArrangementDemoWindow::getTabLabels() const
{
  std::vector<QString> res;
  for (int i = 0; i < this->ui->tabWidget->count(); ++i)
    res.push_back(this->ui->tabWidget->tabText(i));
  return res;
}

std::vector<CGAL::Object> ArrangementDemoWindow::getArrangements() const
{
  std::vector<CGAL::Object> res;
  for (auto& tab : this->tabs) res.push_back(tab.first->getArrangement());
  return res;
}

// TODO: Clean up all following parts!  (serialization/preferences/overlay)
void ArrangementDemoWindow::on_actionOverlay_triggered()
{
  OverlayDialog overlayDialog{this};
  if (overlayDialog.exec() == QDialog::Accepted)
  {
    std::vector<CGAL::Object> arrs = overlayDialog.selectedArrangements();
    if (arrs.size() == 2)
    {
      Seg_arr* seg_arr;
      Seg_arr* seg_arr2;
      Pol_arr* pol_arr;
      Pol_arr* pol_arr2;

#ifdef CGAL_USE_CORE
      Conic_arr* conic_arr;
      Conic_arr* conic_arr2;
#endif
      Lin_arr* lin_arr;
      Lin_arr* lin_arr2;
      if (CGAL::assign(seg_arr, arrs[0]) && CGAL::assign(seg_arr2, arrs[1]))
        this->makeOverlayTab(seg_arr, seg_arr2, SEGMENT_TRAITS);
      if (CGAL::assign(pol_arr, arrs[0]) && CGAL::assign(pol_arr2, arrs[1]))
        this->makeOverlayTab(pol_arr, pol_arr2, POLYLINE_TRAITS);

#ifdef CGAL_USE_CORE
      if (CGAL::assign(conic_arr, arrs[0]) && CGAL::assign(conic_arr2, arrs[1]))
        this->makeOverlayTab(conic_arr, conic_arr2, CONIC_TRAITS);
#endif
      if (CGAL::assign(lin_arr, arrs[0]) && CGAL::assign(lin_arr2, arrs[1]))
        this->makeOverlayTab(lin_arr, lin_arr2, LINEAR_TRAITS);
    }
  }
}

template <class ArrType>
ArrangementDemoTabBase* ArrangementDemoWindow::makeTab(
  std::unique_ptr<ArrType> arr, QString tabLabel, TraitsType tt)
{
  auto demoTab = new ArrangementDemoTab<ArrType>(this, std::move(arr));
  this->addTab(demoTab, tabLabel, tt);
  return demoTab;
}

template <class ArrType>
void ArrangementDemoWindow::makeOverlayTab(
  ArrType* arr1, ArrType* arr2, TraitsType tt)
{
  QString tabLabel = QString("Overlay Tab");

  auto overlayArr = std::make_unique<ArrType>();
  CGAL::Arr_default_overlay_traits<ArrType> defaultTraits;

  CGAL::overlay(*arr1, *arr2, *overlayArr, defaultTraits);
  makeTab(std::move(overlayArr), tabLabel, tt);
}

void ArrangementDemoWindow::on_actionSaveAs_triggered()
{
  auto tab_tt = this->getActiveTab();
  auto activeTab = tab_tt.first;
  auto tt = tab_tt.second;
  if (!activeTab)
  {
    QMessageBox::information(this, "Oops", "Create a new tab first");
    return;
  }

  QString filename = QFileDialog::getSaveFileName(
    this, tr("Save file"), "", "Arrangement (*.arr)");
  if (filename.isNull()) return;

  QByteArray ba = filename.toLocal8Bit();
  std::ofstream ofs(ba.data());
  CGAL::Object arr = activeTab->getArrangement();
  Seg_arr* seg;
  Pol_arr* pol;
  Lin_arr* lin;
  Alg_seg_arr* alg_seg;

#ifdef CGAL_USE_CORE
  Conic_arr* conic;
#endif

  // write type info
  ofs << static_cast<int>(tt) << std::endl;
  if (CGAL::assign(seg, arr))
  {
    typedef CGAL::Arr_text_formatter<Seg_arr> Seg_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Seg_text_formatter>
      ArrFormatter;
    ArrFormatter arrFormatter;
    CGAL::write(*seg, ofs, arrFormatter);
  }
  else if (CGAL::assign(pol, arr))
  {
    typedef CGAL::Arr_text_formatter<Pol_arr> Pol_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Pol_text_formatter>
      ArrFormatter;
    ArrFormatter arrFormatter;
    CGAL::write(*pol, ofs, arrFormatter);
  }
  else if (CGAL::assign(lin, arr))
  {
    typedef CGAL::Arr_text_formatter<Lin_arr> Lin_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Lin_text_formatter>
      ArrFormatter;
    ArrFormatter arrFormatter;
    CGAL::write(*lin, ofs, arrFormatter);
  }
  else if (CGAL::assign(alg_seg, arr))
  {
    typedef CGAL::Arr_text_formatter<Alg_seg_arr> Arc_text_formatter;
    typedef CGAL::Arr_with_history_text_formatter<Arc_text_formatter>
      ArrFormatter;
    ArrFormatter arrFormatter;
    CGAL::write(*alg_seg, ofs, arrFormatter);
  }
#ifdef CGAL_USE_CORE
  else if (CGAL::assign(conic, arr))
  {
    Conic_reader<Conic_arr::Geometry_traits_2> conicReader;
    conicReader.write_data(ofs, conic->curves_begin(), conic->curves_end());
  }
#endif

  ofs.close();
}

void ArrangementDemoWindow::on_actionOpen_triggered()
{
  QString filename = QFileDialog::getOpenFileName(
    this, tr("Open file"), "", "Arrangement files (*.arr)");
  if (filename.isNull()) return;

  if (filename.endsWith(".arr"))
  {
    auto tab = this->openArrFile(filename);
    tab->adjustViewport();
  }
  else
  {
    QMessageBox::information(this, "Oops", "Unsupported file format");
  }
}

ArrangementDemoTabBase* ArrangementDemoWindow::openArrFile(QString filename)
{
  if (filename.isNull()) return nullptr;

  QByteArray filename_ba = filename.toLocal8Bit();
  std::ifstream ifs(filename_ba.data());

  // read type info
  int tt_int;
  ifs >> tt_int;
  auto tt = static_cast<TraitsType>(tt_int);

  ArrangementDemoTabBase* createdTab = nullptr;
  switch (tt)
  {
#ifdef CGAL_USE_CORE
  case CONIC_TRAITS: {
    Conic_reader<Conic_arr::Geometry_traits_2> conicReader;
    std::vector<Conic_arr::Curve_2> curve_list;
    CGAL::Bbox_2 bbox;
    conicReader.read_data(ifs, std::back_inserter(curve_list), bbox);
    auto conic = std::make_unique<Conic_arr>();
    CGAL::insert(*conic, curve_list.begin(), curve_list.end());
    createdTab = this->makeTab(std::move(conic), this->makeTabLabel(tt), tt);
    break;
  }
#endif
  default:
    visitTraitsType(tt, [&](auto type_holder) {
      using Arr_ = typename decltype(type_holder)::type;
      // Conic has no text_formatter
      using Arr =
        std::conditional_t<std::is_same<Arr_, Conic_arr>::value, Seg_arr, Arr_>;
      using Text_formatter = CGAL::Arr_text_formatter<Arr>;
      using ArrFormatter =
        CGAL::Arr_with_history_text_formatter<Text_formatter>;

      ArrFormatter arrFormatter;
      auto arr = std::make_unique<Arr>();
      CGAL::read(*arr, ifs, arrFormatter);
      createdTab = this->makeTab(std::move(arr), this->makeTabLabel(tt), tt);
    });
  }
  return createdTab;
}

void ArrangementDemoWindow::on_actionPreferences_triggered()
{
  unsigned int currentTabIndex = this->ui->tabWidget->currentIndex();
  if (currentTabIndex == static_cast<unsigned int>(-1)) return;
  ArrangementDemoTabBase* currentTab = this->tabs[currentTabIndex].first;
  CGAL::Qt::ArrangementGraphicsItemBase* agi =
    currentTab->getArrangementGraphicsItem();
  ArrangementDemoGraphicsView* view =
    dynamic_cast<ArrangementDemoGraphicsView*>(currentTab->getView());
  EnvelopeCallbackBase* envelopeCallback = currentTab->getEnvelopeCallback();
  VerticalRayShootCallbackBase* verticalRayShootCallback =
    currentTab->getVerticalRayShootCallback();
  SplitEdgeCallbackBase* splitEdgeCallback = currentTab->getSplitEdgeCallback();

  ArrangementDemoPropertiesDialog* dialog =
    new ArrangementDemoPropertiesDialog(this);
  if (dialog->exec() == QDialog::Accepted)
  {
    typedef ArrangementDemoPropertiesDialog Dialog;

    QColor edgeColor = dialog->property(Dialog::EDGE_COLOR_KEY).value<QColor>();
    unsigned int edgeWidth =
      dialog->property(Dialog::EDGE_WIDTH_KEY).value<unsigned int>();
    QColor vertexColor =
      dialog->property(Dialog::VERTEX_COLOR_KEY).value<QColor>();
    unsigned int vertexRadius =
      dialog->property(Dialog::VERTEX_RADIUS_KEY).value<unsigned int>();
    QColor envelopeEdgeColor =
      dialog->property(Dialog::ENVELOPE_EDGE_COLOR_KEY).value<QColor>();
    unsigned int envelopeEdgeWidth =
      dialog->property(Dialog::ENVELOPE_EDGE_WIDTH_KEY).value<unsigned int>();
    QColor envelopeVertexColor =
      dialog->property(Dialog::ENVELOPE_VERTEX_COLOR_KEY).value<QColor>();
    unsigned int envelopeVertexRadius =
      dialog->property(Dialog::ENVELOPE_VERTEX_RADIUS_KEY)
        .value<unsigned int>();
    QColor verticalRayEdgeColor =
      dialog->property(Dialog::VERTICAL_RAY_EDGE_COLOR_KEY).value<QColor>();
    unsigned int verticalRayEdgeWidth =
      dialog->property(Dialog::VERTICAL_RAY_EDGE_WIDTH_KEY)
        .value<unsigned int>();
    DeleteCurveMode mode =
      dialog->property(Dialog::DELETE_CURVE_MODE_KEY).value<DeleteCurveMode>();
    unsigned int gridSize =
      dialog->property(Dialog::GRID_SIZE_KEY).value<unsigned int>();
    QColor gridColor = dialog->property(Dialog::GRID_COLOR_KEY).value<QColor>();
    // end new for Qt5 version !

    QPen edgesPen(QBrush(edgeColor), edgeWidth);
    QPen verticesPen(QBrush(vertexColor), vertexRadius);
    agi->setEdgesPen(edgesPen);
    agi->setVerticesPen(verticesPen);
    agi->modelChanged();
    view->setGridSize(gridSize);
    view->setGridColor(gridColor);
    envelopeCallback->setEnvelopeEdgeColor(envelopeEdgeColor);
    envelopeCallback->setEnvelopeEdgeWidth(envelopeEdgeWidth);
    envelopeCallback->setEnvelopeVertexColor(envelopeVertexColor);
    envelopeCallback->setEnvelopeVertexRadius(envelopeVertexRadius);
    verticalRayShootCallback->setEdgeColor(verticalRayEdgeColor);
    verticalRayShootCallback->setEdgeWidth(verticalRayEdgeWidth);
    splitEdgeCallback->setColor(edgeColor);
  }
}

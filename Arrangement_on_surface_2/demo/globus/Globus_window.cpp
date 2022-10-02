// Copyright (c) 2022 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Efi Fogel <efifogel@gmail.com>

#include <QActionGroup>
#include <QColorDialog>
#include <QFileDialog>
#include <QInputDialog>
#include <QMessageBox>
#include <QString>
#include <QGraphicsView>
#include <QLabel>

#include <CGAL/Qt/GraphicsViewNavigation.h>

#include "Globus_window.h"
#include "ui_Globus_window.h"

//! \brief constructs
Globus_window::Globus_window(QWidget* parent) :
  CGAL::Qt::DemosMainWindow(parent), ui(new Ui::Globus_window) {
  this->setup_ui();
  // this->setupStatusBar();
  // this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/Arrangement_on_surface_2/about.html");
  this->addAboutCGAL();
}

//! \brief destructs
Globus_window::~Globus_window() {}

//! \brief
void Globus_window::setup_ui() {
  this->ui->setupUi(this);

  this->mode_group = new QActionGroup(this);
  this->mode_group->addAction(this->ui->action_drag);
}

//! \brief
// void Globus_window::resetActionGroups(ArrangementDemoTab* tab, TraitsType tt) {
//   this->hideInsertMethods();
//   this->ui->actionInsert->setChecked(false);
//   this->ui->actionDrag->setChecked(false);
//   this->ui->actionLowerEnvelope->setChecked(false);
//   this->ui->actionUpperEnvelope->setChecked(false);
//   this->ui->actionShowGrid->setChecked(tab->isGridVisible());
//   this->ui->actionGridSnapMode->setChecked(tab->isSnapToGridEnabled());
//   this->ui->actionArrangementSnapMode->setChecked(
//     tab->isSnapToArrangementEnabled());
//   this->ui->actionLowerEnvelope->setChecked(tab->isLowerEnvelopeShown());
//   this->ui->actionUpperEnvelope->setChecked(tab->isUpperEnvelopeShown());

//   if (
//     tt != TraitsType::SEGMENT_TRAITS && tt != TraitsType::POLYLINE_TRAITS &&
//     tt != TraitsType::LINEAR_TRAITS)
//     this->ui->actionArrangementSnapMode->setVisible(false);
//   else
//     this->ui->actionArrangementSnapMode->setVisible(true);

//   // default action group is scrolling
//   this->ui->actionDrag->activate(QAction::Trigger);
// }

//! \brief
// void Globus_window::reset_callback_state(ArrangementDemoTab* tab) {
//   if (tab) {
//     tab->getView()->setDragMode(QGraphicsView::NoDrag);
//     tab->unhookCallbacks();
//   }
// }

//! \brief
void Globus_window::on_action_quit_triggered() { qApp->exit(); }

//! \brief
void Globus_window::on_action_open_triggered() {
  const QString filename =
    QFileDialog::getOpenFileName(this, tr("Open file"), "",
                                 "Globus files (*.gdb)");
  if (filename.isNull()) return;

  if (filename.endsWith(".gdb")) {
    QMessageBox::information(this, "Oops", "Not implemented yet");
    // auto tab = this->openArrFile(filename);
    // if (tab) {
    //   tab->setParent(this);
    //   this->addTab(tab, this->makeTabLabel(tab->traitsType()));
    //   tab->adjustViewport();
    // }
  }
  else {
    QMessageBox::information(this, "Oops", "Unsupported file format");
  }
}

//! \brief
void Globus_window::on_action_preferences_triggered() {
  // ArrangementDemoPropertiesDialog dialog{this};

  // if (dialog.exec() == QDialog::Accepted) {
  //   typedef ArrangementDemoPropertiesDialog Dialog;
  //   ArrangementDemoTab::Preferences pref;

  //   pref.edgeColor = dialog.property(Dialog::EDGE_COLOR_KEY).value<QColor>();
  //   pref.edgeWidth =
  //     dialog.property(Dialog::EDGE_WIDTH_KEY).value<unsigned int>();
  //   pref.vertexColor =
  //     dialog.property(Dialog::VERTEX_COLOR_KEY).value<QColor>();
  //   pref.vertexRadius =
  //     dialog.property(Dialog::VERTEX_RADIUS_KEY).value<unsigned int>();
  //   pref.envelopeEdgeColor =
  //     dialog.property(Dialog::ENVELOPE_EDGE_COLOR_KEY).value<QColor>();
  //   pref.envelopeEdgeWidth =
  //     dialog.property(Dialog::ENVELOPE_EDGE_WIDTH_KEY).value<unsigned int>();
  //   pref.envelopeVertexColor =
  //     dialog.property(Dialog::ENVELOPE_VERTEX_COLOR_KEY).value<QColor>();
  //   pref.envelopeVertexRadius =
  //     dialog.property(Dialog::ENVELOPE_VERTEX_RADIUS_KEY).value<unsigned int>();
  //   pref.verticalRayEdgeColor =
  //     dialog.property(Dialog::VERTICAL_RAY_EDGE_COLOR_KEY).value<QColor>();
  //   pref.verticalRayEdgeWidth =
  //     dialog.property(Dialog::VERTICAL_RAY_EDGE_WIDTH_KEY)
  //       .value<unsigned int>();
  //   pref.axesColor = dialog.property(Dialog::GRID_COLOR_KEY).value<QColor>();
  //   pref.gridColor = pref.axesColor;
  //   pref.gridColor.setAlphaF(0.5);
  //   currentTab->updatePreferences(pref);
  // }
}

//! \brief
void Globus_window::on_action_save_as_triggered() {
    QMessageBox::information(this, "Save as", "Not implemented yet");
}

//! \brief
void Globus_window::on_action_zoom_in_triggered() {
  // QGraphicsView* view = currentTab->getView();
  // view->scale(2.0, 2.0);
}

//! \brief
void Globus_window::on_action_zoom_out_triggered() {
  // QGraphicsView* view = currentTab->getView();
  // view->scale(0.5, 0.5);
}

//! \brief
void Globus_window::on_action_zoom_reset_triggered() {
  // currentTab->adjustViewport();
}

void Globus_window::on_action_drag_toggled(bool checked) {
  // TODO: Move this to DemoTab
  // QGraphicsView* activeView = currentTab->getView();
  // if (!checked)
  //   activeView->setDragMode(QGraphicsView::NoDrag);
  // else
  //   activeView->setDragMode(QGraphicsView::ScrollHandDrag);
}

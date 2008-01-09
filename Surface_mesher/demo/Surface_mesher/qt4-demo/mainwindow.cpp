#include "mainwindow.h"
#include "ui_optionsdialog.h"

#include <QFileDialog>
#include <QUrl>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QtDebug>
#include <QSpinBox>

#include "get_polyhedral_surface.h"
#include <QGLViewer/vec.h>

#include <algorithm> // std::max
#include <cmath> // std::sqrt
#include <boost/format.hpp>


MainWindow::MainWindow() : 
  QMainWindow(),
  sharp_edges_angle_lower_bound(60.),
  sharp_edges_angle_upper_bound(180.)
{
  setupUi(this);
  setAcceptDrops(true);

  viewer_ptr = qFindChild<QGLViewer*>(this, "viewer");

  viewer_ptr->restoreStateFromFile();
  
  surface = get_polyhedral_surface(this,
				   sharp_edges_angle_lower_bound,
				   sharp_edges_angle_upper_bound);
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent *event)
{
  QString filename = event->mimeData()->urls().at(0).path();
  surface_open(filename);
  event->acceptProposedAction();
}

void MainWindow::surface_open(const QString& filename)
{
  surface->open(filename);
}


void MainWindow::on_action_Open_triggered()
{
  QString filename = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                  "",
                                                  tr("OFF Files (*.off)"));
  if(!filename.isNull())
    surface_open(filename);
}

void MainWindow::on_action_Quit_triggered()
{
  qApp->exit();
}

void MainWindow::on_action_Options_triggered()
{
  QDialog *options_dialog = new QDialog(this);
  Ui::OptionDialog ui;
  ui.setupUi(options_dialog);

  QDoubleSpinBox* sb_upper = options_dialog->findChild<QDoubleSpinBox*>("angle_upper_bound");
  QDoubleSpinBox* sb_lower = options_dialog->findChild<QDoubleSpinBox*>("angle_lower_bound");

  if(!sb_lower || !sb_upper) 
    return;

  sb_lower->setValue(sharp_edges_angle_lower_bound);
  sb_upper->setValue(sharp_edges_angle_upper_bound);
  if(options_dialog->exec() == QDialog::Accepted)
  {
    sharp_edges_angle_upper_bound = sb_upper->value();
    sharp_edges_angle_lower_bound = sb_lower->value();
    emit new_sharp_edges_angle_bounds(sharp_edges_angle_lower_bound,
				      sharp_edges_angle_upper_bound);
  }
}

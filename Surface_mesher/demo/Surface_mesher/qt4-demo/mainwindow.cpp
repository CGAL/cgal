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

  connect(viewer_ptr, SIGNAL(drawNeeded()), surface, SLOT(draw()));
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
  qWarning() << QString("Opening file \"%1\"...\n").arg(filename);
  QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));
  surface->open(filename);
  QApplication::restoreOverrideCursor();
  float xmin, ymin, zmin, xmax, ymax, zmax;
  surface->get_bbox(xmin, ymin, zmin, xmax, ymax, zmax);
  const float xcenter = (xmin + xmax) / 2;
  const float ycenter = (ymin + ymax) / 2;
  const float zcenter = (zmin + zmax) / 2;
  const float xdelta = (-xmin + xmax);
  const float ydelta = (-ymin + ymax);
  const float zdelta = (-zmin + zmax);
  const float radius = std::max(std::max(xdelta, ydelta), zdelta) * std::sqrt(3.)/ 2.;
  std::cerr << boost::format("Bounding box: xmin=%1%, ymin=%2%, zmin=%3%\n"
                             "              xmax=%4%, ymax=%5%, zmax=%6%\n"
                             "              center=(%7%, %8%, %9%)")
    % xmin % ymin % zmin % xmax % ymax % zmax
    % xcenter % ycenter % zcenter;
  viewer_ptr->camera()->setSceneCenter(qglviewer::Vec(xcenter, ycenter, zcenter));
  viewer_ptr->camera()->setSceneRadius(radius);                                
  viewer_ptr->setBackgroundColor(Qt::white);
  viewer_ptr->showEntireScene();
  
  QAction* actionInverse_normals = qFindChild<QAction*>(this, "actionInverse_normals");
  if(actionInverse_normals) actionInverse_normals->setChecked(false);
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

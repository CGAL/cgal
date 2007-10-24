#ifndef _MAINWINDOW_H
#define _MAINWINDOW_H

#include "ui_mainwindow.h"

#include <QMainWindow>
#include <QDragEnterEvent>
#include <QDropEvent>

#include "get_polyhedral_surface.h"
#include <QGLViewer/qglviewer.h>
#include <QGLViewer/vec.h>

#include <QtDebug>
#include <QUrl>
#include <QFileDialog>

#include <algorithm> // std::max
#include <cmath> // std::sqrt
#include <boost/format.hpp>

class MainWindow : public  QMainWindow, private Ui::MainWindow 
{
  Q_OBJECT
public:
  MainWindow() :
    surface(get_polyhedral_surface())
  { 
    setupUi(this);
    setAcceptDrops(true);

    viewer = qFindChild<QGLViewer*>(this, "viewer");

    viewer->restoreStateFromFile();

    viewer->setBackgroundColor(Qt::white);

    connect(viewer, SIGNAL(drawNeeded()), surface, SLOT(draw()));
  }

  void dragEnterEvent(QDragEnterEvent *event)
  {
    if (event->mimeData()->hasFormat("text/uri-list"))
      event->acceptProposedAction();
  }

  void dropEvent(QDropEvent *event)
  {
    QString filename = event->mimeData()->urls().at(0).path();
    surface_open(filename);
    event->acceptProposedAction();
  }

  void surface_open(const QString& filename)
  {
    qWarning() << QString("Opening file \"%1\"...\n").arg(filename);
    surface->open(filename);
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
    viewer->camera()->setSceneCenter(qglviewer::Vec(xcenter, ycenter, zcenter));
    viewer->camera()->setSceneRadius(radius);				     
    viewer->setBackgroundColor(Qt::white);
    viewer->showEntireScene();
  }

private slots:
  void on_action_Open_triggered()
  {
    QString filename = QFileDialog::getOpenFileName(this, tr("Open File"),
						    "",
						    tr("OFF Files (*.off)"));
    surface_open(filename);
  }

  void on_action_Quit_triggered()
  {
    qApp->exit();
  }

private:
  Surface* surface;
  QGLViewer* viewer;
};


#endif // _MAINWINDOW_H

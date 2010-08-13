#include "mainwindow.h"

#include <QFileDialog>
#include <QUrl>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QtDebug>
#include <QMenu>
#include <QAction>
#include <QVariant>
#include <QStringList>
#include <QtGlobal>
#include <QToolBar>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QSettings>

#include <QGLViewer/vec.h>

#include <algorithm> // std::max
#include <cmath> // std::sqrt
#include <boost/format.hpp>

#include "ui_mainwindow.h"
#include "volume.h"
#ifndef CGAL_DO_NOT_USE_POLYHEDRAL_SURFACE
#  include "polyhedral_surface.h"
#endif

MainWindow::MainWindow(MainWindow* other_window /* = 0 */) : 
  CGAL::Qt::DemosMainWindow(),
  surface(0),
  sharp_edges_angle_lower_bound(60.),
  sharp_edges_angle_upper_bound(180.)
{
  setupUi(this);
  setAcceptDrops(true);

  if(other_window != 0)
  {
    viewer->setCamera(other_window->viewer->camera());
    connect(other_window, SIGNAL(destroyed()),
            this, SLOT(close()));
  }

  this->addAboutCGAL();
  this->addRecentFiles(this->menu_File,
		       this->action_Quit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(surface_open(QString)));

  this->readState();

  show_only("");
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("text/uri-list"))
    event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent *event)
{
  QString filename = event->mimeData()->urls().at(0).toLocalFile();
  surface_open(filename);
  event->acceptProposedAction();
}

void MainWindow::surface_open(const QString& filename)
{
  if(surface != 0) {
    delete surface;
    surface = 0;
  }
#ifndef CGAL_DO_NOT_USE_POLYHEDRAL_SURFACE
  surface = new Polyhedral_surface(this);
  if(surface->open(filename)) {
    this->addToRecentFiles(filename);
    return;
  }
  delete surface;
  surface = 0;
#endif
  surface = new Volume(this);
  if(surface->open(filename)) {
    this->addToRecentFiles(filename);
  }
}

void MainWindow::show_only(QString tag)
{
#if 0
  QTextStream err(stderr);
#else
  QString dummy;
  QTextStream err(&dummy);
#endif
  err << "** Show only in \"" << tag << "\"\n";
  Q_FOREACH(QObject* object, 
            this->findChildren<QObject*>())
  {
    QStringList show_only_in = object->property("show_only_in").toStringList();
    if(!show_only_in.isEmpty())
    {
      err << object->metaObject()->className()
          << " \"" << object->objectName() << "\" only in: ";
      foreach(QString s, show_only_in)
        err << s << " ";
      const bool visible = show_only_in.contains(tag);
      err << (visible ? "(enabled)\n" : "(disabled)\n");
      if(QMenu* menu = qobject_cast<QMenu*>(object)) {
        menu->menuAction()->setVisible(visible);
      }
      else {
	object->setProperty("visible", QVariant::fromValue<bool>(visible));
      }
    }
  }
}

void MainWindow::on_action_Open_triggered()
{
  QSettings settings;
  QString directory = settings.value("Open directory",
				     QDir::current().dirName()).toString();
  QString filename = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                  directory,
                                                  tr("all Files (*.*)"));
  if(!filename.isEmpty()) {
    QFileInfo fileinfo(filename);
    if(fileinfo.isFile() && fileinfo.isReadable()) {
      settings.setValue("Open directory",
			fileinfo.absoluteDir().absolutePath());
      surface_open(filename);
    }
  }
}

void MainWindow::on_action_Quit_triggered()
{
  this->writeState();
  qApp->exit();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
  this->writeState();
  event->accept();
}

void MainWindow::on_action_Clone_triggered()
{
  MainWindow* other = new MainWindow(this);
  other->show();
}

#include "mainwindow.moc"

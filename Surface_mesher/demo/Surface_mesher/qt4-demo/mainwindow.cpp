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

#include <QGLViewer/vec.h>

#include "ui_meshing_bar.h"

#include <algorithm> // std::max
#include <cmath> // std::sqrt
#include <boost/format.hpp>

#include "ui_mainwindow.h"
#include "volume.h"

MainWindow::MainWindow(MainWindow* other_window /* = 0 */) : 
  CGAL::Qt::DemosMainWindow(),
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

  QToolBar* tb_meshing = qFindChild<QToolBar*>(this, "toolBar_meshing");
//   tb_meshing->setVisible(false);
  QAction* action_mc = qFindChild<QAction*>(this, "actionMarching_cubes");

  if(tb_meshing && action_mc) {
    QWidget* meshing_bar = new QWidget;
    Ui::meshing_bar ui;
    ui.setupUi(meshing_bar);
    tb_meshing->insertWidget(action_mc, 
                             meshing_bar);
    tb_meshing->insertSeparator(action_mc);
  }

  show_only("");
  surface = new Volume(this);

  addAboutCGAL();
  this->addRecentFiles(this->menu_File,
		       this->action_Quit);
  connect(this, SIGNAL(openRecentFile(QString)),
	  this, SLOT(surface_open(QString)));
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
  this->addToRecentFiles(filename);
}

void MainWindow::show_only(QString tag)
{
  QTextStream err(stderr);
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
      object->setProperty("visible", QVariant::fromValue<bool>(visible));
      if(QMenu* menu = qobject_cast<QMenu*>(object))
        menu->menuAction()->setVisible(visible);
    }
  }
}

void MainWindow::on_action_Open_triggered()
{
  QString filename = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                  "",
                                                  tr("all Files (*.*)"));
  if(!filename.isEmpty())
    surface_open(filename);
}

void MainWindow::on_action_Quit_triggered()
{
  qApp->exit();
}

void MainWindow::on_action_Clone_triggered()
{
  MainWindow* other = new MainWindow(this);
  other->show();
}

#include "mainwindow.moc"

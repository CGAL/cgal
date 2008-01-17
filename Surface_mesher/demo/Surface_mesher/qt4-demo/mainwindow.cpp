#include "mainwindow.h"

#include <QFileDialog>
#include <QUrl>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QtDebug>
#include <QMenu>
#include <QAction>
#include <QtGlobal>

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
  fix_menus_visibility();
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

void MainWindow::fix_menus_visibility()
{
  fix_one_menu_visibility(findChild<QMenu*>("menuEdit"));
  fix_one_menu_visibility(findChild<QMenu*>("menuOptions"));
}

void MainWindow::fix_one_menu_visibility(QMenu* menu)
{
  if(menu) {
    bool is_non_empty = false;
    Q_FOREACH(QAction* action, menu->actions()) {
      is_non_empty = is_non_empty || action->isVisible();
      if(action->isVisible())
        QTextStream(stderr) << action->text() << "\n";
    }
    menu->setEnabled(is_non_empty);
  }
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

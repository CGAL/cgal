#include "mainwindow.h"

#include <QFileDialog>
#include <QUrl>
#include <QDragEnterEvent>
#include <QDropEvent>
#include <QtDebug>
#include <QMenu>
#include <QAction>
#include <QtGlobal>
#include <QToolBar>
#include <QDoubleSpinBox>
#include <QLabel>

#include "get_polyhedral_surface.h"
#include <QGLViewer/vec.h>

#include "ui_meshing_bar.h"

#include <algorithm> // std::max
#include <cmath> // std::sqrt
#include <boost/format.hpp>

#include "volume.h"

MainWindow::MainWindow(MainWindow* other_window /* = 0 */) : 
  QMainWindow(),
  sharp_edges_angle_lower_bound(60.),
  sharp_edges_angle_upper_bound(180.)
{
  setupUi(this);
  setAcceptDrops(true);

  viewer_ptr = qFindChild<QGLViewer*>(this, "viewer");
  Q_ASSERT_X(viewer_ptr, "MainWindow constructor", "cannot find widget \"viewer\"");

  if(other_window != 0)
  {
    QGLViewer* other_viewer_ptr = qFindChild<QGLViewer*>(other_window, "viewer");
    viewer_ptr->setCamera(other_viewer_ptr->camera());
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

  surface = new Volume(this);
  
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
//       if(action->isVisible())
//         QTextStream(stderr) << action->text() << "\n";
    }
    menu->setVisible(is_non_empty);
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

void MainWindow::on_action_Clone_triggered()
{
  MainWindow* other = new MainWindow(this);
  other->show();
}

#include "mainwindow.moc"

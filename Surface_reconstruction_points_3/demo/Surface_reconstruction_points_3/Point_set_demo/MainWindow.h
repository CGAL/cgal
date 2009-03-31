#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include "config.h"

#include <QtOpenGL/qgl.h>
#include <CGAL/Qt/DemosMainWindow.h>

#include <QVector>

class QDragEnterEvent;
class QDropEvent;
class Scene;
class Viewer;
class QTreeView;
class QMenu;
class Point_set_demo_io_plugin_interface;

namespace Ui {
  class MainWindow;
}

#include "Polyhedron_type_fwd.h"

#include "Messages_interface.h"

class MainWindow : 
  public CGAL::Qt::DemosMainWindow,
  public Messages_interface
{
  Q_OBJECT
  Q_INTERFACES(Messages_interface)
public:
  MainWindow(QWidget* parent = 0);
  ~MainWindow();

public slots:
  void updateViewerBBox();
  void open(QString filename);

  void selectPolyhedron(int i);

  void setAddKeyFrameKeyboardModifiers(Qt::KeyboardModifiers);

  void clearMenu(QMenu*);
  void enableAction(QAction*);

  void information(QString);
  void warning(QString);
  void error(QString);

protected slots:
  void selectionChanged();
  void updateInfo();

  // settings
  void quit();
  void readSettings();
  void writeSettings();

  // load, erase, duplicate
  void on_actionEraseAll_triggered();
  void on_actionLoadPolyhedron_triggered();
  bool on_actionErasePolyhedron_triggered();
  void on_actionDuplicatePolyhedron_triggered();

  // Show/Hide
  void on_actionShowHide_triggered();

  // save selected polyhedron as...
  void on_actionSaveAs_triggered(); 

protected:
  void message(QString, QString, QString = QString("normal"));
  void loadPlugins();
  bool initPlugin(QObject*);
  bool initIOPlugin(QObject*);

  void dragEnterEvent(QDragEnterEvent *event);
  void dropEvent(QDropEvent *event);
  void closeEvent(QCloseEvent *event);

  bool onePolygonIsSelected() const;
  int getSelectedPolygonIndex() const;

  Polyhedron* getSelectedPolygon();
private:
  QString strippedName(const QString &fullFileName);

  Scene* scene;
  Viewer* viewer;
  QTreeView* treeView;
  Ui::MainWindow* ui;
  QVector<Point_set_demo_io_plugin_interface*> io_plugins;
};

#endif // ifndef MAINWINDOW_H

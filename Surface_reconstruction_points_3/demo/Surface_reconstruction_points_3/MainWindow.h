#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include "config.h"

#include <QtOpenGL/qgl.h>
#include <CGAL/Qt/DemosMainWindow.h>

#include <QVector>

class Scene;
class Viewer;
class QTreeView;
class QMenu;
class Polyhedron_demo_io_plugin_interface;

class Scene_item;

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

  void selectSceneItem(int i);

  void setAddKeyFrameKeyboardModifiers(Qt::KeyboardModifiers);

  void clearMenu(QMenu*);
  void addAction(QAction*);

  void information(QString);
  void warning(QString);
  void error(QString);

protected slots:
  void selectionChanged();
  void updateInfo();
  void removeManipulatedFrame(Scene_item*);

  // settings
  void quit();
  void readSettings();
  void writeSettings();

  // load, erase, duplicate
  void on_actionFileCloseAll_triggered();
  void on_actionFileOpen_triggered();
  bool on_actionFileClose_triggered();
  void on_actionDuplicate_triggered();
  void on_actionConvertToPointSet_triggered();

  // selection
  void on_actionDeleteSelection_triggered(); 
  void on_actionResetSelection_triggered(); 

  // Show/Hide
  void on_actionShowHide_triggered();

  // save as...
  void on_actionSaveAs_triggered(); 

protected:
  void message(QString, QString, QString = QString("normal"));
  void loadPlugins();
  bool initPlugin(QObject*);
  bool initIOPlugin(QObject*);

  void closeEvent(QCloseEvent *event);

  int getSelectedSceneItemIndex() const;

private:
  QString strippedName(const QString &fullFileName);

  Scene* scene;
  Viewer* viewer;
  QTreeView* treeView;
  Ui::MainWindow* ui;
  QVector<Polyhedron_demo_io_plugin_interface*> io_plugins;
};

#endif // ifndef MAINWINDOW_H

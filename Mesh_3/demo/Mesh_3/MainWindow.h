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
class Io_plugin_interface;

class Scene_item;

namespace Ui {
  class MainWindow;
}

#include <CGAL_demo/Messages_interface.h>

class MainWindow : 
  public CGAL::Qt::DemosMainWindow,
  public Messages_interface
{
  Q_OBJECT
  Q_INTERFACES(Messages_interface)
public:
  MainWindow(QWidget* parent = 0);
  ~MainWindow();

public Q_SLOTS:
  void updateViewerBBox();
  void open(QString filename);

  void selectSceneItem(int i);
  void selectSceneItem();

  void setAddKeyFrameKeyboardModifiers(Qt::KeyboardModifiers);

  void clearMenu(QMenu*);
  void addAction(QAction*);

  void information(QString);
  void warning(QString);
  void error(QString);

protected Q_SLOTS:
  void selectionChanged();
  void updateInfo();
  void updateDisplayInfo();
  void removeManipulatedFrame(Scene_item*);

  // settings
  void quit();
  void readSettings();
  void writeSettings();

	// snapshot
	void on_actionCopy_snapshot_triggered();
	void on_actionSave_snapshot_triggered();

  // load, erase, duplicate
  void on_actionEraseAll_triggered();
  void on_actionLoad_triggered();
  bool on_actionErase_triggered();
  void on_actionDuplicate_triggered();

  // Show/Hide
  void on_actionShowHide_triggered();

  // Select A/B
  void on_actionSetPolyhedronA_triggered();
  void on_actionSetPolyhedronB_triggered();

  // save as...
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
  int getSelectedSceneItemIndex() const;

private:
  QString strippedName(const QString &fullFileName);

  Scene* scene;
  Viewer* viewer;
  QTreeView* treeView;
  Ui::MainWindow* ui;
  QVector<Io_plugin_interface*> io_plugins;
};

#endif // ifndef MAINWINDOW_H

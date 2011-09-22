#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include "config.h"

#include <QtOpenGL/qgl.h>
#include <CGAL/Qt/DemosMainWindow.h>
#ifdef QT_SCRIPT_LIB
#  include  <QScriptEngine>
#endif

#include <QVector>
#include <QList>
#include <QFileInfo>
#include <QStringList>

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
  void open(QString filename, bool no_popup = false);
  Scene_item* load_item(QFileInfo) const;
  void reload_item();

  void selectSceneItem(int i);
  void showSelectedPoint(double, double, double);
  void unSelectSceneItem(int i);
  void selectAll();
  void addSceneItemInSelection(int i);
  void removeSceneItemFromSelection(int i); // same as unSelectSceneItem

  void setAddKeyFrameKeyboardModifiers(Qt::KeyboardModifiers);

  void clearMenu(QMenu*);
  void addAction(QAction*);
  void addAction(QString actionName,
                 QString actionText,
                 QString menuName);
  void viewerShow(float, float, float);

  void information(QString);
  void warning(QString);
  void error(QString);
  void message(QString, QString, QString = QString("normal"));

  bool hasPlugin(QString);
  void enableScriptDebugger(bool = true);

protected slots:
  void selectionChanged();

  void contextMenuRequested(const QPoint& global_pos);
  void showSceneContextMenu(int selectedItemIndex,
                            const QPoint& global_pos);
  void showSceneContextMenu(const QPoint& local_pos_of_treeview);

  void updateInfo();
  void updateDisplayInfo();
  void removeManipulatedFrame(Scene_item*);

  // settings
  void quit();
  void readSettings();
  void writeSettings();

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
  void save(QString filename, Scene_item* item);

  void on_actionSetBackgroundColor_triggered();

  void on_action_Look_at_triggered();

  QString camera_string() const;
  void on_actionDumpCamera_triggered();
  void on_action_Copy_camera_triggered();
  void on_action_Paste_camera_triggered();

protected:
  void loadPlugins();
  bool initPlugin(QObject*);
  bool initIOPlugin(QObject*);

  void closeEvent(QCloseEvent *event);

  bool onePolygonIsSelected() const;
  int getSelectedSceneItemIndex() const;
  QList<int> getSelectedSceneItemIndices() const;

private:
  QString strippedName(const QString &fullFileName);

  Scene* scene;
  Viewer* viewer;
  QTreeView* treeView;
  Ui::MainWindow* ui;
  QVector<Polyhedron_demo_io_plugin_interface*> io_plugins;
  QStringList plugins;
#ifdef QT_SCRIPT_LIB
  QScriptEngine* script_engine;
public:
  void evaluate_script(QString script, 
                       const QString & fileName = QString(),
                       const bool quiet = false);
  void evaluate_script_quiet(QString script, 
                             const QString & fileName = QString());
#endif
};

#endif // ifndef MAINWINDOW_H

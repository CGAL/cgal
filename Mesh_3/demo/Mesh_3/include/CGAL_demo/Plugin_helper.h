#ifndef POLYHEDRON_DEMO_OPERATION_HELPER_H
#define POLYHEDRON_DEMO_OPERATION_HELPER_H

#include <CGAL_demo/Scene_item_config.h> //defines SCENE_ITEM_EXPORT

#include <QString>
#include <QStringList>
#include <QMap>

class QAction;
struct QMetaObject;
class QMainWindow;
class Scene_interface;

#include <CGAL_demo/Plugin_interface.h>

class SCENE_ITEM_EXPORT Plugin_helper
  : public Plugin_interface
{
  typedef Plugin_interface Base;
  
public:
  // get action object from its name
  static QAction* getActionFromMainWindow(QMainWindow*, QString action_name);
  
  // Init plugin
  using Base::init;
  virtual void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  
  // Get list of actions supported by this plugin
  virtual QStringList actionsNames() const;
  virtual QList<QAction*> actions() const;

  // Auto-connect actions to slots. Called by init().
  void autoConnectActions();
  
protected:
  QMap<QString, QAction*> actions_map;
  Scene_interface* scene;
  QMainWindow* mw;
};

#endif // POLYHEDRON_DEMO_OPERATION_HELPER_H
